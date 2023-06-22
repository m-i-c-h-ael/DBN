rm(list=ls())
cat("\014")  #is the code to send ctrl+L to the console and therefore will clear the screen
plot.new()

# to add
#allow diff. timesteps
#dephosphorylation
#estimation of SD
#replicates
#no negative concentration
#relative phosphorylation -> max=1
#JAGS
#noise
#fixed vs. relative noise

library('deSolve')
library('animation')
library('rjags')
library('tidyverse')
library('plotly')
source('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/DBDA2E-utilities.R')
# `$.error` <- function(x, name) {   #not working
#   if (!name %in% colnames(x)) {
#     stop(paste("Column", name, "does not exist in the dataframe."))
#   }
#   x[[name]]
# }

# True:
 #A is phosphorylated by X, which is present in constant concentration
 #B is phosphorylated by Ap
 #C is phosphorylated by Bp
 #all are dephosphorylated at constant rate
 #all start from unphosphorylated state
 #assuming mass action kinetics
 #what do data need to look like for parameters to be identifiable; not too fast phosphorylation; S/N?

# rates
k_XA= .2 #phosphorylation of A by X
k_AB= .3
k_BC= 1
k_WZ= 1
k_XW= 0
k_WU= 1.5
k_WV= 1.5

k_dA= 1 #dephosphorylation
k_dB= 1
k_dC= 1
k_dZ= 0
k_dW= 0
k_dU= .05
k_dV= 0

# initial values
Y0= c(Ap=0,Bp=0,Zp=0.01,Wp=0.01,Up=0.1,Vp=0.1)

parmsX= matrix(c(k_XA,k_AB,k_BC,k_WZ,k_XW, k_WU, k_WV,
                 k_dA,k_dB,k_dC,k_dZ,k_dW, k_dU, k_dV),ncol=2,byrow=FALSE)  #col1= phosphorylation; col2: dephosphorylation
dYp_dt= function(t,Y,parms) {    #Y= Ap, Bp; parmsX: col1: k_XA, k_AB, k_BC; col2: k_dA, k_dB, k_dC
  dAp_dt= parms[1,1]*(1-Y[1]) - parms[1,2]*Y[1]  #saturated 1-order phosphorylation (w/o explicit kinase); 
                        # 0-order dephosphorylation; assuming Atot = 1
  dBp_dt= Y[1]*parms[2,1]*(1-Y[2]) - parms[2,2]*Y[2]
  dZp_dt= parms[4,1]*Y[4]  #0-order phosphorylation rate
  dWp_dt= 0 #constant kinase
  dUp_dt= parms[6,1]*Y[4] - parms[6,2]*Y[5] #0-order phosphorylation and first-order dephosphorylation
  dVp_dt= Y[4]*parms[7,1]*(1-Y[6])
  return(list(c(dAp_dt, dBp_dt, dZp_dt, dWp_dt, dUp_dt, dVp_dt))) 
}

par(mfrow=c(1,2))
concs_out= data.frame(ode(func= dYp_dt,y=Y0,times=0:1000,parms=parmsX))
plot(concs_out$t,concs_out$Ap,ylim=c(0,1),main='A')
plot(concs_out$t,concs_out$Bp,ylim=c(0,1),main='B')
plot(concs_out$t,concs_out$Zp,ylim=c(0,1),main='Z')
lines(concs_out$t,concs_out$Wp,col='blue',lty=2)
plot(concs_out$t,concs_out$Up,ylim=c(0,1),main='U')
lines(concs_out$t,concs_out$Wp,col='blue',lty=2)
plot(concs_out$t,concs_out$Vp,ylim=c(0,1),main='V')
#concs_out -> columns time, concA, concB




##############################
#model k_XZ (NO DEPHOSPHORYLATION; 0-order phosphorylation w/o saturation) - using JAGS
#Prior: assume gamma-distribution with low mean and low effective sample size for k_X

#prepare model
modelString2= "
 model {
  for (i in 2:Ntotal) {
    y[i] ~ dnorm(m[i], 1/sigma^2)  #likelihood
    m[i]= y[i-1]+ k1dt
  } 
  k1dt ~ dnorm(0, 1/0.05^2) #prior
  sigma ~ dgamma(1,2)  #mean: a/b; sd= sqrt(a)/b
 }
"

initList= function()  { #initial values generated via function: mu centered around 0; sigma in range 0-0.1
  set.seed(123)
  initMu= rnorm(1,0,0.1)
  initSigma= runif(1,0,0.1)
  return( list( k1dt= initMu, sigma=initSigma))
}


#relative noise
# real_relSD= 0.01
# Zp= concs_out$Zp
# Zp_wNoise= Zp+ sapply(concs_out$Zp, function(x){rnorm(n=1,mean=0,sd= x*real_relSD)})
# plot(concs_out$t,Zp_wNoise)

sigma_real_vec= 0.01 # c(0.01, 0.1,0.2,0.5,1)  
n_rounds= 1
postMean_mtx= matrix(NA,nrow=n_rounds*length(sigma_real_vec),ncol=2)
colnames(postMean_mtx)= c('sigma_real','postMean_k1dt')
for (sr in 1:length(sigma_real_vec)) {
  sigma_real= sigma_real_vec[sr]
  for (r in 1:n_rounds) {
      #absolute noise
    set.seed(234)
    Zp_wNoise= concs_out$Zp + rnorm(length(concs_out$Zp),mean=0,sd= sigma_real)
    Zp_wNoise= ifelse(Zp_wNoise<0,0,Zp_wNoise)
    print(plot(concs_out$t,Zp_wNoise,main=paste('noise:',sigma_real)))
    
    #prepare data
    #deltat= diff(concs_out$t)
    
    Ntotal= length(Zp_wNoise)
    dataList2= list(
      y= Zp_wNoise,
      Ntotal= Ntotal
    )
    
    jagsModel= jags.model( file= textConnection(modelString2), data=dataList2, inits=initList,
                           n.chains=1, n.adapt=500)  #find MCMC sampler
    #update( jagsModel,n.iter=500)  #burn in
    
    codaSamples= coda.samples(jagsModel, variable.names=c('k1dt','sigma'),n.iter=10000)
    
    # par(mfrow= c(1,1))
    # plotPost(codaSamples[,'k1dt'],main='k1dt',xlab= 'k1 * delta_t')
    #   #true parameter should be k_WZ * Wp= 0.01
    # 
    # diagMCMC(codaObject = codaSamples, parName = 'k1dt')
    # dev.off()
    # 
    # plotPost(codaSamples[,'sigma'],main='sigma',xlab= 'sigma')
    # diagMCMC(codaObject = codaSamples, parName = 'sigma')
    # dev.off()
    # 
    # plot(type='n',x=0,y=0,xlim=c(1,dim(codaSamples[[1]])[1]),xlab='iteration',ylab='k1 * deltat' )
    # cols= c('black','blue','green')
    # for (i in 1:length(codaSamples)){
    #   lines(codaSamples[[i]][,1], col= cols[i])
    # }
    # plot(type='n',x=0,y=0,xlim=c(1,dim(codaSamples[[1]])[1]),xlab='iteration',ylab='sigma' )
    # cols= c('black','blue','green')
    # for (i in 1:length(codaSamples)){
    #   lines(codaSamples[[i]][,2], col= cols[i])
    # }
    # plot(type='n',x=0,y=0,xlab='k1 * deltat', ylab='sigma' )
    # for (i in 1:length(codaSamples)){
    #   points(codaSamples[[i]][,1],codaSamples[[i]][,2])
    # }
    
    # plot(x=as.numeric(codaSamples[[1]][,1]), y= as.numeric(codaSamples[[1]][,2]),type='l',
    #      xlab='mu',ylab='sigma')
    
    #plot(codaSamples)

    r_idx= (sr-1)*n_rounds+r
    postMean_mtx[r_idx,]= c(sigma_real,mean(codaSamples[[1]][,1]))
    
    #based on the simulated parameters, simulate traces
    sim= data.frame(codaSamples[[1]])
    head(sim)
    mu_sim= matrix(NA,nrow= dim(codaSamples[[1]])[1],ncol=dim(concs_out)[1])
    mu_sim[,1]= rep(Y0[3],dim(mu_sim)[1]) #starting value
    for(i in 2:dim(mu_sim)[2]) { #one TP after another
      mu_sim[,i]= mu_sim[,i-1]+ sim$k1dt
    }
    par(mfrow=c(1,1))
    plot(x=concs_out$time,y=concs_out$Zp,main=paste('Real sigma:',sigma_real_vec[sr]))
    for(i in 9000:9100){   #
      lines(concs_out$time,mu_sim[i,],col='red')
    }
    
    kde= kde2d(sim$k1dt,sim$sigma)
    contour(kde,xlab= 'k1dt',ylab='sigma',xlim=c(0,max(sim$k1dt)),ylim=c(0,max(c(sim$k1dt,sim$sigma))))
    
    realSteps_Zp= diff(Zp_wNoise)
    plot(concs_out$time[2:length(concs_out$time)],realSteps_Zp,main="Difference between steps",
          xlab='time')
    abline(h=mean(sim$k1dt),col='red',lwd=3)
    abline(h=quantile(sim$k1dt,.1),col='red',lwd=1,lty=2)
    abline(h=quantile(sim$k1dt,.9),col='red',lwd=1,lty=2)
    }
}

postMean_stats= data.frame(postMean_mtx) %>%
  group_by(sigma_real) %>%
  summarize(m= mean(postMean_k1dt), q10=quantile(postMean_k1dt,0.1), q90=quantile(postMean_k1dt,0.9))

ggplot(data= postMean_stats)+
  geom_point(aes(x=sigma_real,y=m))+
  geom_errorbar(aes(x=sigma_real,ymin=q10,ymax=q90))+
  geom_hline(yintercept= k_WZ*Y0['Wp'],col='blue',lty=2)+
  coord_cartesian(xlim=c(0,1),ylim=c(0,0.015))+
  scale_y_continuous(name= 'mean (Posterior k1dt)')


##############################
#model (NO DEPHOSPHORYLATION; first order phosphorylation with saturation) - using JAGS

#specify model
modelString3= "
 model {
  for (i in 2:length(y)) {
    y[i] ~ dnorm(m[i], 1/sigma^2)  #likelihood
    m[i]= y[i-1]- k1dt*w[i-1]*y[i-1]+ k1dt*w[i-1]*y_tot
  } 
  k1dt ~ dnorm(0, 1/0.05^2) #prior
  sigma ~ dgamma(1,2)  #mean: a/b; sd= sqrt(a)/b
  y_tot ~ dgamma(1,2)  #to be optimized; should be positive
 }
"

initList3= function()  { #initial values generated via function: mu centered around 0; sigma in range 0-0.1
  set.seed(123)
  initMu= rnorm(1,0,0.1)
  initSigma= runif(1,0,0.1)
  init_ytot= 1
  return( list( k1dt= initMu, sigma=initSigma, y_tot= init_ytot))
}

#relative noise
real_relSD= 0.01
Wp= concs_out$Wp
Vp= concs_out$Vp
Vp_wNoise= Vp+ sapply(concs_out$Vp, function(x){rnorm(n=1,mean=0,sd= x*real_relSD)})
plot(concs_out$t,Vp_wNoise)

dataList3= list(
  y= Vp_wNoise,
  w= concs_out$Wp
)

jagsModel3= jags.model( file= textConnection(modelString3), data=dataList3, inits=initList3,
                       n.chains=2, n.adapt=500)  #find MCMC sampler

codaSamples3= coda.samples(jagsModel3, variable.names=c('k1dt','sigma','y_tot'),n.iter=10000)
plot(codaSamples3)

diagMCMC(codaObject = codaSamples3, parName = 'k1dt')
diagMCMC(codaObject = codaSamples3, parName = 'y_tot')
 dev.off()

#based on the simulated parameters, simulate traces
sim3= data.frame(codaSamples3[[1]])
head(sim3)
mu_sim= matrix(NA,nrow= dim(codaSamples3[[1]])[1],ncol=dim(concs_out)[1])
mu_sim[,1]= rep(Y0[6],dim(mu_sim)[1]) #starting value
for(i in 2:dim(mu_sim)[2]) { #one TP after another
  mu_sim[,i]= mu_sim[,i-1]- sim3$k1dt*concs_out$Wp[i-1]*mu_sim[,i-1] + 
    sim3$k1dt*concs_out$Wp[i-1]*sim3$y_tot
}
par(mfrow=c(1,1))
plot(concs_out$time,concs_out$Vp)
for(i in 9000:9100){   #
  lines(concs_out$time,mu_sim[i,],col='red')
}


kde3= kde2d(sim3$k1dt,sim3$y_tot)
contour(kde3,xlab= 'k1dt',ylab='y_tot',xlim=c(min(c(0,sim3$k1dt)),max(sim3$k1dt)),
        ylim=c(0,max(c(sim3$k1dt,sim3$y_tot))))

summary(sim3)

# k1dt= delta_y / (y_tot - y) * 1/w
realSteps_Vp= diff(Vp_wNoise)  #delta_y
real_ytot_y= 1-Vp_wNoise   #y_tot - y
real_rel_dy= realSteps_Vp/real_ytot_y[2:length(real_ytot_y)] * 1/concs_out$Wp[2:length(real_ytot_y)]

plot(concs_out$time[2:length(concs_out$time)],real_rel_dy,main="Difference between steps",
     xlab='time',ylab='delta_y / (y_tot - y) * 1/Wp',ylim=c(-500,500))
abline(h=mean(sim3$k1dt),col='red',lwd=3)
abline(h=quantile(sim3$k1dt,.1),col='red',lwd=1,lty=2)
abline(h=quantile(sim3$k1dt,.9),col='red',lwd=1,lty=2)


##############################
#model (NO DEPHOSPHORYLATION; first order phosphorylation with saturation); sigma is relative

#specify model
modelString4= "
 model {
  for (i in 2:length(y)) {
    y[i] ~ dnorm(y_model[i], 1/sigma_rel[i]^2)  #likelihood   #instead of relative sigma, you could also normalize m by its distance from y_max
    y_model[i]= y[i-1]- k1dt*w[i-1]*y[i-1]+ k1dt*w[i-1]*y_tot
    #sigma_rel[i]= sigma/(y_tot-y[i-1])   ####################################
    sigma_rel[i]= sigma/(max(10^-6,y_tot-y[i-1])) #y_tot may be estimated < y[i-1]
  } 
  k1dt ~ dnorm(0, 1/0.05^2) #prior
  sigma ~ dgamma(1,2)  #mean: a/b; sd= sqrt(a)/b
  y_tot ~ dgamma(1,2)  #to be optimized; should be positive
 }
"

jagsModel4= jags.model( file= textConnection(modelString4), data=dataList3, inits=initList3,
                        n.chains=2, n.adapt=500)  #find MCMC sampler

codaSamples4= coda.samples(jagsModel4, variable.names=c('k1dt','sigma','y_tot'),n.iter=10000)
plot(codaSamples4)

diagMCMC(codaObject = codaSamples4, parName = 'k1dt')
diagMCMC(codaObject = codaSamples4, parName = 'y_tot')
dev.off()

#based on the simulated parameters, simulate traces
sim4= data.frame(codaSamples4[[1]])
head(sim4)
mu_sim= matrix(NA,nrow= dim(codaSamples4[[1]])[1],ncol=dim(concs_out)[1])
mu_sim[,1]= rep(Y0[6],dim(mu_sim)[1]) #starting value
rel_deltaY_sim= matrix(NA,nrow= dim(codaSamples4[[1]])[1],ncol=dim(concs_out)[1])
sigma_rel_sim= matrix(NA,nrow= dim(codaSamples4[[1]])[1],ncol=dim(concs_out)[1])
for(tpidx in 2:dim(mu_sim)[2]) { #one TP after another
  mu_sim[,tpidx]= mu_sim[,tpidx-1]- sim4$k1dt*concs_out$Wp[tpidx-1]*mu_sim[,tpidx-1] + 
    sim4$k1dt*concs_out$Wp[tpidx-1]*sim4$y_tot
  rel_deltaY_sim[,tpidx]= (mu_sim[,tpidx]-mu_sim[,tpidx-1])/ (sim4$y_tot-mu_sim[,tpidx])* 1/concs_out$Wp[tpidx-1]
  sigma_rel_sim[,tpidx]= sim4$sigma/(sim4$y_tot-mu_sim[,tpidx])  # -concs_out$Vp[tpidx]
}
par(mfrow=c(1,1))
plot(concs_out$time,concs_out$Vp)
for(i in 9000:9100){   #
  lines(concs_out$time,mu_sim[i,],col='red')
}

#DFs:
 #concs_out: input data  (rows=TPs; variables in columns)
 #sim4: parameter chains (rows= chain steps; parameters in columns)
 #mu_sim: simulated PLAUSIBLE VALUE of variables over time(rows: chain steps; columns: TPs)  (variation from sigma is not incorporated)


kde4= kde2d(sim4$k1dt,sim4$y_tot)
contour(kde4,xlab= 'k1dt',ylab='y_tot',xlim=c(min(c(0,sim4$k1dt)),max(sim4$k1dt)),
        ylim=c(0,max(c(sim4$k1dt,sim4$y_tot))))

summary(sim4)

# k1dt= delta_y / (y_tot - y) * 1/w
realSteps_Vp= diff(Vp_wNoise)  #delta_y
real_ytot_y= 1-Vp_wNoise   #y_tot - y
real_rel_dy= realSteps_Vp/real_ytot_y[2:length(real_ytot_y)] * 1/concs_out$Wp[2:length(real_ytot_y)]

plot(concs_out$time[2:length(concs_out$time)],real_rel_dy,main="Difference between steps",
     xlab='time',ylab='delta_y / (y_tot - y) * 1/Wp',ylim=c(-500,500))
lines(concs_out$time[2:length(concs_out$time)],sigma_rel_sim[1,2:length(concs_out$time)],col='red')
# abline(h=mean(sim4$k1dt),col='red',lwd=3)
# abline(h=quantile(sim4$k1dt,.1),col='red',lwd=1,lty=2)
# abline(h=quantile(sim4$k1dt,.9),col='red',lwd=1,lty=2)

plot(concs_out$time[2:length(concs_out$time)],realSteps_Vp,main="Difference between steps",
     xlab='time',ylab='delta_y ')
lines(concs_out$time[2:length(concs_out$time)],sigma_rel_sim[1,2:length(concs_out$time)],col='red')









###################### not used ##########################
#allow parameter range from 0 to 10
parRes= 100 #resolution of parameter Space
ran= seq(0,3,length.out=parRes)
l= layout(mat=matrix(c(0,1,0,2,3,4),nrow=2,ncol=3,byrow=TRUE),heights=c(3,3),widths=c(3,3,3))
par(oma=c(0,0,3,1)) #mfrow=c(1,3),

##############################
#model k_XZ (NO DEPHOSPHORYLATION)
#Prior: assume gamma-distribution with low mean and low effective sample size for k_X
k_XZ_prior= dgamma(x= ran,shape=1,rate=5)
PriLikPost_modes= matrix(NA,nrow=dim(concs_out)[1],ncol=3)
PriLikPost_modes[1,3]= mean(k_XZ_prior) #initial prior mean
j=4
for(i in 2:dim(concs_out)[1]) {
  delta_t= concs_out[i,1]-concs_out[i-1,1]
  
  #Plot showing parameter-choices within data; different from likelihood, which is data within distributions based on parameters
  plot(x= seq(0,3,.1),dnorm(x= seq(0,3,.1),mean=concs_out[i,j],sd=0.2),xlab='c(Z)',main='conc.') #distribution around observed value
  abline(v=c(concs_out[i-1,j]+ran[1]*concs_out[i-1,5]*delta_t,   #range covered by parameter range
             concs_out[i-1,j]+ran[parRes/2]*concs_out[i-1,5]*delta_t,
             concs_out[i-1,j]+ran[parRes]*concs_out[i-1,5]*delta_t),
         col='blue',lty=2)
  
  plot(ran,k_XZ_prior,type='l',main='Prior k_XZ') #,xlim=c(0,1)
  
  #Likelihood: Normal distribution around calculated value with arbitrary SD 0.02
  #x= observed value; mean= Calculated value for each value on parameter range
  #k_XZ_lik= dnorm(x= concs_out[i-1,j]+ran*concs_out[i-1,5]*delta_t,mean=concs_out[i,j],sd=0.02)
  k_XZ_lik= dnorm(x= concs_out[i,j],mean=concs_out[i-1,j]+ran*concs_out[i-1,5]*delta_t,sd=0.02)
  plot(ran,k_XZ_lik,type='l',main='Likelihood k_XZ')
  
  #Posterior
  k_XZ_post= k_XZ_prior * k_XZ_lik /sum(k_XZ_prior * k_XZ_lik)
  plot(ran,k_XZ_post,type='l',main='Posterior k_XZ')
  
  mtext(text=paste('round',i),side=3,line=1,outer=TRUE)
  
  priorMode= ran[ which(max(k_XZ_prior)==k_XZ_prior)[1] ]
  likMode= ran[ which(max(k_XZ_lik)==k_XZ_lik)[1] ]
  postMode= ran[ which(max(k_XZ_post)==k_XZ_post)[1] ]
  PriLikPost_modes[i,]= c(priorMode,likMode,postMode)
  k_XZ_prior= k_XZ_post
}

par(mfrow=c(1,2))
plot(1:dim(PriLikPost_modes)[1],PriLikPost_modes[,3],type='l',xlab='iteration',
     ylab='k_XZ Posterior mode',main='Posterior')
plot(1:dim(PriLikPost_modes)[1],PriLikPost_modes[,2],type='l',xlab='iteration',
     ylab='k_XZ Likelihood mode',main='Likelihood')

###############################
#model "Up" (constant phosphorylation, linear dephosphorylation rate)
#Prior: assume gamma-distribution with low mean and low effective sample size for k_X
ran1= seq(0,3,length.out=parRes)
ran2= seq(0,3,length.out=parRes)

k_WU_prior= dgamma(x= ran1,shape=1,rate=5)
k_dU_prior= dgamma(x= ran2,shape=1,rate=5)
prior_2D= 
  matrix(rep(k_WU_prior,times=parRes),nrow=parRes, byrow=TRUE) *
  matrix(rep(k_dU_prior,times=parRes),ncol=parRes, byrow=FALSE)
#contour(x= ran1, y=ran2, z=prior_2D, main= 'Prior',xlab='k_WU',ylab='k_dU')
k_WU_LikPost_modes= matrix(NA,nrow=dim(concs_out)[1],ncol=2)
k_WU_LikPost_modes[1,2]= ran1[which(k_WU_prior==max( k_WU_prior))[1]]
k_dU_LikPost_modes= matrix(NA,nrow=dim(concs_out)[1],ncol=2)
k_dU_LikPost_modes[1,2]= ran2[which(k_dU_prior==max( k_dU_prior))[1]]

j=6  #colNo of U in concs_out
#ani.record(reset = TRUE) #reset
for(i in 2:dim(concs_out)[1]) {
  delta_t= concs_out[i,1]-concs_out[i-1,1]
  
  #Plot showing parameter-choices within data; different from likelihood, which is data within distributions based on parameters
  plot(x= seq(0,3,.1),dnorm(x= seq(0,3,.1),mean=concs_out[i,j],sd=0.2),xlab='c(U)',main='conc.') #distribution around observed value
  abline(v= c(concs_out[i-1,j]+ ran1[1]*concs_out[i-1,5]*delta_t - ran2[1]*concs_out[i-1,j]*delta_t,  #kinase parameter range at low PPTase 
              concs_out[i-1,j]+ ran1[parRes]*concs_out[i-1,5]*delta_t - ran2[1]*concs_out[i-1,j]*delta_t), col='blue',lty=2)
  abline(v= c(concs_out[i-1,j]+ ran1[1]*concs_out[i-1,5]*delta_t - ran2[parRes]*concs_out[i-1,j]*delta_t,  #kinase parameter range at high PPTase
              concs_out[i-1,j]+ ran1[parRes]*concs_out[i-1,5]*delta_t - ran2[parRes]*concs_out[i-1,j]*delta_t), col='red',lty=2)
  
  contour(x= ran1, y=ran2, z=prior_2D, main= 'Prior',xlab='k_WU',ylab='k_dU')
  
  #Likelihood: Normal distribution around calculated value with arbitrary SD 0.02
  #x= observed value; mean= Calculated value for each value on parameter range
  lik_2D= matrix(NA,nrow=parRes, ncol=parRes)
  for(a1 in 1:parRes){
    for(a2 in 1:parRes){
      lik_2D[a1,a2]= dnorm(x= concs_out[i,j],
                           concs_out[i-1,j]+ ran1[a1]*concs_out[i-1,5]*delta_t - ran2[a2]*concs_out[i-1,j]*delta_t,
                           sd=0.02)
    }
  }
  contour(x= ran1, y=ran2, z=lik_2D, main= 'Likelihood',xlab='k_WU',ylab='k_dU')
  
  #Posterior
  post_2D= prior_2D * lik_2D /sum(prior_2D * lik_2D)
  contour(x= ran1, y=ran2, z=post_2D, main= 'Posterior',xlab='k_WU',ylab='k_dU')
  
  mtext(text=paste('round',i),side=3,line=1,outer=TRUE)
  #ani.record()
  prior_2D= post_2D
  
  #marginalize:
  k_WU_likMarg= apply(lik_2D,1,sum)
  k_dU_likMarg= apply(lik_2D,2,sum)
  k_WU_postMarg= apply(post_2D,1,sum)
  k_dU_postMarg= apply(post_2D,2,sum)
  
  k_WU_likMode= ran1[which(k_WU_likMarg==max( k_WU_likMarg))[1]]
  k_dU_likMode= ran2[which(k_dU_likMarg==max( k_dU_likMarg))[1]]
  k_WU_postMode= ran1[which(k_WU_postMarg==max( k_WU_postMarg))[1]]
  k_dU_postMode= ran2[which(k_dU_postMarg==max( k_dU_postMarg))[1]]
  
  k_WU_LikPost_modes[i,]= c(k_WU_likMode,k_WU_postMode)
  k_dU_LikPost_modes[i,]= c(k_dU_likMode,k_dU_postMode)
} 
# oopt= ani.options(interval = 2, loop=1)   #loop decides how often GIF will loop
# saveGIF(ani.replay(),movie.name=paste(plotDir,plotname'.gif',sep=''))

# Plot mode of marginal likelihood and marginal posterior for each parameter  
par(mfrow=c(2,2))
plot(ran1, k_WU_likMarg,main='Marginal likelihood',xlab='k_WU',ylab='P( D | k_WU)')
plot(ran2, k_dU_likMarg,main='Marginal likelihood',xlab='k_dU',ylab='P( D | k_dU)')
plot(ran1, k_WU_postMarg,main='Marginal posterior',xlab='k_WU',ylab='P( k_WU |D )')
plot(ran2, k_dU_postMarg,main='Marginal posterior',xlab='k_dU',ylab='P( k_dU |D )')

# Plot mode of marginal likelihoods over rounds of adding data
par(mfrow=c(1,2))
plot(type='n',x=0,y=0,xlim=c(0,dim(concs_out)[1]),ylim=c(0,max(k_WU_LikPost_modes,na.rm=TRUE)),
     xlab='Iteration',ylab='k_WU Marg. Mode')
lines(1:dim(concs_out)[1], k_WU_LikPost_modes[,1])
lines(1:dim(concs_out)[1], k_WU_LikPost_modes[,2],col='blue')
plot(type='n',x=0,y=0,xlim=c(0,dim(concs_out)[1]),ylim=c(0,max(k_dU_LikPost_modes,na.rm=TRUE)),
     xlab='Iteration',ylab='k_dU Marg. Mode')
lines(1:dim(concs_out)[1], k_dU_LikPost_modes[,1])
lines(1:dim(concs_out)[1], k_dU_LikPost_modes[,2],col='blue')


###############################
#model "Ap" (saturating 1-order phosphorylation (w/o explicit kinase), assuming Atotal= 1; 1-order dephosphorylation)
#Prior: assume gamma-distribution with low mean and low effective sample size for k_X
ran1= seq(0,3,length.out=parRes) #kinase parameter range
ran2= seq(0,3,length.out=parRes) #PPTase parameter range

k_XA_prior= dgamma(x= ran1,shape=1,rate=5)
k_dA_prior= dgamma(x= ran2,shape=1,rate=5)
prior_2D= 
  matrix(rep(k_XA_prior,times=parRes),nrow=parRes, byrow=TRUE) *
  matrix(rep(k_dA_prior,times=parRes),ncol=parRes, byrow=FALSE)
k_XA_LikPost_modes= matrix(NA,nrow=dim(concs_out)[1],ncol=2)
k_XA_LikPost_modes[1,2]= ran1[which(k_XA_prior==max( k_XA_prior))[1]]
k_dA_LikPost_modes= matrix(NA,nrow=dim(concs_out)[1],ncol=2)
k_dA_LikPost_modes[1,2]= ran2[which(k_dA_prior==max( k_dA_prior))[1]]

j=2  #colNo of U in concs_out
l1= layout(mat=matrix(c(0,1,2,3),nrow=2,ncol=2,byrow=TRUE),heights=c(3,3),widths=c(3,3))
ani.record(reset = TRUE) #reset
for(i in 2:dim(concs_out)[1]) {
  delta_t= concs_out[i,1]-concs_out[i-1,1]
  
  #Plot showing parameter-choices within data; different from likelihood, which is data within distributions based on parameters
  #showing lowest and highest allowed parameter
  plot(x= seq(0,3,.1),dnorm(x= seq(0,3,.1),mean=concs_out[i,j],sd=0.2),xlab='c(A)',main='conc.') #distribution around observed value
  abline(v= c(concs_out[i-1,j]+ ran1[1]*(1-concs_out[i-1,j])*delta_t      - ran2[1]*concs_out[i-1,j]*delta_t,  #kinase parameter range at low PPTase 
              concs_out[i-1,j]+ ran1[parRes]*(1-concs_out[i-1,j])*delta_t - ran2[1]*concs_out[i-1,j]*delta_t), col='blue',lty=2)
  abline(v= c(concs_out[i-1,j]+ ran1[1]*(1-concs_out[i-1,j])*delta_t      - ran2[parRes]*concs_out[i-1,j]*delta_t,  #kinase parameter range at high PPTase
              concs_out[i-1,j]+ ran1[parRes]*(1-concs_out[i-1,j])*delta_t - ran2[parRes]*concs_out[i-1,j]*delta_t), col='red',lty=2)
  
  # contour(x= ran1, y=ran2, z=prior_2D, main= 'Prior',xlab='k_XA',ylab='k_dA',
  #         xlim=c(0,1),ylim=c(0,1))
  # persp(x= ran1, y=ran2, z=prior_2D, main= 'Prior',xlab='k_XA',ylab='k_dA',
  #         xlim=c(0,1),ylim=c(0,1))
  
  #Likelihood: Normal distribution around calculated value with arbitrary SD 0.02
  #x= observed value; mean= Calculated value for each value on parameter range
  lik_2D= matrix(NA,nrow=parRes, ncol=parRes)
  for(a1 in 1:parRes){
    for(a2 in 1:parRes){
      lik_2D[a1,a2]= dnorm(x= concs_out[i,j],
                           concs_out[i-1,j]+ ran1[a1]*(1-concs_out[i-1,j])*delta_t - ran2[a2]*concs_out[i-1,j]*delta_t,
                           sd=0.02)
    }
  }
  contour(x= ran1, y=ran2, z=lik_2D, main= 'Likelihood',xlab='k_XA',ylab='k_dA',
          xlim=c(0,1),ylim=c(0,1))
  # persp(x= ran1, y=ran2, z=lik_2D, main= 'Likelihood',xlab='k_XA',ylab='k_dA',
  #         xlim=c(0,1),ylim=c(0,1),theta=30,phi=30,,col='lightblue')  #expand=0.5
  
  #Posterior
  rel_post_2D= prior_2D * lik_2D #relative posterior w/o division by P(D)
  contour(x= ran1, y=ran2, z=rel_post_2D, main= 'Posterior',xlab='k_XA',ylab='k_dA',
          xlim=c(0,1),ylim=c(0,1))
  
  mtext(text=paste('round',i),side=3,line=1,outer=TRUE)
  ani.record()
  prior_2D= rel_post_2D
  
  #marginalize:
  k_XA_likMarg= apply(lik_2D,1,sum)
  k_dA_likMarg= apply(lik_2D,2,sum)
  k_XA_postMarg= apply(rel_post_2D,1,sum)
  k_dA_postMarg= apply(rel_post_2D,2,sum)
  (postFullMarg= sum(k_XA_postMarg))  #full Posterior
  (postFullMarg2= sum(k_dA_postMarg))
  
  k_XA_likMode= ran1[which(k_XA_likMarg==max( k_XA_likMarg))[1]]
  k_dA_likMode= ran2[which(k_dA_likMarg==max( k_dA_likMarg))[1]]
  k_XA_postMode= ran1[which(k_XA_postMarg==max( k_XA_postMarg))[1]]
  k_dA_postMode= ran2[which(k_dA_postMarg==max( k_dA_postMarg))[1]]
  
  k_XA_LikPost_modes[i,]= c(k_XA_likMode,k_XA_postMode)
  k_dA_LikPost_modes[i,]= c(k_dA_likMode,k_dA_postMode)
} 
oopt= ani.options(interval = .5, loop=1)   #loop decides how often GIF will loop
saveGIF(ani.replay(),movie.name=paste('params2D.gif',sep=''))

par(mfrow=c(1,1))
persp(x= ran1, y=ran2, z=lik_2D, main= 'Likelihood',xlab='k_XA',ylab='k_dA',
      xlim=c(0,1),ylim=c(0,1),theta=30,phi=30,expand=0.5,col='lightblue')  #
persp(x= ran1, y=ran2, z=rel_post_2D, main= 'Relative Posterior',xlab='k_XA',ylab='k_dA',
      xlim=c(0,1),ylim=c(0,1),theta=30,phi=30,expand=0.5,col='lightblue')

# Plot mode of marginal likelihood and marginal posterior for each parameter  
par(mfrow=c(2,2))
plot(ran1, k_XA_likMarg,main='Marginal likelihood',xlab='k_XA',ylab='P( D | k_XA)')
plot(ran2, k_dA_likMarg,main='Marginal likelihood',xlab='k_dA',ylab='P( D | k_dA)')
plot(ran1, k_XA_postMarg,main='Marginal posterior',xlab='k_XA',ylab='P( k_XA |D )')
plot(ran2, k_dA_postMarg,main='Marginal posterior',xlab='k_dA',ylab='P( k_dA |D )')

# Plot mode of marginal likelihoods over rounds of adding data
par(mfrow=c(1,2))
plot(type='n',x=0,y=0,xlim=c(0,dim(concs_out)[1]),ylim=c(0,max(k_XA_LikPost_modes,na.rm=TRUE)),
     xlab='Iteration',ylab='k_XA Marg. Mode')
lines(1:dim(concs_out)[1], k_XA_LikPost_modes[,1])
lines(1:dim(concs_out)[1], k_XA_LikPost_modes[,2],col='blue')
plot(type='n',x=0,y=0,xlim=c(0,dim(concs_out)[1]),ylim=c(0,max(k_dA_LikPost_modes,na.rm=TRUE)),
     xlab='Iteration',ylab='k_dA Marg. Mode')
lines(1:dim(concs_out)[1], k_dA_LikPost_modes[,1])
lines(1:dim(concs_out)[1], k_dA_LikPost_modes[,2],col='blue')

#############################################
#add noise
t_vec= concs_out[,1]
Y_concs= as.matrix(concs_out[,c(2,3)])  #needs to be matrix
#Y_concs= matrix(1:60,nrow=30,ncol=2,byrow = TRUE)
set.seed(123)
sigm= 0.05
Y_3D= array(rep(NA,dim(Y_concs)[1]*2*3), dim=c(dim(Y_concs)[1],2,3))   #data for 3 replicates; replicates on third dimension
#Y_3D[,,1]= Y_concs

for (i in 1:3){
  Y_curr= cbind.data.frame(
              A= Y_concs[,1]+ rnorm(dim(concs_out)[1],mean=0,sd=sigm),
              B= Y_concs[,2]+ rnorm(dim(concs_out)[1],mean=0,sd=sigm) )
  Y_3D[,,i]= as.matrix(Y_curr)
}
p1<- plot(type='n',x=0,y=0,xlim=c(0,max(t_vec)),ylim=c(0,1))
p2<- plot(type='n',x=0,y=0,xlim=c(0,max(t_vec)),ylim=c(0,1))
for (i in 1:3){
  p1<- p1+points(t_vec,Y_3D[,1,i])
  p2<- p2+points(t_vec,Y_3D[,1,i])
}

k= 1   #k to explore
head(concs_out)
deltaYB_t4= concs_out[concs_out$time ==4,2] - concs_out[concs_out$time ==3,2]   #change in B btw TPs
YA_t3= concs_out[concs_out$time==3,1]


deltaYB_t4= deltaYB_t4+ rnorm(3,mean=0,sd=sigm) 
YA_t3= YA_t3+ rnorm(3,mean=0,sd= k*sigm)

#probability of observed data under model
P= prod( dnorm(deltaYB_t4,mean= k*YA_t3,sd=sigm) )

#20230502: one vectorized instead of multiple scalar ODE functions
