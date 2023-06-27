rm(list=ls())
cat("\014")  #is the code to send ctrl+L to the console and therefore will clear the screen
plot.new()

# to add
#allow diff. timesteps
#do you arrive at the same result if you only take at subset of TPs -> is prior overwhelmed
#fixed vs. relative noise
#link function other than normal
#check agreement between chains
#relative phosphorylation -> max=1

#replicates
#dephosphorylation
#estimation of SD
#no negative concentration
#JAGS
#noise

library('deSolve')
library('animation')
library('rjags')
library('tidyverse')
library('plotly')
library('MASS')
source('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/DBDA2E-utilities.R')


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
k_W2V= 0.5

k_dA= 1 #dephosphorylation
k_dB= 1
k_dC= 1
k_dZ= 0
k_dW= 0
k_dU= .05
k_dV= 0
k_dV2= 1

# initial values
Y0= c(Ap=0,Bp=0,Zp=0.01,Wp=0.01,Up=0.1,Vp=0.1,W2p=0.01,V2p=0.01,V3p=0.01)  #

parmsX= matrix(c(k_XA,k_AB,k_BC,k_WZ,k_XW, k_WU, k_WV, k_W2V,
                 k_dA,k_dB,k_dC,k_dZ,k_dW, k_dU, k_dV, k_dV2),ncol=2,byrow=FALSE)  #col1= phosphorylation; col2: dephosphorylation
dYp_dt= function(t,Y,parms) {    #Y= Ap, Bp; parmsX: col1: k_XA, k_AB, k_BC; col2: k_dA, k_dB, k_dC
  dAp_dt= parms[1,1]*(1-Y[1]) - parms[1,2]*Y[1]  #saturated 1-order phosphorylation (w/o explicit kinase); 
                        # 0-order dephosphorylation; assuming Atot = 1
  dBp_dt= Y[1]*parms[2,1]*(1-Y[2]) - parms[2,2]*Y[2]
  dZp_dt= parms[4,1]*Y[4]  #0-order phosphorylation rate
  dWp_dt= 0 #constant kinase
  dUp_dt= parms[6,1]*Y[4] - parms[6,2]*Y[5] #0-order phosphorylation and first-order dephosphorylation
  dVp_dt= Y[4]*parms[7,1]*(1-Y[6])
  dW2p_dt= 0 #constant kinase
  dV2p_dt= Y[4]*parms[7,1]*(1-Y[8]) + Y[7]*parms[8,1]*(1-Y[8])
  dV3p_dt= Y[4]*parms[7,1]*(1-Y[9]) - Y[7]*parms[8,2]*Y[9]   #kinase and PPTase
  return(list(c(dAp_dt, dBp_dt, dZp_dt, dWp_dt, dUp_dt, dVp_dt,dW2p_dt,dV2p_dt,dV3p_dt))) #
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
plot(concs_out$t,concs_out$V2p,ylim=c(0,1),main='V2')
plot(concs_out$t,concs_out$V3p,ylim=c(0,1),main='V3')
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
    Zp_wNoise= concs_out$Zp + rnorm(length(concs_out$Zp),mean=0,sd= sigma_real) #absolute noise
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

#relative noise
real_relSD= 0.001
Wp= concs_out$Wp
Vp= concs_out$Vp
Vp_wNoise= Vp+ sapply(concs_out$Vp, function(x){rnorm(n=1,mean=0,sd= x*real_relSD)})  #relative
plot(concs_out$t,Vp_wNoise)

#simulate as difference equation, rather than differential equation
V_fromDiff= rep(NA,1001)
V_fromDiff= Y0[6]
for(tidx in 2:1001) {
  V_fromDiff[tidx]= V_fromDiff[tidx-1]+ Y0[4]* parmsX[7,1] * (1-V_fromDiff[tidx-1]) #Wp is constant; delta_t=1
}
lines(concs_out$t,V_fromDiff,col='green')

dataList3= list(
  y= Vp_wNoise,
  w= concs_out$Wp,
  w2= concs_out$W2p
)

#specify model
modelString3= "
 model {
  for (i in 2:length(y)) {
    y[i] ~ dnorm(m[i], 1/sigma^2)  #likelihood
    m[i]= y[i-1]- k1dt*w[i-1]*y[i-1]+ k1dt*w[i-1]*y_tot
  } 
  k1dt ~ dnorm(0, 1/1^2) #prior
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

jagsModel3= jags.model( file= textConnection(modelString3), data=dataList3, inits=initList3,
                       n.chains=2, n.adapt=500)  #find MCMC sampler

codaSamples3= coda.samples(jagsModel3, variable.names=c('k1dt','sigma','y_tot'),n.iter=10000)
plot(codaSamples3)

diagMCMC(codaObject = codaSamples3, parName = 'k1dt')
dev.off()
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


kde3= MASS::kde2d(sim3$k1dt,sim3$y_tot)
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
    y[i] ~ dnorm(y_model[i], 1/SD_abs[i]^2)  #likelihood   #instead of relative sigma, you could also normalize m by its distance from y_max
    y_model[i]= y[i-1]- k1dt*w[i-1]*y[i-1]+ k1dt*w[i-1]*y_tot
    SD_abs[i]= SD_rel * y[i-1]
  } 
  k1dt ~ dnorm(0, 1/1^2) #prior
  SD_rel ~ dgamma(1,2)  #mean: a/b; sd= sqrt(a)/b
  y_tot ~ dgamma(1,2)  #to be optimized; should be positive
 }
"

initList4= function()  { #initial values generated via function: mu centered around 0; sigma in range 0-0.1
  set.seed(123)
  initMu= rnorm(1,0,0.1)
  initSigma= runif(1,0,0.1)
  init_ytot= 1
  return( list( k1dt= initMu, SD_rel=initSigma, y_tot= init_ytot))
}

jagsModel4= jags.model( file= textConnection(modelString4), data=dataList3, inits=initList4,
                        n.chains=2, n.adapt=500)  #find MCMC sampler

codaSamples4= coda.samples(jagsModel4, variable.names=c('k1dt','SD_rel','y_tot'),n.iter=10000)
plot(codaSamples4)

diagMCMC(codaObject = codaSamples4, parName = 'k1dt')
dev.off()
diagMCMC(codaObject = codaSamples4, parName = 'y_tot')
dev.off()

#based on the simulated parameters, simulate traces
sim4= data.frame(codaSamples4[[1]])
head(sim4)
mu_sim= matrix(NA,nrow= dim(codaSamples4[[1]])[1],ncol=dim(concs_out)[1])
mu_sim[,1]= rep(Y0[6],dim(mu_sim)[1]) #starting value
deltaY_sim= matrix(NA,nrow= dim(codaSamples4[[1]])[1],ncol=dim(concs_out)[1])
rel_deltaY_sim= matrix(NA,nrow= dim(codaSamples4[[1]])[1],ncol=dim(concs_out)[1])
SD_abs_sim= matrix(NA,nrow= dim(codaSamples4[[1]])[1],ncol=dim(concs_out)[1])
for(tpidx in 2:dim(mu_sim)[2]) { #one TP after another
  mu_sim[,tpidx]= mu_sim[,tpidx-1]- sim4$k1dt*concs_out$Wp[tpidx-1]*mu_sim[,tpidx-1] + 
    sim4$k1dt*concs_out$Wp[tpidx-1]*sim4$y_tot
  deltaY_sim[,tpidx]= (mu_sim[,tpidx]-mu_sim[,tpidx-1])
  rel_deltaY_sim[,tpidx]= deltaY_sim[,tpidx]/ (sim4$y_tot-mu_sim[,tpidx])* 1/concs_out$Wp[tpidx-1]
  SD_abs_sim[,tpidx]= sim4$SD_rel * mu_sim[,tpidx]  # -concs_out$Vp[tpidx]
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
# abline(h=mean(sim4$k1dt),col='red',lwd=3)
# abline(h=quantile(sim4$k1dt,.1),col='red',lwd=1,lty=2)
# abline(h=quantile(sim4$k1dt,.9),col='red',lwd=1,lty=2)


deltaY_plusSD_sim= deltaY_sim + SD_abs_sim
deltaY_minusSD_sim= deltaY_sim - SD_abs_sim
plot(concs_out$time[2:length(concs_out$time)],realSteps_Vp,main="Difference between steps",
     xlab='time',ylab='delta_y ')
for(i in 9000:9100) {   #dim(deltaY_sim)[1]
    lines(concs_out$time[2:length(concs_out$time)],deltaY_plusSD_sim[i,2:length(concs_out$time)],col='red')
    lines(concs_out$time[2:length(concs_out$time)],deltaY_minusSD_sim[i,2:length(concs_out$time)],col='red')
}
legend('topright',legend=c( paste('Rel. SD:',real_relSD),'b'))
#lines(concs_out$time[2:length(concs_out$time)],deltaY_sim[1,2:length(concs_out$time)],col='red')
#lines(concs_out$time[2:length(concs_out$time)],SD_abs_sim[1,2:length(concs_out$time)],col='red')



##############################
#model5 (first order phosphorylation with saturation, two kinases/phosphatases); sigma is relative

#relative noise
real_relSD= 0.001
Wp= concs_out$Wp
W2p= concs_out$W2p
V2p= concs_out$V2p
V2p_wNoise= V2p+ sapply(concs_out$V2p, function(x){rnorm(n=1,mean=0,sd= x*real_relSD)})  #relative
plot(concs_out$t,V2p_wNoise)

#simulate as difference equation, rather than differential equation
V2_fromDiff= rep(NA,1001)
V2_fromDiff= Y0[6]
for(tidx in 2:1001) {
  V2_fromDiff[tidx]= V2_fromDiff[tidx-1]+ (Y0[4]*parmsX[7,1]+Y0[7]*parmsX[8,1]) * (1-V2_fromDiff[tidx-1]) #Wp, W2p are constant; delta_t=1
}
lines(concs_out$t,V2_fromDiff,col='green')

dataList5= list(
  y= V2p_wNoise,
  w= concs_out$Wp,
  w2= concs_out$W2p
)

#specify model
modelString5= "
 model {
  for (i in 2:length(y)) {
    y[i] ~ dnorm(y_model[i], 1/SD_abs[i]^2)  #likelihood   #instead of relative sigma, you could also normalize m by its distance from y_max
    y_model[i]= y[i-1]- (k1dt*w[i-1]+k2dt*w2[i-1])*y[i-1]+ (k1dt*w[i-1]+k2dt*w2[i-1])*y_tot
    SD_abs[i]= SD_rel * y[i-1]
  } 
  k1dt ~ dnorm(0, 1/1^2) #prior
  k2dt ~ dnorm(0, 1/1^2)
  SD_rel ~ dgamma(1,2)  #mean: a/b; sd= sqrt(a)/b
  y_tot ~ dgamma(1,2)  #to be optimized; should be positive
 }
"

initList5= function()  { #initial values generated via function: mu centered around 0; sigma in range 0-0.1
  set.seed(123)
  init_k1dt= rnorm(1,0,0.1)
  init_k2dt= rnorm(1,0,0.1)
  initSigma= runif(1,0,0.1)
  init_ytot= 1
  return( list( k1dt= init_k1dt, k2dt= init_k2dt, SD_rel=initSigma, y_tot= init_ytot))
}

jagsModel5= jags.model( file= textConnection(modelString5), data=dataList5, inits=initList5,
                        n.chains=2, n.adapt=500)  #find MCMC sampler

codaSamples5= coda.samples(jagsModel5, variable.names=c('k1dt', 'k2dt','SD_rel','y_tot'),n.iter=10000)
plot(codaSamples5)

diagMCMC(codaObject = codaSamples5, parName = 'k1dt')
dev.off()
diagMCMC(codaObject = codaSamples5, parName = 'y_tot')
dev.off()

#based on the simulated parameters, simulate traces
sim5= data.frame(codaSamples5[[1]])
head(sim5)
mu_sim= matrix(NA,nrow= dim(codaSamples5[[1]])[1],ncol=dim(concs_out)[1])
mu_sim[,1]= rep(Y0[8],dim(mu_sim)[1]) #starting value
deltaY_sim= matrix(NA,nrow= dim(codaSamples5[[1]])[1],ncol=dim(concs_out)[1])
SD_abs_sim= matrix(NA,nrow= dim(codaSamples5[[1]])[1],ncol=dim(concs_out)[1])
for(tpidx in 2:dim(mu_sim)[2]) { #one TP after another
  mu_sim[,tpidx]= mu_sim[,tpidx-1]- 
    (sim5$k1dt*concs_out$Wp[tpidx-1]+sim5$k2dt*concs_out$W2p[tpidx-1]) * mu_sim[,tpidx-1] + 
    (sim5$k1dt*concs_out$Wp[tpidx-1]+sim5$k2dt*concs_out$W2p[tpidx-1]) *sim5$y_tot
  deltaY_sim[,tpidx]= (mu_sim[,tpidx]-mu_sim[,tpidx-1])
  SD_abs_sim[,tpidx]= sim5$SD_rel * mu_sim[,tpidx]  # -concs_out$Vp[tpidx]
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


kde5= kde2d(sim5$k1dt,sim5$y_tot)
contour(kde5,xlab= 'k1dt',ylab='y_tot',xlim=c(min(c(0,sim5$k1dt)),max(sim5$k1dt)),
        ylim=c(0,max(c(sim5$k1dt,sim5$y_tot))))

summary(sim5)

# k1dt= delta_y / (y_tot - y) * 1/w
realSteps_Vp= diff(V2p_wNoise)  #delta_y

deltaY_plusSD_sim= deltaY_sim + SD_abs_sim
deltaY_minusSD_sim= deltaY_sim - SD_abs_sim
plot(concs_out$time[2:length(concs_out$time)],realSteps_Vp,main="Difference between steps",
     xlab='time',ylab='delta_y ')
for(i in 9000:9100) {   #dim(deltaY_sim)[1]
  lines(concs_out$time[2:length(concs_out$time)],deltaY_plusSD_sim[i,2:length(concs_out$time)],col='red')
  lines(concs_out$time[2:length(concs_out$time)],deltaY_minusSD_sim[i,2:length(concs_out$time)],col='red')
}
legend('topright',legend=c( paste('Rel. SD:',real_relSD),'+/- SD'),col=c(NA,'red'),pch=c(NA,'-'))
#lines(concs_out$time[2:length(concs_out$time)],deltaY_sim[1,2:length(concs_out$time)],col='red')
#lines(concs_out$time[2:length(concs_out$time)],SD_abs_sim[1,2:length(concs_out$time)],col='red')



##############################
#model6 (first order phosphorylation with saturation, first order dephosphorylation); sigma is relative

#relative noise
real_relSD= 0.001
Wp= concs_out$Wp
W2p= concs_out$W2p #PPTase
V3p= concs_out$V3p
V3p_wNoise= V3p+ sapply(concs_out$V3p, function(x){rnorm(n=1,mean=0,sd= x*real_relSD)})  #relative
plot(concs_out$t,V3p_wNoise,ylim=c(0,1))

dataList6= list(
  y= V3p_wNoise,
  w= concs_out$Wp,
  w2= concs_out$W2p
)

#specify model
modelString6= "
 model {
  for (i in 2:length(y)) {
    y[i] ~ dnorm(y_model[i], 1/SD_abs[i]^2)  #likelihood   #instead of relative sigma, you could also normalize m by its distance from y_max
    y_model[i]= y[i-1]- k1dt*w[i-1]*y[i-1] + k1dt*w[i-1]*y_tot - k2dt*w2[i-1]*y[i-1]
    SD_abs[i]= SD_rel * y[i-1]
  } 
  k1dt ~ dnorm(0, 1/1^2) #prior
  k2dt ~ dnorm(0, 1/1^2)
  SD_rel ~ dgamma(1,2)  #mean: a/b; sd= sqrt(a)/b
  y_tot ~ dgamma(1,2)  #to be optimized; should be positive
 }
"

initList6= function()  { #initial values generated via function: mu centered around 0; sigma in range 0-0.1
  set.seed(123)
  init_k1dt= rnorm(1,0,0.1)
  init_k2dt= rnorm(1,0,0.1)
  initSigma= runif(1,0,0.1)
  init_ytot= 1
  return( list( k1dt= init_k1dt, k2dt= init_k2dt, SD_rel=initSigma, y_tot= init_ytot))
}

jagsModel6= jags.model( file= textConnection(modelString6), data=dataList6, inits=initList6,
                        n.chains=2, n.adapt=500)  #find MCMC sampler

codaSamples6= coda.samples(jagsModel6, variable.names=c('k1dt', 'k2dt','SD_rel','y_tot'),n.iter=10000)
plot(codaSamples6)

diagMCMC(codaObject = codaSamples6, parName = 'k1dt')
dev.off()
diagMCMC(codaObject = codaSamples6, parName = 'y_tot')
dev.off()

#based on the simulated parameters, simulate traces
sim6= data.frame(codaSamples6[[1]])
head(sim6)
mu_sim= matrix(NA,nrow= dim(codaSamples6[[1]])[1],ncol=dim(concs_out)[1])
mu_sim[,1]= rep(Y0[8],dim(mu_sim)[1]) #starting value
deltaY_sim= matrix(NA,nrow= dim(codaSamples6[[1]])[1],ncol=dim(concs_out)[1])
SD_abs_sim= matrix(NA,nrow= dim(codaSamples6[[1]])[1],ncol=dim(concs_out)[1])
for(tpidx in 2:dim(mu_sim)[2]) { #one TP after another
  mu_sim[,tpidx]= mu_sim[,tpidx-1]-    
    sim6$k1dt*concs_out$Wp[tpidx-1] * mu_sim[,tpidx-1] + 
    sim6$k1dt*concs_out$Wp[tpidx-1] *sim6$y_tot -
    sim6$k2dt*concs_out$W2p[tpidx-1] *mu_sim[,tpidx-1]
  deltaY_sim[,tpidx]= (mu_sim[,tpidx]-mu_sim[,tpidx-1])
  SD_abs_sim[,tpidx]= sim6$SD_rel * mu_sim[,tpidx]
}
par(mfrow=c(1,1))
plot(concs_out$time,concs_out$V3p, ylim= c(0,1))
for(i in 9000:9100){   #
  lines(concs_out$time,mu_sim[i,],col='red')
}

kde6= kde2d(sim6$k1dt,sim6$y_tot)
contour(kde6,xlab= 'k1dt',ylab='y_tot',xlim=c(min(c(0,sim6$k1dt)),max(sim6$k1dt)),
        ylim=c(0,max(c(sim6$k1dt,sim6$y_tot))))

kde6b= kde2d(sim6$k1dt,sim6$k2dt)
contour(kde6b,xlab= 'k1dt',ylab='k2dt',xlim=c(min(c(0,sim6$k1dt)),max(sim6$k1dt)),
        ylim=c(min(c(0,sim6$k2dt)),max(sim6$k2dt)))
summary(sim6)

# k1dt= delta_y / (y_tot - y) * 1/w
realSteps_Vp= diff(V3p_wNoise)  #delta_y
