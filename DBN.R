rm(list=ls())
cat("\014")  #is the code to send ctrl+L to the console and therefore will clear the screen
plot.new()

# to add
 #dephosphorylation
 #estimation of SD
 #replicates
 #no negative concentration
 #relative phosphorylation -> max=1
 #JAGS
 #noise

library('deSolve')
library('animation')
library('rjags')
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

k_dA= 1 #dephosphorylation
k_dB= 1
k_dC= 1
k_dZ= 0
k_dW= 0
k_dU= .05

# initial values
Y0= c(Ap=0,Bp=0,Zp=0.01,Wp=0.01,Up=0.1)

parmsX= matrix(c(k_XA,k_AB,k_BC,k_WZ,k_XW, k_WU,
                 k_dA,k_dB,k_dC,k_dZ,k_dW, k_dU),ncol=2,byrow=FALSE)  #col1= phosphorylation; col2: dephosphorylation
dYp_dt= function(t,Y,parms) {    #Y= Ap, Bp; parmsX: col1: k_XA, k_AB, k_BC; col2: k_dA, k_dB, k_dC
  dAp_dt= parms[1,1]*(1-Y[1]) - parms[1,2]*Y[1]  #saturated 1-order phosphorylation (w/o explicit kinase); 
                        # 0-order dephosphorylation; assuming Atot = 1
  dBp_dt= Y[1]*parms[2,1]*(1-Y[2]) - parms[2,2]*Y[2]
  dZp_dt= parms[4,1]*Y[4]  #0-order phosphorylation rate
  dWp_dt= 0 #constant kinase
  dUp_dt= parms[6,1]*Y[4] - parms[6,2]*Y[5] #0-order phosphorylation and first-order dephosphorylation
  return(list(c(dAp_dt, dBp_dt, dZp_dt, dWp_dt, dUp_dt))) 
}

par(mfrow=c(1,2))
Y_out= data.frame(ode(func= dYp_dt,y=Y0,times=0:100,parms=parmsX))
plot(Y_out$t,Y_out$Ap,ylim=c(0,1),main='A')
plot(Y_out$t,Y_out$Bp,ylim=c(0,1),main='B')
plot(Y_out$t,Y_out$Zp,ylim=c(0,1),main='Z')
lines(Y_out$t,Y_out$Wp,col='blue',lty=2)
plot(Y_out$t,Y_out$Up,ylim=c(0,1),main='U')
lines(Y_out$t,Y_out$Wp,col='blue',lty=2)
#Y_out -> columns time, concA, concB




##############################
#model k_XZ (NO DEPHOSPHORYLATION) - using JAGS
#Prior: assume gamma-distribution with low mean and low effective sample size for k_X

# As random variable, use the difference between subsequent timepoints
 #the mu-parameter of the normal distribution is k1*deltat 

#prepare data
deltaZp= diff(Y_out$Zp)
deltat= diff(Y_out$t)
Ntotal= length(deltaZp)
dataList= list(
  y=deltaZp,
  Ntotal= length(deltaZp)
)

#prepare model
modelString= "
 model {
  for (i in Ntotal) {
    y[i] ~ dnorm(k1dt, sigma)  #likelihood
  } 
  k1dt ~ dnorm(0,0.05) #prior
  sigma ~ dgamma(1,2)  #mean: a/b; sd= sqrt(a)/b
 }
"
writeLines(modelString, con="TEMPmodel.txt")


initList= function()  { #initial values generated via function: mu centered around 0; sigma in range 0-0.1
  initMu= rnorm(1,0,0.1)
  initSigma= runif(1,0,0.1)
  return( list( k1dt= initMu, sigma=initSigma))
}

jagsModel= jags.model( file= "TEMPmodel.txt", data=dataList, inits=initList,
                       n.chains=3, n.adapt=500)  #find MCMC sampler
update( jagsModel,n.iter=500)  #burn in

codaSamples= coda.samples(jagsModel, variable.names=c('k1dt','sigma'),n.iter=3334)

par(mfrow= c(1,1))
plotPost(codaSamples[,'k1dt'],main='k1dt',xlab= 'k1 * delta_t')
  #true parameter should be k_WZ * Wp= 0.01

diagMCMC(codaObject = codaSamples, parName = 'k1dt')

plotPost(codaSamples[,'sigma'],main='sigma',xlab= 'sigma')
diagMCMC(codaObject = codaSamples, parName = 'sigma')

plot(type='n',x=0,y=0,xlim=c(1,dim(codaSamples[[1]])[1]),xlab='iteration',ylab='k1 * deltat' )
cols= c('black','blue','green')
for (i in 1:length(codaSamples)){
  lines(codaSamples[[i]][,1], col= cols[i])
}
# plot(type='n',x=0,y=0,xlim=c(1,dim(codaSamples[[1]])[1]),xlab='iteration',ylab='sigma' )
# cols= c('black','blue','green')
# for (i in 1:length(codaSamples)){
#   lines(codaSamples[[i]][,2], col= cols[i])
# }
# plot(type='n',x=0,y=0,xlab='k1 * deltat', ylab='sigma' )
# for (i in 1:length(codaSamples)){
#   points(codaSamples[[i]][,1],codaSamples[[i]][,2])
# }









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
PriLikPost_modes= matrix(NA,nrow=dim(Y_out)[1],ncol=3)
PriLikPost_modes[1,3]= mean(k_XZ_prior) #initial prior mean
j=4
for(i in 2:dim(Y_out)[1]) {
  delta_t= Y_out[i,1]-Y_out[i-1,1]
  
  #Plot showing parameter-choices within data; different from likelihood, which is data within distributions based on parameters
  plot(x= seq(0,3,.1),dnorm(x= seq(0,3,.1),mean=Y_out[i,j],sd=0.2),xlab='c(Z)',main='conc.') #distribution around observed value
  abline(v=c(Y_out[i-1,j]+ran[1]*Y_out[i-1,5]*delta_t,   #range covered by parameter range
             Y_out[i-1,j]+ran[parRes/2]*Y_out[i-1,5]*delta_t,
             Y_out[i-1,j]+ran[parRes]*Y_out[i-1,5]*delta_t),
         col='blue',lty=2)
  
  plot(ran,k_XZ_prior,type='l',main='Prior k_XZ') #,xlim=c(0,1)
  
  #Likelihood: Normal distribution around calculated value with arbitrary SD 0.02
  #x= observed value; mean= Calculated value for each value on parameter range
  #k_XZ_lik= dnorm(x= Y_out[i-1,j]+ran*Y_out[i-1,5]*delta_t,mean=Y_out[i,j],sd=0.02)
  k_XZ_lik= dnorm(x= Y_out[i,j],mean=Y_out[i-1,j]+ran*Y_out[i-1,5]*delta_t,sd=0.02)
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
k_WU_LikPost_modes= matrix(NA,nrow=dim(Y_out)[1],ncol=2)
k_WU_LikPost_modes[1,2]= ran1[which(k_WU_prior==max( k_WU_prior))[1]]
k_dU_LikPost_modes= matrix(NA,nrow=dim(Y_out)[1],ncol=2)
k_dU_LikPost_modes[1,2]= ran2[which(k_dU_prior==max( k_dU_prior))[1]]

j=6  #colNo of U in Y_out
#ani.record(reset = TRUE) #reset
for(i in 2:dim(Y_out)[1]) {
  delta_t= Y_out[i,1]-Y_out[i-1,1]
  
  #Plot showing parameter-choices within data; different from likelihood, which is data within distributions based on parameters
  plot(x= seq(0,3,.1),dnorm(x= seq(0,3,.1),mean=Y_out[i,j],sd=0.2),xlab='c(U)',main='conc.') #distribution around observed value
  abline(v= c(Y_out[i-1,j]+ ran1[1]*Y_out[i-1,5]*delta_t - ran2[1]*Y_out[i-1,j]*delta_t,  #kinase parameter range at low PPTase 
              Y_out[i-1,j]+ ran1[parRes]*Y_out[i-1,5]*delta_t - ran2[1]*Y_out[i-1,j]*delta_t), col='blue',lty=2)
  abline(v= c(Y_out[i-1,j]+ ran1[1]*Y_out[i-1,5]*delta_t - ran2[parRes]*Y_out[i-1,j]*delta_t,  #kinase parameter range at high PPTase
              Y_out[i-1,j]+ ran1[parRes]*Y_out[i-1,5]*delta_t - ran2[parRes]*Y_out[i-1,j]*delta_t), col='red',lty=2)
  
  contour(x= ran1, y=ran2, z=prior_2D, main= 'Prior',xlab='k_WU',ylab='k_dU')
  
  #Likelihood: Normal distribution around calculated value with arbitrary SD 0.02
  #x= observed value; mean= Calculated value for each value on parameter range
  lik_2D= matrix(NA,nrow=parRes, ncol=parRes)
  for(a1 in 1:parRes){
    for(a2 in 1:parRes){
      lik_2D[a1,a2]= dnorm(x= Y_out[i,j],
                           Y_out[i-1,j]+ ran1[a1]*Y_out[i-1,5]*delta_t - ran2[a2]*Y_out[i-1,j]*delta_t,
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
plot(type='n',x=0,y=0,xlim=c(0,dim(Y_out)[1]),ylim=c(0,max(k_WU_LikPost_modes,na.rm=TRUE)),
     xlab='Iteration',ylab='k_WU Marg. Mode')
lines(1:dim(Y_out)[1], k_WU_LikPost_modes[,1])
lines(1:dim(Y_out)[1], k_WU_LikPost_modes[,2],col='blue')
plot(type='n',x=0,y=0,xlim=c(0,dim(Y_out)[1]),ylim=c(0,max(k_dU_LikPost_modes,na.rm=TRUE)),
     xlab='Iteration',ylab='k_dU Marg. Mode')
lines(1:dim(Y_out)[1], k_dU_LikPost_modes[,1])
lines(1:dim(Y_out)[1], k_dU_LikPost_modes[,2],col='blue')


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
k_XA_LikPost_modes= matrix(NA,nrow=dim(Y_out)[1],ncol=2)
k_XA_LikPost_modes[1,2]= ran1[which(k_XA_prior==max( k_XA_prior))[1]]
k_dA_LikPost_modes= matrix(NA,nrow=dim(Y_out)[1],ncol=2)
k_dA_LikPost_modes[1,2]= ran2[which(k_dA_prior==max( k_dA_prior))[1]]

j=2  #colNo of U in Y_out
l1= layout(mat=matrix(c(0,1,2,3),nrow=2,ncol=2,byrow=TRUE),heights=c(3,3),widths=c(3,3))
ani.record(reset = TRUE) #reset
for(i in 2:dim(Y_out)[1]) {
  delta_t= Y_out[i,1]-Y_out[i-1,1]
  
  #Plot showing parameter-choices within data; different from likelihood, which is data within distributions based on parameters
  #showing lowest and highest allowed parameter
  plot(x= seq(0,3,.1),dnorm(x= seq(0,3,.1),mean=Y_out[i,j],sd=0.2),xlab='c(A)',main='conc.') #distribution around observed value
  abline(v= c(Y_out[i-1,j]+ ran1[1]*(1-Y_out[i-1,j])*delta_t      - ran2[1]*Y_out[i-1,j]*delta_t,  #kinase parameter range at low PPTase 
              Y_out[i-1,j]+ ran1[parRes]*(1-Y_out[i-1,j])*delta_t - ran2[1]*Y_out[i-1,j]*delta_t), col='blue',lty=2)
  abline(v= c(Y_out[i-1,j]+ ran1[1]*(1-Y_out[i-1,j])*delta_t      - ran2[parRes]*Y_out[i-1,j]*delta_t,  #kinase parameter range at high PPTase
              Y_out[i-1,j]+ ran1[parRes]*(1-Y_out[i-1,j])*delta_t - ran2[parRes]*Y_out[i-1,j]*delta_t), col='red',lty=2)
  
  # contour(x= ran1, y=ran2, z=prior_2D, main= 'Prior',xlab='k_XA',ylab='k_dA',
  #         xlim=c(0,1),ylim=c(0,1))
  # persp(x= ran1, y=ran2, z=prior_2D, main= 'Prior',xlab='k_XA',ylab='k_dA',
  #         xlim=c(0,1),ylim=c(0,1))
  
  #Likelihood: Normal distribution around calculated value with arbitrary SD 0.02
  #x= observed value; mean= Calculated value for each value on parameter range
  lik_2D= matrix(NA,nrow=parRes, ncol=parRes)
  for(a1 in 1:parRes){
    for(a2 in 1:parRes){
      lik_2D[a1,a2]= dnorm(x= Y_out[i,j],
                           Y_out[i-1,j]+ ran1[a1]*(1-Y_out[i-1,j])*delta_t - ran2[a2]*Y_out[i-1,j]*delta_t,
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
plot(type='n',x=0,y=0,xlim=c(0,dim(Y_out)[1]),ylim=c(0,max(k_XA_LikPost_modes,na.rm=TRUE)),
     xlab='Iteration',ylab='k_XA Marg. Mode')
lines(1:dim(Y_out)[1], k_XA_LikPost_modes[,1])
lines(1:dim(Y_out)[1], k_XA_LikPost_modes[,2],col='blue')
plot(type='n',x=0,y=0,xlim=c(0,dim(Y_out)[1]),ylim=c(0,max(k_dA_LikPost_modes,na.rm=TRUE)),
     xlab='Iteration',ylab='k_dA Marg. Mode')
lines(1:dim(Y_out)[1], k_dA_LikPost_modes[,1])
lines(1:dim(Y_out)[1], k_dA_LikPost_modes[,2],col='blue')

#############################################
#add noise
t_vec= Y_out[,1]
Y_concs= as.matrix(Y_out[,c(2,3)])  #needs to be matrix
#Y_concs= matrix(1:60,nrow=30,ncol=2,byrow = TRUE)
set.seed(123)
sigm= 0.05
Y_3D= array(rep(NA,dim(Y_concs)[1]*2*3), dim=c(dim(Y_concs)[1],2,3))   #data for 3 replicates; replicates on third dimension
#Y_3D[,,1]= Y_concs

for (i in 1:3){
  Y_curr= cbind.data.frame(
              A= Y_concs[,1]+ rnorm(dim(Y_out)[1],mean=0,sd=sigm),
              B= Y_concs[,2]+ rnorm(dim(Y_out)[1],mean=0,sd=sigm) )
  Y_3D[,,i]= as.matrix(Y_curr)
}
p1<- plot(type='n',x=0,y=0,xlim=c(0,max(t_vec)),ylim=c(0,1))
p2<- plot(type='n',x=0,y=0,xlim=c(0,max(t_vec)),ylim=c(0,1))
for (i in 1:3){
  p1<- p1+points(t_vec,Y_3D[,1,i])
  p2<- p2+points(t_vec,Y_3D[,1,i])
}

k= 1   #k to explore
head(Y_out)
deltaYB_t4= Y_out[Y_out$time ==4,2] - Y_out[Y_out$time ==3,2]   #change in B btw TPs
YA_t3= Y_out[Y_out$time==3,1]


deltaYB_t4= deltaYB_t4+ rnorm(3,mean=0,sd=sigm) 
YA_t3= YA_t3+ rnorm(3,mean=0,sd= k*sigm)

#probability of observed data under model
P= prod( dnorm(deltaYB_t4,mean= k*YA_t3,sd=sigm) )

#20230502: one vectorized instead of multiple scalar ODE functions
