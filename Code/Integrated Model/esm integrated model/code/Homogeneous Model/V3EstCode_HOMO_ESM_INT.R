#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#   Error scaling through tracking model homogeneous
#                 Estimation Code
#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


library(gtools)
library(bayesm)

setwd("/Users/rogerbailey/Dropbox/ESM Integrated Model/Code")

set.seed(77)

#Load data
load("Simdata_HOMO_ESM_Int.RData")

#Set MCMC variables
R=10000 #number of iterations
keep=10 #thinning
accup=100#interval size for updating RW step size
space=10

#Construct stacked design matrix
Xstacked= NULL
for(i in 1:nresp){
  Xstacked=rbind(Xstacked,data[[i]]$X)
}

#construct stacked tracking data
Trstacked=NULL
for(i in 1:nresp){
  Trstacked=rbind(Trstacked,data[[i]]$Tr)
}

#Set priors
betabar=matrix(double(nbeta),nr=1)#mean of homo partworths
Abeta=matrix(.01)#prec for homo partworths
Vbeta=diag(nbeta)#loc for cov of homo partworths
Vbetai=backsolve(chol(Vbeta),diag(nbeta))#inv chol Vbeta
alphabar=matrix(0)#mean of engagement slope
Aalpha=matrix(.01)#prec for engagement slope
Valpha=matrix(1)#loc for var of engagement
Valphai=backsolve(chol(Valpha),diag(1))#inv chol Vbeta
ATr=diag(2)*.01#prec for tracking variables
nuTr=nTr*2+3#df for cov matrix of tracking variables
VTr=nuTr*diag(nTr)#loc for cov matrix  of tracking variables
lambdabar=matrix(double(nTr*2),nc=3)  #mean of variables that transform engagment into tracking

#Set initial values
oldbeta=matrix(double(nbeta),nc=1)#initial beta
oldalpha=talpha#initial engagement multiplier
stepsizeb=.02 #beta stepsize
acceptpropb=double(accup)+.23 #proportion of beta draws accepted
stepsizea=.005
oldVTr=.5*diag(nTr*2)#initial VTr
oldVTri=backsolve(chol(oldVTr),diag(nTr*2))#initial inv chol VTr
acceptpropa=double(accup)+.23 #proportion of alpha draws accepted
llike=matrix(double(nresp),nc=1) #log likelihood
oldbetasigmat=t(oldbeta%*%(1+matrix(rep(rep(c(0:(ntask-1)),each=nalt),nresp),nr=1)*oldalpha)^(-1))
olduts=(Xstacked*oldbetasigmat)%*%matrix(double(nbeta)+1,nc=1)
oldlambda=matrix(double(nTr*2),nr=2)

#Setup functions
#function that returns the value of the prior for betas
logprior=function(beta,betabar,Vbi){
  return((beta-betabar)%*%Vbi%*%t(beta-betabar)*(-.5))
}

#loglike function
loglike=function(xb,y){
  probs=log(exp(xb))-log(matrix(rep(colSums(exp(xb)),length(
    xb[,1])),byrow=T,nc=length(xb[1,])))
  loc=cbind(y,c(1:(ncol(xb))))
  return(sum(probs[loc]))
}
#function for determining the change in the step size
stepupdate=function(accprop){
  step=1
  if(is.na(accprop)){return(step)}else{
    if(accprop<.21) {step=.99}
    if(accprop<.19) {step=.95}
    if(accprop<.15) {step=.85}
    if(accprop<.10) {step=.7}
    if(accprop>.25) {step=1.01}
    if(accprop>.27) {step=1.05}
    if(accprop>.3) {step=1.15}
    if(accprop>.4) {step=1.35}
    return(step)}
}

#Setup storage
betadraws=matrix(double(nbeta*R/keep),nr=R/keep)
acceptpropbs=matrix(double(R/keep),nr=R/keep)
stepsizebdraws=matrix(double(R/keep),nr=R/keep)
alphadraws=matrix(double(R/keep),nr=R/keep)
acceptpropas=matrix(double(R/keep),nr=R/keep)
stepsizeadraws=matrix(double(R/keep),nr=R/keep)
lambdadraws=matrix(double(nTr*2*R/keep), nr=R/keep)
llikes=matrix(double(R/keep))

#timer
itime = proc.time()[3]

for(r in 1:R){
  #current engagement vector 
  oldeng=c(0:(ntask-1))*oldalpha
  #Draw new betas and calculate likelihood
  oldbetasigmat=t(oldbeta%*%(1+matrix(rep(rep(oldeng,each=nalt),nresp),nr=1))^(-1))
  olduts=matrix((Xstacked*oldbetasigmat)%*%matrix(double(nbeta)+1,nc=1),nr=nalt)
  newbeta=oldbeta+rnorm(nbeta)*stepsizeb
  newbetasigmat=t(newbeta%*%(1+matrix(rep(rep(oldeng,each=nalt),nresp),nr=1))^(-1))
  newuts=matrix((Xstacked*newbetasigmat)%*%matrix(double(nbeta)+1,nc=1),nr=nalt)
  oldllikeb=loglike(olduts,ycomb)
  newllikeb=loglike(newuts,ycomb)
  oldlpriorb=logprior(t(oldbeta),betabar,Vbetai)
  newlpriorb=logprior(t(newbeta),betabar,Vbetai) 
  diffvecb=newllikeb+newlpriorb-(oldllikeb+oldlpriorb)
  alphab=min(exp(diffvecb), 1)
  
  #accept or reject new beta
  draw=runif(1)
  acceptb=0
  if(alphab>draw){acceptb=1}
  if(acceptb==1){
    oldbeta=newbeta
    oldbetasigmat= newbetasigmat
    olduts=newuts
  }
  
  #Draw new lambdas regressing tracking on current alpha
  Xlambda=cbind(double(nresp*ntask)+1,rep(oldeng,nresp))
  outlambda=rmultireg(Trstacked,Xlambda,lambdabar,ATr,nuTr,VTr)
  oldlambda=matrix(outlambda$B,nc=nTr)
  oldVlambda=outlambda$Sigma
  oldVlambdai=backsolve(chol(oldVlambda), diag(nTr*2))
  
  #Draw new alpha
  newalpha=oldalpha+rnorm(1)*stepsizea
  while(newalpha<0){
    newalpha=oldalpha+rnorm(1)*stepsizea
  }
  #calculate likelihood of choices given new and old alpha
  #setup for likelihood of choices by calculating utilities
  oldbetasigmat=t(oldbeta%*%(1+matrix(rep(rep(oldeng,each=nalt),nresp),nr=1))^(-1))
  oldutsa=matrix((Xstacked*oldbetasigmat)%*%matrix(double(nbeta)+1,nc=1),nr=nalt)
  neweng=c(0:(ntask-1))*newalpha
  newbetasigmat=t(oldbeta%*%(1+matrix(rep(rep(neweng,each=nalt),nresp),nr=1))^(-1))
  newutsa=matrix((Xstacked*newbetasigmat)%*%matrix(double(nbeta)+1,nc=1),nr=nalt)
  #get likelihood of tracking given old and new alphas
  tracklikeold=0
  tracklikenew=0
  for(i in 1:nresp){
    #calculate density value for tracking variables given alpha and lambda
    muvecold=cbind(double(ntask)+1,1+oldeng)%*%oldlambda
    muvecnew=cbind(double(ntask)+1,1+neweng)%*%oldlambda
    for(j in 1:ntask){
      tracklikeold=tracklikeold+lndMvn(Tr[i,j,],muvecold[j],oldVlambdai)
      tracklikenew=tracklikenew+lndMvn(Tr[i,j,],muvecnew[j],oldVlambdai)
    }
  }
  
  #get likelihood and prior values
  oldllikea=loglike(oldutsa,ycomb)+tracklikeold
  newllikea=loglike(newutsa,ycomb)+tracklikenew
  oldlpriora=logprior(oldalpha,alphabar,Valphai)
  newlpriora=logprior(newalpha,alphabar,Valphai)
  diffveca=newllikea+newlpriora-(oldllikea+oldlpriora)
  alphaa=min(exp(diffveca), 1)
  
  #accept or reject
  draw=runif(1)
  accepta=0
  if(alphaa>draw){accepta=1}
  llike=loglike(oldutsa,ycomb)
  if(acceptb==1){
    oldalpha=newalpha
    llike=loglike(newutsa,ycomb)
  }
  
  #Store acceptance proportions for beta and alpha draws
  acceptpropb=c(acceptpropb[2:accup],acceptb)
  acceptpropa=c(acceptpropa[2:accup],accepta)
  
  #Store values
  if(r%%keep==0){
    betadraws[r/keep,]=oldbeta
    alphadraws[r/keep]=oldalpha
    lambdadraws[r/keep,]=oldlambda
    acceptpropbs[r/keep]=mean(acceptpropb)
    acceptpropas[r/keep]=mean(acceptpropa)
    llikes[r/keep]=llike
    stepsizebdraws[r/keep]=stepsizeb
    stepsizeadraws[r/keep]=stepsizea
  }
  
  #print progress
  #print tte and chart current draw progress
  if(r%%(keep*space)==0){
    par(mfrow=c(3,1))
    ctime = proc.time()[3]
    tuntilend = ((ctime - itime)/r) * (R + 1 - r)
    cat(" ", r, " (", round(tuntilend/60, 1), ")", 
        fill = TRUE)
    plot(llikes,type="l",ylab="Log Likelihood")
    matplot(betadraws,type="l",ylab="Betabar Draws")
    matplot(lambdadraws,type="l",ylab="Lambda Draws")
    fsh()
  }
  
  #update stepsizes
  if(r%%accup==0&&r<(.3*R)){
    stepsizeb=stepsizeb*stepupdate(mean(acceptpropb))
    stepsizea=stepsizea*stepupdate(mean(acceptpropa))
  }
  
}  

plot.bayesm.mat(betadraws,tvalue=tbeta,burnin=200)
plot.bayesm.mat(lambdadraws,tvalue=tOmega,burnin=200)
plot.bayesm.mat(alphaadraws,tvalue=talpha,burnin=200)

