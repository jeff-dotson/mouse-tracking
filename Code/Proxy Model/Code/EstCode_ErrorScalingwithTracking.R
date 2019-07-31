#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#   Error Scaling with Mousetracking  - Estimation Code
#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
library(gtools)
library(bayesm)

#Set working drive
setwd("/Users/rogerbailey/Desktop/Projects/ESM/Code")

set.seed(77)

#Load data
load("Simdata_BasicErrScaling.RData")


#Set MCMC variables
R=10000 #number of iterations
keep=5#thinning
accup=100#interval size for updating RW step size
space=10

#number of betas
nbeta=natt*nlvl+1


#Set priors
nubeta=nbeta+3
Vb=nubeta*diag(nbeta)
Abeta=matrix(.01)
betadoublebar=matrix(double(nbeta),nc=1)
Omegabar=double(nM)
AOmega=1

#Set initial values
oldbetabar=matrix(double(nbeta),nc=1)#initial betabar
oldbetamat=matrix(double(nbeta*nresp),nr=nresp)#initial betas
stepsizeb=.2 #beta stepsize
stepsizeO=.025 #Omega stepsize
oldVbeta=.5*diag(nbeta)#initial Vbeta
oldVbetai=backsolve(chol(oldVbeta),diag(nbeta))#initial inv chol Vbeta
oldOmega=tOmega#matrix(double(nM),nr=1)#Initial omega vector
oldVOmega=.01*diag(nM)
oldVOmegai=backsolve(chol(oldVOmega),diag(nM))
acceptpropb=double(accup)+.23 #proportion of beta draws accepted
acceptpropO=double(accup)+.23 #proportion of Omega draws accepted
llike=matrix(double(nresp),nc=1) #log likelihood




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
betadraws=array(double(nbeta*nresp*R/keep),dim=c(R/keep,nresp,nbeta))
betabardraws=matrix(double(nbeta*R/keep),nr=R/keep)
Omegadraws=matrix(double(nM*R/keep),nr=R/keep)
llikes=matrix(double(R/keep),nr=R/keep)
acceptpropbs=matrix(double(R/keep),nr=R/keep)
acceptpropOs=matrix(double(R/keep),nr=R/keep)
stepsizebdraws=matrix(double(R/keep),nr=R/keep)
stepsizeOdraws=matrix(double(R/keep),nr=R/keep)


#timer
itime = proc.time()[3]

for(r in 1:R){
  accept=matrix(double(nresp),nc=nresp)#iter accept storage
  oldgammas=oldOmega%*%M
  olduts=NULL
  #Draw new betas and screen cutoffs for each respondent
  for(i in 1:nresp){
    #get new and old paramter vectors for RW
    oldbeta=matrix(oldbetamat[i,],nc=1)
    newbeta=oldbeta+t(chol(oldVbeta))%*%rnorm(nbeta)*stepsizeb
    scaler=matrix((1+oldgammas[i]*(c(0:(ntask-1))%x%rep(1,nalt)))^(-1),nr=nalt)
    oldutmat=matrix(data[[i]]$X%*%oldbeta,nr=nalt)
    newutmat=matrix(data[[i]]$X%*%newbeta,nr=nalt)
    #get likelihood and prior values
    oldllikeb=loglike(oldutmat*scaler,as.matrix(data[[i]]$y))
    newllikeb=loglike(newutmat*scaler,as.matrix(data[[i]]$y))
    oldlpriorb=logprior(t(oldbeta),t(oldbetabar),oldVbetai)
    newlpriorb=logprior(t(newbeta),t(oldbetabar),oldVbetai)
    diffvecb=newllikeb+newlpriorb-(oldllikeb+oldlpriorb)
    alphab=min(exp(diffvecb), 1)
    
    #accept or reject
    draw=runif(1)
    acceptb=0
    if(alphab>draw){acceptb=1}
    accept[i]=acceptb
    if(acceptb==1){
      oldbetamat[i,]=newbeta
      oldutmat=newutmat
    }
    olduts=cbind(olduts,oldutmat)#save unscaled utilities for later use
  }
  
  #Draw new Omega
  #create new omega matrix
  newOmega=oldOmega#+rnorm(nM)*stepsizeO
  newgammas=newOmega%*%M
  while((min(newgammas))*ntask<(-1)){
    newOmega=oldOmega+rnorm(nM)*stepsizeO
    newgammas=newOmega%*%M
  }
  scalerold=matrix((1+oldgammas%x%(c(0:(ntask-1))%x%rep(1,nalt)))^(-1),nr=nalt)
  scalernew=matrix((1+newgammas%x%(c(0:(ntask-1))%x%rep(1,nalt)))^(-1),nr=nalt)
  
  #get likelihood and prior values
  oldllikeO=loglike(olduts*scalerold,ycomb)
  newllikeO=loglike(olduts*scalernew,ycomb)
  oldlpriorO=logprior(oldOmega,Omegabar,oldVOmegai)
  newlpriorO=logprior(newOmega,Omegabar,oldVOmegai)
  diffvecO=newllikeO+newlpriorO-(oldllikeO+oldlpriorO)
  alphaO=min(exp(diffvecO), 1)
  
  #accept or reject
  draw=runif(1)
  acceptO=0
  if(alphaO>draw){acceptO=1}
  llike=oldllikeO
  if(acceptb==1){
    oldOmega=newOmega
    llike=newllikeO
  }
  
  
  #Store acceptance proportions
  acceptpropb=c(acceptpropb[2:accup],mean(accept))
  acceptpropO=c(acceptpropb[2:accup],acceptO)
  
  #Draw new values for betabar
  outbetaup=rmultireg(oldbetamat,matrix(1,nr=nresp,nc=1),matrix(betadoublebar,nr=1),Abeta,nubeta,Vb)
  oldbetabar=matrix(outbetaup$B)
  oldVbeta=outbetaup$Sigma
  oldVbetai=chol2inv(chol(oldVbeta))
  
  
  #Store values
  if(r%%keep==0){
    betadraws[r/keep,,]=oldbetamat
    betabardraws[r/keep,]=oldbetabar
    Omegadraws[r/keep,]=oldOmega
    acceptpropbs[r/keep]=mean(acceptpropb)
    acceptpropOs[r/keep]=mean(acceptpropO)
    llikes[r/keep]=llike
    stepsizebdraws[r/keep]=stepsizeb
    stepsizeOdraws[r/keep]=stepsizeO
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
    matplot(betabardraws,type="l",ylab="Betabar Draws")
    matplot(Omegadraws,type="l",ylab="Omega Draws")
    fsh()
  }
  
  #update stepsizes
  if(r%%accup==0&&r<(.3*R)){
    stepsizeb=stepsizeb*stepupdate(mean(acceptpropb))
    stepsizeO=stepsizeO*stepupdate(mean(acceptpropO))
  }
  
}  



