#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#   Poisson race basic reponse time model - Estimation Code
#   February 2016
#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

library(gtools)
library(bayesm)

setwd("/Users/rogerbailey/Desktop/Projects/ESM/RaceTracking/Basic Model Code")

set.seed(77)

#Load data
load("Simdata_PoisRaceBasic.RData")

#Set MCMC variables
R=60000 #number of iterations
keep=20 #thinning
accup=50#interval size for updating RW step size
space=10

#Set priors for beta
nub=nbeta+3#dof for cov matrix of partworths
Vb=nub*diag(nbeta)#loc for cov of partworths
Ab=matrix(.01)#prec for betabar
bdoublebar=matrix(double(nbeta),nc=1)#mean of betabar distribution

#Set priors for delta and lambda
nudl=4+3#dof for cov matrix of deltas and lambdas
Vdl=nudl*diag(4)#loc for cov of deltas and lambdas
Adl=matrix(.01)#prec for deltalambdabar
dldoublebar=matrix(double(4),nc=1)#mean of deltalambdbar distribution

#Set priors for hit thresholds, K
lnthetabar=matrix(0)
Ath=matrix(.01)

#Set initial values
oldbbar=matrix(double(nbeta),nc=1)#initial betabar
oldbmat=matrix(double(nbeta*nresp),nr=nresp)#initial betas
olddlbar=matrix(double(4),nc=1)#initial deltalambdabar
olddlmat=matrix(double(4*nresp),nr=nresp)#initial deltas and lambdas
stepbdl=.2 #stepsize for beta, delta, and lambda RW steps
oldVb=diag(nbeta)#initial loc for cov of partworths
oldVbi=backsolve(chol(oldVb),diag(nbeta))#initial inv chol Vbeta
oldVdl=diag(4)#initial loc for cov of deltas and lambdas
oldVdli=backsolve(chol(oldVdl),diag(4))#initial inv chol Vdeltalambda
oldKmat=matrix(double(nresp)+1,nr=nresp)#initial hit thresholds
oldtheta=.01#initial rate for draw of hit thresholds
stepth=.1 #initial stepsize for theta RW
acceptpropbdl=double(accup)+.23 #running record of proportion of beta,lambda,delta draws accepted
acceptpropK=double(accup)+.23 #running record of proportion of K draws accepted
acceptpropth=double(accup)+.23 #running record of proportion of theta draws accepted
llike=matrix(double(nresp),nc=1) #log likelihood

#Setup storage
betadraws=array(double(nbeta*nresp*R/keep),dim=c(R/keep,nresp,nbeta))
betabardraws=matrix(double(nbeta*R/keep),nr=R/keep)
Vbetadraws=matrix(double(nbeta*nbeta*R/keep),nr=R/keep)
deltalambdadraws=array(double(4*nresp*R/keep),dim=c(R/keep,nresp,4))
deltalambdabardraws=matrix(double(4*R/keep),nr=R/keep)
Vdeltalambdadraws=matrix(double(16*R/keep),nr=R/keep)
Kdraws=matrix(double(nresp*R/keep),nr=R/keep)
thetadraws=matrix(double(R/keep),nr=R/keep)
llikes=matrix(double(nresp*R/keep),nr=R/keep)
accpropbdl=matrix(double(R/keep),nr=R/keep)
accpropK=matrix(double(R/keep),nr=R/keep)
accpropth=matrix(double(R/keep),nr=R/keep)
stsizebdl=matrix(double(R/keep),nr=R/keep)
stsizeth=matrix(double(R/keep),nr=R/keep)

#*************************************************************
#Setup functions
#function that returns the value of the prior for norm upper
logprior=function(beta,betabar,Vbi){
  return((beta-betabar)%*%Vbi%*%t(beta-betabar)*(-.5))
}
#loglike function for choice/response time
loglikeyt=function(bdl,K,X,y,t){
  baserates=exp(X%*%matrix(bdl[1:nbeta],nc=1))
  temp1=colSums(matrix(baserates,nr=nalt))^-1#inverted "attractiveness"
  temp2=exp(bdl[(nbeta+1)] + bdl[(nbeta+2)]*c(0:(ntask-1)))#"accessibility"
  temp3=exp(bdl[(nbeta+3)] + bdl[(nbeta+4)]*c(0:(ntask-1)))#scaling
  fvals=temp1*temp2+temp3
  rates=baserates*rep(fvals,each=nalt)
  choiceind=matrix(diag(nalt)[,y])
  tvec=rep(t,each=nalt)
  ll=0
  for(m in 1:(ntask*nalt)){
    if(choiceind[m]==1){
      ll=ll+dgamma(tvec[m],shape=K,scale=rates[m],log=T)
    }else{ll=ll+log(1-pgamma(tvec[m],shape=K,scale=rates[m]))}
  }
  return(ll)
}
#logprior function for K given theta
logpriorK=function(K,th){
  return(dpois(K-1,th,log=T))
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

#*************************************************************
#begin MCMC routine

#set timer
itime = proc.time()[3]

for(r in 1:R){
  accept=matrix(double(nresp*2),nr=nresp)
  llikevec=matrix(double(nresp))
  for(i in 1:nresp){
    #draw proposal for betas, deltas, and lambdas
    oldbdl=c(oldbmat[i,],olddlmat[i,])
    newbdl=c(oldbdl[1:nbeta]+t(chol(Vb))%*%rnorm(nbeta)*stepbdl,
          oldbdl[(nbeta+1):(nbeta+4)]+t(chol(Vdl))%*%rnorm(4)*stepbdl)
    #calculate likelihood of choices/response times and priors
    oldllikebdl=loglikeyt(oldbdl,oldKmat[i,],data[[i]]$X,data[[i]]$y,data[[i]]$time)
    newllikebdl=loglikeyt(newbdl,oldKmat[i,],data[[i]]$X,data[[i]]$y,data[[i]]$time)
    oldlprb=logprior(t(oldbdl[1:nbeta]),t(oldbbar),oldVbi)
    newlprb=logprior(t(newbdl[1:nbeta]),t(oldbbar),oldVbi)
    oldlprdl=logprior(t(oldbdl[(nbeta+1):(nbeta+4)]),t(olddlbar),oldVdli)
    newlprdl=logprior(t(newbdl[(nbeta+1):(nbeta+4)]),t(olddlbar),oldVdli)
    diffvecbdl=newllikebdl+newlprb+newlprdl-(oldllikebdl+oldlprb+newlprdl)
    if(is.nan(diffvecbdl)){diffvecbdl=-Inf}
    alphabdl=min(exp(diffvecbdl), 1)
    #accept or reject new draw of beta,lambda and delta
    drawbdl=runif(1)
    acceptbdl=0
    if(alphabdl>drawbdl){acceptbdl=1}
    accept[i,1]=acceptbdl
    if(acceptbdl==1){
      oldbmat[i,]=newbdl[1:nbeta]
      olddlmat[i,]=newbdl[(nbeta+1):(nbeta+4)]
      oldbdl=newbdl
    }
  
    #draw proposal for K
    oldK=oldKmat[i,]
    newK=oldKmat[i,]+(rbinom(1,1,.5)*2-1)
    #calculate likelihood of choices/response times and priors
    oldllikeK=loglikeyt(oldbdl,oldK,data[[i]]$X,data[[i]]$y,data[[i]]$time)
    if(newK>0){ #only consider proposals with K>0
    newllikeK=loglikeyt(oldbdl,newK,data[[i]]$X,data[[i]]$y,data[[i]]$time)
    oldlprK=logpriorK(oldK,oldtheta)
    newlprK=logpriorK(newK,oldtheta)
    diffvecK=newllikeK+newlprK-(oldllikeK+oldlprK)
    if(is.nan(diffvecK)){diffvecK=-Inf}
    alphaK=min(exp(diffvecK), 1)
    #accept or reject new draw of beta,lambda and delta
    drawK=runif(1)
    acceptK=0
    if(alphaK>drawK){acceptK=1}
    }else{acceptK=0}
    accept[i,2]=acceptK
    llikevec[i]=oldllikeK
    if(acceptK==1){
      oldKmat[i,]=newK
      llikevec[i]=newllikeK
    }
  }
  
  #draw new proposal for theta
  newtheta=exp(log(oldtheta)+rnorm(1)*stepth)
  #calculate lieklihood and prior for thetas(efficieny can be 
  #increased by doing this as part of the above respondent-level loop)
  newlliketh=0
  oldlliketh=0
  for(i in 1:nresp){
    oldlliketh=oldlliketh+logpriorK(oldKmat[i,],oldtheta)
    newlliketh=newlliketh+logpriorK(oldKmat[i,],newtheta)
  }
  oldlprth=logprior(matrix(log(oldtheta)),matrix(lnthetabar),Ath)
  newlprth=logprior(matrix(log(newtheta)),matrix(lnthetabar),Ath)
  diffvecth=newlliketh+newlprth-(oldlliketh+oldlprth)
  if(is.nan(diffvecth)){diffvecth=-Inf}
  alphath=min(exp(diffvecth), 1)
  #accept or reject new draw of theta
  drawth=runif(1)
  acceptth=0
  if(alphath>drawth){acceptth=1}
  if(acceptth==1){
  oldtheta=newtheta
  }
  
  #draw new values of beta hyperparameters
  outbetahyp=rmultireg(oldbmat,matrix(1,nr=nresp,nc=1),matrix(bdoublebar,nr=1),Ab,nub,Vb)
  oldbbar=matrix(outbetahyp$B)
  oldVb=outbetahyp$Sigma
  oldVbi=chol2inv(chol(oldVb))
  
  #draw new values of delta and lambda hyperparameters
  outdeltalambdahyp=rmultireg(olddlmat,matrix(1,nr=nresp,nc=1),matrix(dldoublebar,nr=1),Adl,nudl,Vdl)
  olddlbar=matrix(outdeltalambdahyp$B)
  oldVdl=outdeltalambdahyp$Sigma
  oldVdli=chol2inv(chol(oldVdl))
  
  #Store acceptance proportions
  acceptpropbdl=c(acceptpropbdl[2:accup],mean(accept[,1]))
  acceptpropK=c(acceptpropK[2:accup],accept[,2])
  acceptpropth=c(acceptpropth[2:accup],acceptth)
  
  #Store values
  if(r%%keep==0){
    #Setup storage
    betadraws[r/keep,,]=oldbmat
    betabardraws[r/keep,]=oldbbar
    Vbetadraws[r/keep,]=oldVb
    deltalambdadraws[r/keep,,]=olddlmat
    deltalambdabardraws[r/keep,]=olddlbar
    Vdeltalambdadraws[r/keep,]=oldVdl
    Kdraws[r/keep,]=oldKmat
    thetadraws[r/keep]=oldtheta
    llikes[r/keep,]=llikevec
    accpropbdl[r/keep]=mean(acceptpropbdl)
    accpropK[r/keep]=mean(acceptpropK)
    accpropth[r/keep]=mean(acceptpropth)
    stsizebdl[r/keep]=stepbdl
    stsizeth[r/keep]=stepth
  }
  
  #print progress
  #print tte and chart current draw progress
  if(r%%(keep*space)==0){
    par(mfrow=c(4,1))
    ctime = proc.time()[3]
    tuntilend = ((ctime - itime)/r) * (R + 1 - r)
    cat(" ", r, " (", round(tuntilend/60, 1), ")", 
        fill = TRUE)
    plot(rowSums(llikes),type="l",ylab="Log Likelihood")
    matplot(betabardraws,type="l",ylab="Betabar Draws")
    matplot(deltalambdabardraws,type="l",ylab="Delta Lambda Draws")
    plot(thetadraws,type="l", col="blue",ylab="Theta Draws")
    fsh()
  }
  
  #update stepsizes
  if(r%%accup==0&&r<(.3*R)){
    stepbdl=stepbdl*stepupdate(mean(acceptpropbdl))
    stepth=stepth*stepupdate(mean(acceptpropth))
  }
}
  
plot.bayesm.mat(betabardraws,tvalue=tbetabar,burnin=500)
plot.bayesm.mat(deltalambdabardraws,tvalue=tdeltalambdabar,burnin=500)
plot.bayesm.mat(thetadraws,tvalue=ttheta,burnin=500)
