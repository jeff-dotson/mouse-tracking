#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#   Poisson race tracking and response time- Estimation Code
#   Roger Bailey
#   May 2016
#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

library(gtools)
library(bayesm)

 
setwd("/Users/rogerbailey/Dropbox/Projects/RaceTracking/Tracking Model Code")
 
set.seed(77)

#Load data
load("Simdata_PoisRaceTrackingV6_Balanced.RData")

#Set MCMC variables
R=10000 #number of iterations
keep=2#thinning
accup=40 #interval size for updating RW step size
space=10

#Set priors for beta
nub=nbeta+3#dof for cov matrix of partworths
Vb=nub*diag(nbeta)#loc for cov of partworths
Ab=matrix(.01)#prec for betabar
bdoublebar=matrix(double(nbeta),nc=1)#mean of betabar distribution
 
#Set priors for gammas
Apg=diag(nTr+1)*.01#prec for deltalambdabar
psigammabar=matrix(c(1,double(nTr)),nc=1)#mean of distribution of gammas


#Get large matirx of tracking variables
Omegamat=NULL
for(i in 1:nresp){
  Omegamat=rbind(Omegamat,data[[i]]$Om)
}
#Set initial values
oldbbar=tbetabar#matrix(double(nbeta),nc=1)#initial betabar
oldbmat=tbetas#matrix(double(nbeta*nresp),nr=nresp)#initial betas
stepb=.15 #stepsize for beta RW steps
steppg=.25  #stepsize for gamma RW step
oldVb=diag(nbeta)#initial loc for cov of partworths
oldVbi=backsolve(chol(oldVb),diag(nbeta))#initial inv chol Vbeta
oldKmat=tK#matrix(double(nresp)+1,nr=nresp)#initial hit thresholds
oldpsigamma=c(tpsi,tgamma)#matrix(c(1,double(nTr)))#initial value of gammas
acceptpropb=double(accup)+.23 #running record of proportion of beta,lambda,delta draws accepted
acceptpropdl=double(accup)+.23 #running record of proportion of beta,lambda,delta draws accepted
acceptpropK=double(accup)+.23 #running record of proportion of K draws accepted 
acceptproppg=double(accup)+.23 #running record of proportion of gamma draws accepted
llike=matrix(double(nresp),nc=1) #log likelihood

#Setup storage
betadraws=array(double(nbeta*nresp*R/keep),dim=c(R/keep,nresp,nbeta))
betabardraws=matrix(double(nbeta*R/keep),nr=R/keep)
Vbetadraws=matrix(double(nbeta*nbeta*R/keep),nr=R/keep)
Kdraws=matrix(double(nresp*ntask*R/keep),nr=R/keep)
psidraws=matrix(double(nresp*R/keep),nr=R/keep)
psibardraws=matrix(double(R/keep),nr=R/keep)
Vpsidraws=matrix(double(R/keep),nr=R/keep)
gammadraws=matrix(double(R/keep*nTr),nr=R/keep)
llikes=matrix(double(nresp*R/keep),nr=R/keep)
accpropb=matrix(double(R/keep),nr=R/keep)
accpropK=matrix(double(R/keep),nr=R/keep)
accproppg=matrix(double(R/keep),nr=R/keep)
stsizeb=matrix(double(R/keep),nr=R/keep)
stsizepg=matrix(double(R/keep),nr=R/keep)
thlike=matrix(double(R/keep),nr=R/keep)
 
#*************************************************************
#Setup functions
#function that returns the value of the prior for norm upper
logprior=function(beta,betabar,Vbi){
  return((beta-betabar)%*%Vbi%*%t(beta-betabar)*(-.5))
}
#loglike function for choice/response time
loglikeyt=function(bdl,K,X,y,t){
  Kvec=rep(K,each=nalt)
  baserates=exp(X%*%matrix(bdl[1:nbeta],nc=1))
  rates=baserates
  choiceind=matrix(diag(nalt)[,y])
  tvec=rep(t,each=nalt)
  ll=0
  for(m in 1:(ntask*nalt)){
    if(choiceind[m]==1){
      ll=ll+dgamma(tvec[m],shape=Kvec[m],scale=rates[m],log=T)
    }else{ll=ll+log(1-pgamma(tvec[m],shape=Kvec[m],scale=rates[m]))}
  }
  return(ll)
}
#logprior function for K given theta
logpriorK=function(K,th){
  return(sum(dpois(K-1,th,log=T)))
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
  return(step)
  }
}
 
#*************************************************************
#begin MCMC routine
   
#set timer
itime = proc.time()[3]

for(r in 1:R){
  thetalike=double(nresp)
  if(r/keep==100){
  zeros=double(nresp)
  zeros[llikevec==-Inf]=1
  oldbmat[zeros==1,]=tbetas[zeros==1,]
  }
  accept=matrix(double(nresp*4),nr=nresp)
  llikevec=matrix(double(nresp))
  for(i in 1:nresp){
    #draw proposal for betas
    oldbdl=c(oldbmat[i,])
    newbdl=oldbdl+c(t(chol(Vb))%*%rnorm(nbeta)*stepb)
    #calculate likelihood of choices/response times and priors
    oldllikeb=loglikeyt(oldbdl,oldKmat[i,],data[[i]]$X,data[[i]]$y,data[[i]]$time)
    newllikeb=loglikeyt(newbdl,oldKmat[i,],data[[i]]$X,data[[i]]$y,data[[i]]$time)
    oldlprb=logprior(t(oldbdl[1:nbeta]),t(oldbbar),oldVbi)
    newlprb=logprior(t(newbdl[1:nbeta]),t(oldbbar),oldVbi)
    diffvecb=newllikeb+newlprb-(oldllikeb+oldlprb)
    if(is.nan(diffvecb)){diffvecb=-Inf} #keep old value if new value creates nan
    alphab=min(exp(diffvecb), 1)
    #accept or reject new draw of beta
    drawb=runif(1)
    acceptb=0
    if(alphab>drawb){acceptb=1}
    accept[i,1]=acceptb
    if(acceptb==1){
      oldbmat[i,]=newbdl[1:nbeta]
      oldbdl=newbdl
    }
      
         
    #draw proposal for K
    oldK=oldKmat[i,]
    old1ind=double(ntask)
    old1ind[oldK==1]=1
    randpert=double(ntask)
    randpert[old1ind==0]=sample(c(-1,0,1),size=ntask-sum(old1ind), replace=T, prob=double(3)+1/3)
    randpert[old1ind==1]=sample(c(0,1),size=sum(old1ind), replace=T, prob=double(2)+1/2)
    newK=oldK+randpert

    #Get current theta
    baseth=data[[i]]$Om%*%oldgamma
    oldpsi=oldpmat[i]
    oldtheta=baseth+oldpsi
    #calculate likelihood of choices/response times and priors
    oldllikeK=loglikeyt(oldbdl,oldK,data[[i]]$X,data[[i]]$y,data[[i]]$time)
    newllikeK=loglikeyt(oldbdl,newK,data[[i]]$X,data[[i]]$y,data[[i]]$time)
    oldlprK=logpriorK(oldK,oldtheta)
    newlprK=logpriorK(newK,oldtheta)
    newtooldK=(1/2)^(sum(newK==1))*(1/3)^(ntask-sum(newK==1))
    oldtonewK=(1/2)^(sum(old1ind))*(1/3)^(ntask-sum(old1ind))
    diffvecK=newllikeK+newlprK+newtooldK-(oldllikeK+oldlprK+oldtonewK)
    if(is.nan(diffvecK)){diffvecK=-Inf}
    alphaK=min(exp(diffvecK), 1)
    #accept or reject new draw of beta,lambda and delta
    drawK=runif(1)
    acceptK=0
    if(alphaK>drawK){acceptK=1}
    accept[i,3]=acceptK
    llikevec[i]=oldllikeK
    if(acceptK==1){
      oldKmat[i,]=newK
      oldK=newK
      llikevec[i]=newllikeK
    }
  }

  #draw new psi and gamma
  aggKvec=matrix(t(oldKmat),nc=1)
  acceptpg=0
  newpsigamma=oldpsigamma+rnorm(ntr+1)*steppg
  oldaggthetavec=Omegamat%*%oldpsigamma
  newaggthetavec=Omegamat%*%newpsigamma
  if(min(newaggthetavec)>0){
    oldllikepg=logpriorK(aggKvec,oldaggthetavec)
    newllikepg=logpriorK(aggKvec,newaggthetavec)
    oldlprpg=logprior(t(oldgammapg),t(psigammabar),Apg)
    newlprpg=logprior(t(newgammapg),t(psigammabar),Apg)
    diffvecpg=newllikeg+newlprg+newtooldg-(oldllikeg+oldlprg+oldtonewg)
    if(is.nan(diffvecpg)){diffvecpg=-Inf}
    alphapg=min(exp(diffvecpg), 1)
    #accept or reject new draw of theta
    drawpg=runif(1)
    if(alphag>drawg){acceptpg=1}
  }
  if(acceptpg==1){
    oldpsigamma=newpsigamma
  }
    
  
  #draw new values of beta hyperparameters
  outbetahyp=rmultireg(oldbmat,matrix(1,nr=nresp,nc=1),matrix(bdoublebar,nr=1),Ab,nub,Vb)
  oldbbar=matrix(outbetahyp$B)
  oldVb=outbetahyp$Sigma
  oldVbi=chol2inv(chol(oldVb))
       
  #draw new values of psi hyperparameters
  outpsihyp=rmultireg(matrix(oldpmat),matrix(1,nr=nresp,nc=1),matrix(pdoublebar,nr=1),Ap,nup,Vp)
  oldpbar=matrix(outpsihyp$B)
  oldVp=outpsihyp$Sigma
  oldVpi=chol2inv(chol(oldVp))
       
  #Store acceptance proportions
  acceptpropb=c(acceptpropb[2:accup],mean(accept[,1]))
  acceptpropK=c(acceptpropK[2:accup],mean(accept[,3]))
  acceptpropp=c(acceptpropp[2:accup],mean(accept[,4]))
  acceptpropg=c(acceptpropg[2:accup],mean(acceptgvec))
       
  #Store values
  if(r%%keep==0){
    #Setup storage
    betadraws[r/keep,,]=oldbmat
    betabardraws[r/keep,]=oldbbar
    Vbetadraws[r/keep,]=oldVb
    Kdraws[r/keep,]=matrix(t(oldKmat),nc=1)
    psidraws[r/keep,]=matrix(oldpmat)
    psibardraws[r/keep,]=matrix(oldpbar)
    Vpsidraws[r/keep,]=matrix(oldVp)
    gammadraws[r/keep,]=oldgamma
    llikes[r/keep,]=llikevec
    accpropb[r/keep]=mean(acceptpropb)
    accpropK[r/keep]=mean(acceptpropK)
    accpropp[r/keep]=mean(acceptpropp)
    accpropg[r/keep]=mean(acceptpropg)
    stsizeb[r/keep]=stepb
    stsizep[r/keep]=stepp
    stsizeg[r/keep]=stepg
    thlike=sum(thetalike)
  }
       
  #print progress
  #print tte and chart current draw progress
  if(r%%(keep*space)==0){
    par(mfrow=c(4,1))
    ctime = proc.time()[3]
    tuntilend = ((ctime - itime)/r) * (R + 1 - r)
    test=double(nresp)
    test[llikevec==-Inf]=1
    te=sum(test)
    cat(" ", r, " (", round(tuntilend/60, 1), ")", sum(thetalike), fill = TRUE)
    plot(rowSums(llikes),type="l",ylab="Log Likelihood")
    matplot(cbind(betabardraws,matrix(rep(tbetabar,R/keep),byrow=T,nc=nbeta)),type="l",ylab="Betabar Draws")
    plot(psibardraws,type="l",ylab="Psibar Draws", col="red")
    abline(h=tpsibar,col="blue")
    matplot(cbind(gammadraws,matrix(rep(tgamma,R/keep),byrow=T,nc=nTr)),type="l",ylab="Gamma Draws")
    fsh()
  }
       
  #update stepsizes
  if(r%%accup==0&&r<(.3*R)){
    stepb=stepb*stepupdate(mean(acceptpropb))
    stepp= stepp*stepupdate(mean(acceptpropp))
    stepg= stepg*stepupdate(mean(acceptpropg))
  }
}
plot.bayesm.mat(betabardraws, burnin=0, tvalues=tbetabar)
plot.bayesm.mat(deltalambdabardraws, burnin=0, tvalues=tdeltalambdabar)
for(i in 1:20){
  plot.bayesm.mat(psidraws, burnin=0, tvalues=tpsi)
}
plot.bayesm.mat(gammadraws, burnin=0, tvalues=tgamma)


