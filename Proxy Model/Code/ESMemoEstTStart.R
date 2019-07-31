#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#   Error scaling through tracking/emotional recognition
#   Version 4
#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


library(gtools)
library(bayesm)
library(gdata)
set.seed(22)

ESMemoEstTStart=function(data,Prior,MCMC,sim,tvals){
  # FUNCTION OUTPUT:
  # betadraw - posterior draws of individual-level betas
  # baccept - acceptance rate for the beta step
  # bstepdraw - adjusted beta stepsize draws, fixed after burn-in
  # Vbetadraw - posterior draws of the upper level covariance matrix
  # Deltadraw - posterior draws of the upper level coefficient matrix
  # deltadraw - posterior draws of the upper level coefficients for resp-level error scale deviations 
  # omegadraw - posterior draws of the upper level coefficients for task-level error scale deviations 
  # oaccept - acceptance rate for the beta step
  # lambdadraw - posterior draws of individual-task-level error scale residuals
  # laccept - acceptance rate for the lambda step
  # lstepdraw - adjusted lambda stepsize draws, fixed after burn-in
  # Vlambdadraw - posterior draws of the upper level variance for the error scale 
  # llikedraw - draws of the log likelihood

  #get variables
  nresp=length(data)# number of respondents
  ntask=length(data[[1]]$y)# number of alternatives per task (inc outside)
  nalt=nrow(data[[1]]$X)/ntask#number of choice tasks
  nbeta=ncol(data[[1]]$X) #number of partworths to be estimated
  nZ=ncol(data[[1]]$Z)#number of upper level covariates(including unity)
  nTr=ncol(data[[1]]$Tr1)#number of tracking variables
  Z=data[[1]]$Z
  Tr1=data[[1]]$Tr1
  Tr2=data[[1]]$Tr2
  
  #get priors
  Deltabar=Prior$Deltabar
  ADelta=Prior$ADelta
  omegabar=Prior$omegabar
  Aomega=Prior$Aomega
  deltabar=Prior$deltabar
  Adelta=Prior$Adelta
  nub=Prior$nub
  Vb=Prior$Vb
  nul=Prior$nul
  Vl=Prior$Vl
  R=MCMC$R
  accav=MCMC$accav
  keep=MCMC$keep
  
  #* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  #functions
  #function returning the value of the prior for betas
  logprior=function(beta,betabar,Vbi){
    return((beta-betabar)%*%Vbi%*%t(beta-betabar)*(-.5))
  }
  #loglike function
  loglike=function(xb,y){
    probs=log(exp(xb))-log(matrix(rep(colSums(exp(xb)),length(
      xb[,1])),byrow=T,nc=length(xb[1,])))
    loc=cbind(y,c(1:(ncol(xb))))
    out=sum(probs[loc])
    if(is.nan(out)){out=-Inf}
    return(out)
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
  
  #Set initial values
  oldDelta=t(Tvals$tDelta)#matrix(double(nZ*nbeta),nc=nbeta)#initial betabar
  oldbetamat=Tvals$tbetas#matrix(double(nbeta*nresp),nr=nresp)#initial betas
  stepb=.5 #beta stepsize
  oldVbeta=Tvals$tSigmab#diag(nbeta)#initial Vbeta
  oldVbetai=backsolve(chol(oldVbeta),diag(nbeta))#initial inv chol Vbeta
  acceptpropb=double(accav)+.23 #proportion of beta draws accepted
  oldlambdavec=matrix(rowSums(Tvals$tscalecomp[,1:2]))#matrix(double(nresp),nr=nresp)#initial lambdas
  stepl=.7 #lambda stepsize
  oldVlambda=diag(1)#initial Vbeta
  oldVlambdai=backsolve(chol(oldVlambda),diag(1))#initial inv chol Vbeta
  acceptpropl=double(accav)+.23 #proportion of lambda draws accepted
  olddelta=Tvals$tdelta#matrix(double(nTr),nc=1) #initial delta
  oldomega=Tvals$tomega#matrix(double(nTr),nc=1) #initial omega
  stepo=.1 #omega stepsize
  acceptpropo=double(accav)+.23 #proportion of lambda draws accepted
  llike=matrix(double(nresp),nr=1) #log likelihood
  
  #Setup storage
  betadraw=array(double(nbeta*nresp*R/keep),dim=c(R/keep,nresp,nbeta))
  baccept=matrix(double(R/keep),nr=R/keep)
  bstepdraw=matrix(double(R/keep),nr=R/keep)
  Vbetadraw=matrix(double(nbeta^2*R/keep),nr=R/keep)
  Deltadraw=matrix(double(nZ*nbeta*R/keep),nr=R/keep)
  lambdadraw=matrix(double(R/keep*nresp),nr=R/keep)
  laccept=matrix(double(R/keep),nr=R/keep)
  lstepdraw=matrix(double(R/keep),nr=R/keep)
  Vlambdadraw=matrix(double(R/keep),nr=R/keep)
  deltadraw=matrix(double(R/keep*nTr),nr=R/keep)
  omegadraw=matrix(double(R/keep*nTr),nr=R/keep)
  oaccept=matrix(double(R/keep),nr=R/keep)
  ostepdraw=matrix(double(R/keep),nr=R/keep)
  llikes=matrix(double(R/keep),nr=R/keep)

  
  #start timing
  itime = proc.time()[3]
  for(r in 1:R){
    acceptmat=matrix(double(nresp*2),nr=nresp)#iter accept storage
    
    #Draw new betas for each respondent
    for(resp in 1:nresp){
      track=t(Tr2[((resp-1)*nTr+1):(resp*nTr),])
      oldlambda=oldlambdavec[resp]
      oldgamma=track%*%oldomega
      #get new and old paramter vectors for RW
      oldbeta=matrix(oldbetamat[resp,],nc=1)
      newbeta=oldbeta+t(chol(oldVbeta))%*%rnorm(nbeta)*stepb
      oldutmat=matrix(data[[resp]]$X%*%oldbeta/exp(rep(oldgamma+oldlambda,each=nalt)),nr=nalt)
      newutmat=matrix(data[[resp]]$X%*%newbeta/exp(rep(oldgamma+oldlambda,each=nalt)),nr=nalt)
      #get likelihood and prior values
      oldllike=loglike(oldutmat,as.matrix(data[[resp]]$y))
      newllike=loglike(newutmat,as.matrix(data[[resp]]$y))
      oldlprior=logprior(t(oldbeta),Z[resp,]%*%oldDelta,oldVbetai)
      newlprior=logprior(t(newbeta),Z[resp,]%*%oldDelta,oldVbetai)
      diffvec=newllike+newlprior-(oldllike+oldlprior)
      alpha=min(exp(diffvec), 1)
      #accept or reject
      draw=runif(1)
      accept=0
      if(alpha>draw){accept=1}
      acceptmat[resp,1]=accept
      if(accept==1){
        oldbetamat[resp,]=newbeta
        oldbeta=newbeta
      }

      
      #draw new lambda for each respondent
      newlambda=oldlambda+rnorm(1)*stepl
      #get new and old paramter vectors for RW
      oldutmat=matrix(data[[resp]]$X%*%oldbeta/exp(rep(oldgamma+oldlambda,each=nalt)),nr=nalt)
      newutmat=matrix(data[[resp]]$X%*%oldbeta/exp(rep(oldgamma+newlambda,each=nalt)),nr=nalt)
      #get likelihood and prior values
      oldllike=loglike(oldutmat,as.matrix(data[[resp]]$y))
      newllike=loglike(newutmat,as.matrix(data[[resp]]$y))
      oldlprior=logprior(t(oldlambda),Tr1[resp,]%*%olddelta,oldVlambdai)
      newlprior=logprior(t(newlambda),Tr1[resp,]%*%olddelta,oldVlambdai)
      diffvec=newllike+newlprior-(oldllike+oldlprior)
      alpha=min(exp(diffvec), 1)
      #accept or reject
      draw=runif(1)
      accept=0
      if(alpha>draw){accept=1}
      acceptmat[resp,2]=accept
      if(accept==1){
        oldlambdavec[resp]=newlambda
      }
    }
    
    #Store acceptance proportions
    acceptpropb=c(acceptpropb[2:accav],mean(acceptmat[,1]))
    acceptpropl=c(acceptpropl[2:accav],mean(acceptmat[,2]))
    acceptb=mean(acceptpropb)
    acceptl=mean(acceptpropl)
    
    #Draw new values for Delta and Vbeta 
    outbetaup=rmultireg(oldbetamat,Z,Deltabar,ADelta,nub,Vb)
    oldDelta=matrix(outbetaup$B,nr=nZ)
    oldVbeta=outbetaup$Sigma
    oldVbetai=chol2inv(chol(oldVbeta))
    
    #Draw new values for delta and Vlambda
    outlambdaup=rmultireg(oldlambdavec,Tr1,deltabar,Adelta,nul,Vl)
    olddelta=matrix(outlambdaup$B)
    oldVlambda=outlambdaup$Sigma
    oldVlambdai=chol2inv(chol(oldVlambda))
    
    #Draw new values for omega via RW metrop
    newomega=oldomega+rnorm(nTr)*stepo
    oldllike=0
    newllike=0
    for(resp in 1:nresp){
      track=t(Tr2[((resp-1)*nTr+1):(resp*nTr),])
      oldlambda=oldlambdavec[resp]
      oldgamma=track%*%oldomega
      newgamma=track%*%newomega
      oldbeta=matrix(oldbetamat[resp,],nc=1)
      oldutmat=matrix(data[[resp]]$X%*%oldbeta/exp(rep(oldgamma+oldlambda,each=nalt)),nr=nalt)
      newutmat=matrix(data[[resp]]$X%*%oldbeta/exp(rep(newgamma+oldlambda,each=nalt)),nr=nalt)
      oldllike=oldllike+loglike(oldutmat,as.matrix(data[[resp]]$y))
      newllike=newllike+loglike(newutmat,as.matrix(data[[resp]]$y))
    }
    oldlprior=logprior(t(oldomega),t(omegabar),Aomega)
    newlprior=logprior(t(newomega),t(omegabar),Aomega)
    diffvec=newllike+newlprior-(oldllike+oldlprior)
    alpha=min(exp(diffvec), 1)
    #accept or reject
    draw=runif(1)
    accept=0
    if(alpha>draw){accept=1}
    llike=oldllike
    if(accept==1){
      oldomega=newomega
      llike=newllike
    }
    acceptpropo=c(acceptpropo[2:accav],accept)
    accepto=mean(acceptpropo)
    
    
    
    #Store values
    if(r%%keep==0){
      betadraw[r/keep,,]=oldbetamat
      baccept[r/keep]=acceptb
      bstepdraw[r/keep]=stepb
      Deltadraw[r/keep,]=matrix(oldDelta)
      Vbetadraw[r/keep,]=matrix(oldVbeta)
      lambdadraw[r/keep,]=oldlambdavec
      laccept[r/keep]=acceptl
      lstepdraw[r/keep]=stepl
      Vlambdadraw[r/keep,]=matrix(oldVlambda)
      deltadraw[r/keep,]=olddelta
      omegadraw[r/keep,]=oldomega
      ostepdraw[r/keep]=stepo
      oaccept[r/keep,]=accepto
      llikes[r/keep]=llike
    }
    
    #print progress
    #print tte and chart current draw progress
    if(r%%(keep*MCMC$space)==0){
      if(sim==0){
        par(mfrow=c(4,1))
        ctime = proc.time()[3]
        tuntilend = ((ctime - itime)/r) * (R + 1 - r)
        cat(" ", r, " (", round(tuntilend/60, 1), ")", 
            fill = TRUE)
        plot(llikes,type="l",ylab="Log Likelihood")
        matplot(Deltadraw,type="l",ylab="Delta Draws")
        #matplot(rowMeans(lambdadraw),type="l",ylab="lambda Draws")
        matplot(Vlambdadraw,type="l",ylab="lambda Draws")
        matplot(omegadraw,type="l",ylab="omega Draws")
      }
      if(sim==1){
        par(mfrow=c(4,1))
        ctime = proc.time()[3]
        tuntilend = ((ctime - itime)/r) * (R + 1 - r)
        cat(" ", r, " (", round(tuntilend/60, 1), ")",
            fill = TRUE)
        plot(llikes,type="l",ylab="Log Likelihood")
        matplot(Deltadraw,type="l",ylab="Delta Draws"); abline(h=Tvals$tDelta)
        #matplot(rowMeans(lambdadraw),type="l",ylab="lambda Draws")
        #matplot(Vlambdadraw,type="l",ylab="lambda Draws"); abline(h=Tvals$tSigmae)
        matplot(deltadraw,type="l",ylab="delta Draws"); abline(h=Tvals$tdelta)
        matplot(omegadraw,type="l",ylab="omega Draws"); abline(h=Tvals$tomega)
      }
      fsh()
    }
    
    #update stepsizes
    if(r%%accav==0&&r<(.3*R)){
      stepb=stepb*stepupdate(acceptb)
      stepl=stepl*stepupdate(acceptl)
      stepo=stepo*stepupdate(accepto)
    }
    
  } 
  return(list( betadraw=betadraw, baccept=baccept ,bstepdraw=bstepdraw, Deltadraw=Deltadraw, Vbetadraw=Vbetadraw,
               lambdadraw=lambdadraw, laccept=laccept, lstepdraw=lstepdraw, Vlambdadraw=Vlambdadraw,
               deltadraw=deltadraw, omegadraw=omegadraw, oaccept=oaccept, llikedraw=llikes))
}





