#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#   Error scaling through tracking/emotional recognition
#   Version 4
#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


library(gtools)
library(bayesm)
library(gdata)
set.seed(77)

ESMemoEst=function(data,Prior,MCMC){
  # FUNCTION OUTPUT:
  # betadraw - posterior draws of individual-level betas
  # baccept - acceptance rate for the beta step
  # bstepdraw - adjusted beta stepsize draws, fixed after burn-in
  # Vbetadraw - posterior draws of the upper level covariance matrix
  # Deltadraw - posterior draws of the upper level coefficient matrix
  # deltadraw - posterior draws of the upper level coefficients for resp-level error scale deviations 
  # omegadraw - posterior draws of the upper level coefficients for task-level error scale deviations 
  # lambdadraw - posterior draws of individual-task-level error scale residuals
  # Vlambdadraw - posterior draws of the upper level variance for the error scale 
  # llikedraw - draws of the log likelihood

  #get variables
  nresp=length(data)# number of respondents
  ntask=length(data[[1]]$y)# number of alternatives per task (inc outside)
  nalt=nrow(data[[1]]$X)/ntask#number of choice tasks
  nbeta=ncol(data[[1]]$X) #number of partworths to be estimated
  nZ=ncol(data[[1]]$Z)#number of upper level covariates(including unity)
  nTr=length(data[[1]]$Tr1)#number of tracking variables
  
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
  
  #Set initial values
  oldDelta=matrix(double(nZ*nbeta),nc=nbeta)#initial betabar
  oldbetamat=matrix(double(nbeta*nresp),nr=nresp)#initial betas
  stepb=.2 #beta stepsize
  oldVbeta=diag(nbeta)#initial Vbeta
  oldVbetai=backsolve(chol(oldVbeta),diag(nbeta))#initial inv chol Vbeta
  acceptpropb=double(accav)+.23 #proportion of beta draws accepted
  oldlambda=matrix(double(nresp),nr=nresp)#initial lambdas
  stepl=.2 #lambda stepsize
  oldVlambda=diag(1)#initial Vbeta
  oldVlambdai=backsolve(chol(oldVlambda),diag(1))#initial inv chol Vbeta
  acceptpropl=double(accav)+.23 #proportion of lambda draws accepted
  olddelta=matrix(double(nTr),nc=1) #initial delta
  oldomega=matrix(double(nTr),nc=1) #initial omega
  stepo=.2 #omega stepsize
  acceptpropo=double(accav)+.23 #proportion of lambda draws accepted
  llike=matrix(double(nresp),nc=1) #log likelihood
  
  #Setup storage
  betadraw=array(double(nbeta*nresp*R/keep),dim=c(R/keep,nresp,nbeta))
  baccept=matrix(double(R/keep),nr=R/keep)
  bstepdraw=matrix(double(R/keep),nr=R/keep)
  lamdadraw=matrix(double(R/keep),nr=R/keep)
  laccept=matrix(double(R/keep),nr=R/keep)
  deltadraw=matrix(double(R/keep*nTr),nr=R/keep)
  omegadraw==matrix(double(R/keep),nr=R/keep)
  oaccept=matrix(double(R/keep),nr=R/keep)
  llikes=matrix(double(nresp*R/keep),nr=R/keep)
  acceptpropbs=matrix(double(R/keep),nr=R/keep)
  
  #start timing
  itime = proc.time()[3]
  for(r in 1:R){
    accept=matrix(double(nresp),nc=nresp)#iter accept storage
    
    #Draw new betas for each respondent
    for(resp in 1:nresp){
      track=t(Tr2[((resp-1)*nTr+1):(resp*nTr),])
      oldlambda=oldlambdamat[resp]
      oldgamma=track%*%oldomega
      #get new and old paramter vectors for RW
      oldbeta=matrix(oldbetamat[resp,],nc=1)
      newbeta=oldbeta+t(chol(oldVbeta))%*%rnorm(nbeta)*stepb
      oldutmat=matrix(data[[resp]]$X%*%oldbeta/exp(rep(oldgamma+oldlambda,each=nalt)),nr=nalt)
      newutmat=matrix(data[[resp]]$X%*%newbeta/exp(rep(oldgamma+oldlambda,each=nalt)),nr=nalt)
      #get likelihood and prior values
      oldllike=loglike(oldutmat,as.matrix(data[[resp]]$y))
      newllike=loglike(newutmat,as.matrix(data[[resp]]$y))
      oldlprior=logprior(t(oldbeta),t(oldbetabar),oldVbetai)
      newlprior=logprior(t(newbeta),t(oldbetabar),oldVbetai)
      diffvec=newllike+newlprior-(oldllike+oldlprior)
      alpha=min(exp(diffvec), 1)
      #accept or reject
      draw=runif(1)
      accept=0
      if(alpha>draw){accept=1}
      accept[resp,1]=accept
      if(accept==1){
        oldbetamat[resp,]=newbeta
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
      accept[resp,2]=accept
      if(accept==1){
        oldlambda[resp,]=newlambda
      }
    }
    
    #Store acceptance proportions
    acceptb=mean(c(acceptpropb[2:accav],mean(accept[,1])))
    acceptl=mean(c(acceptpropl[2:accav],mean(accept[,2])))
    
    #Draw new values for Delta and Vbeta 
    outbetaup=rmultireg(oldbetamat,Z,Deltabar,ADelta,nub,Vb)
    oldDelta=matrix(outbetaup$B)
    oldVbeta=outbetaup$Sigma
    oldVbetai=chol2inv(chol(oldVbeta))
    
    #Draw new values for delta and Vlambda
    outlambdaup=rmultireg(oldlambdamat,Tr1,deltabar,Adelta,nud,Vd)
    delta=matrix(outlambdaup$B)
    oldVlambda=outlambdaup$Sigma
    oldVlambdai=chol2inv(chol(oldVlambda))
    
    #Draw new values for omega via RW metrop
    newomega=oldomega+rnorm(nTr)*stepg
    oldllike=0
    newllike=0
    for(resp in 1:nresp){
      track=t(Tr2[((resp-1)*nTr+1):(resp*nTr),])
      oldlambda=oldlambdamat[resp]
      oldgamma=track%*%oldomega
      newgamma=track%*%newomega
      oldbeta=matrix(oldbetamat[resp,],nc=1)
      oldutmat=matrix(data[[resp]]$X%*%oldbeta/exp(rep(oldgamma+oldlambda,each=nalt)),nr=nalt)
      newutmat=matrix(data[[resp]]$X%*%oldbeta/exp(rep(newgamma+oldlambda,each=nalt)),nr=nalt)
      oldllike=oldllike+loglike(oldutmat,as.matrix(data[[resp]]$y))
      newllike=newllike+loglike(newutmat,as.matrix(data[[resp]]$y))
    }
    oldlprior=logprior(t(oldomega),omegabar,Aomega)
    newlprior=logprior(t(newomega),omegabar,Aomega)
    diffvec=newllike+newlprior-(oldllike+oldlprior)
    alpha=min(exp(diffvec), 1)
    #accept or reject
    draw=runif(1)
    accept=0
    if(alpha>draw){accept=1}
    llike=oldllike
    if(accept==1){
      oldlambda[resp,]=newlambda
      llike=newllike
    }
    acceptpropo=c(acceptpropo[2:accav],accept)
    accepto=mean(acceptpropo)
    
    
    
    #Store values
    if(r%%keep==0){
      betadraws[r/keep,,]=oldbetamat
      betabardraws[r/keep,]=oldbetabar
      acceptpropbs[r/keep]=mean(acceptpropb)
      llikes[r/keep,]=llike
      stepsizebdraws[r/keep]=stepsizeb
    }
    
    #print progress
    #print tte and chart current draw progress
    if(r%%(keep*space)==0){
      par(mfrow=c(2,1))
      ctime = proc.time()[3]
      tuntilend = ((ctime - itime)/r) * (R + 1 - r)
      cat(" ", r, " (", round(tuntilend/60, 1), ")", 
          fill = TRUE)
      plot(rowSums(llikes),type="l",ylab="Log Likelihood")
      matplot(betabardraws,type="l",ylab="Betabar Draws")
      fsh()
    }
    
    #update stepsizes
    if(r%%accup==0&&r<(.3*R)){
      stepsizeb=stepsizeb*stepupdate(mean(acceptpropb))
    }
    
  }  
}







