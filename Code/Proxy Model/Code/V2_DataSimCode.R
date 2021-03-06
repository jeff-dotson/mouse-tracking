#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#   Error scaling through tracking/emotional recognition
#   Version 4
#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


library(gtools)
library(bayesm)
library(gdata)

setwd("/Users/rogerbailey/Desktop/Projects/ESM")

#Use simulated data
Sim=1

if(Sim==1){
  set.seed(77)
  nresp=200# number of respondents
  nalt=4# number of alternatives per task (inc outside)
  natt=3#number of discretely-leveled attributes (inc dummied price)
  nlvl=4#number of levels for the attributes(differenced)
  ntask=10#number of choice tasks
  ndesn=100#number of different designs
  nbeta=natt*nlvl+1
  nZ=1#number of upper level covariates(including unity)
  nTr=3#number of tracking variables

  data=NULL
  #*************************************************************
  #construct array of design matrices:
  #construct list of all alternative (without price)
  altmat=rbind(diag(nlvl),double(nlvl))
  for(i in 1:(natt-1)){
    ind=rep(1,nrow(altmat))
    altmat=cbind(rbind(diag(nlvl),double(nlvl))%x%ind,
                 rep(1,nlvl+1)%x%altmat)
  }
  #construct array of designs
  altsample=matrix(rbind(matrix(sample(c(1:nrow(altmat)),ntask*
          (nalt-1)*ndesn,replace=T),nr=nalt-1),rep(nrow(altmat),
          ntask*ndesn)),nc=ndesn)
  desnarray=array(0,dim=c(ntask*nalt,ncol(altmat)+1,ndesn))
  for(i in 1:ndesn){
    desnarray[,,i]=cbind(altmat[altsample[,i],],rep(
      c(double(nalt-1),1),ntask))
  }
  #Draw design for each individual
  desnvec=matrix(sample(c(1:ndesn),nresp,replace=T))
  #*************************************************************
  #Draw tracking variables and other covariates
  #Draw other covariates
  Z = matrix(rep(1,nresp),nc=1)
  #Draw tracking covariates (resp and task level)
  tmat=matrix(runif(nresp*ntask*nTr),nc=ntask)
  #construct demeaned tracking variables total and individual
  trackmeanind=matrix(rowMeans(tmat), nc=nresp)
  trackmeantot=rowMeans(trackmeanind)
  Tr1=t(trackmeanind-matrix(rep(trackmeantot,nresp),nc=nresp))
  Tr2=tmat-matrix(rep(trackmeanind,ntask),nc=ntask)
  #*************************************************************
  #Draw true values of parameters:
  #draw distribution parameter error scale (error in lambda)
  tSigmae=.1
  #Omega coeffcients
  tOmega=matrix(runif(nTr),nr=nTr)#converts tracking into gamma
  #delta coefficients
  tdelta=matrix(runif(nTr),nr=nTr)#converts mean tracking into lambda
  #Determine errors
  terrs=rnorm(nresp)*tSigmae
  #calculate lambda, gamma and total error for each respondent/task
  tscalecomp=matrix(double((ntask+2)*nresp),nc=ntask+2)
  tscalecomp[,1]=Tr1%*%tdelta
  tscalecomp[,2]=terrs
  tscalecomp[,3:(ntask+2)]=matrix(colSums(matrix(Tr2*matrix(rep(tOmega,ntask*nresp),nc=ntask),nr=3)),nc=ntask)
  tscale=exp(tscalecomp[,3:(ntask+2)]+matrix(rep(tscalecomp[,1],ntask),nc=ntask)+matrix(rep(tscalecomp[,2],ntask),nc=ntask))
  #Delta coefficients (upper level for betas)
  tDelta=matrix(c(runif(natt*nlvl,-2,3),2),nc=1)#true betabar with out good
  tSigmab=rwishart(natt*nlvl+3,(natt*nlvl+3)*diag(natt*nlvl+1))$IW #true cov
  # Generate individual level betas for each respondent.
  tbetas=matrix(double(nresp*nbeta),nc=nbeta)
  for (resp in 1:nresp){
    # Generate individual-level betas.
    tbetabar=tDelta%*%Z[resp,]
    tindbeta=tbetabar + t(chol(tSigmab))%*%rnorm(nbeta)
    tbetas[resp,]=tindbeta
  }
  #*************************************************************
  #Generate utilities and choices
  utmat=matrix(0,nr=ntask*nalt,nc=nresp)#matrix of alternative utilities
  ycomb=NULL
  Xcomb=NULL
  data=NULL
  for (resp in 1:nresp){
    desn=desnarray[,,desnvec[resp]]
    Xcomb=rbind(Xcomb,desn)
    #Compute utilities scaled by task
    ut=(desn%*%tindbeta)/rep(tscale[resp,],each=nalt) + matrix(-log(-log(runif(ntask*nalt))),nc=1)
    utmat[,i]=ut
    #Construct choices
    y=apply(matrix(ut,nr=nalt),2,which.max)
    ycomb=c(ycomb,y)
    #Store data
    data[[resp]]=list(X=desn,y=y,Z=Z[resp,],Tr1=Tr1[resp,],Tr2=Tr2[((resp-1)*nTr+1):(resp*nTr),])
  }
  keep(data, nresp, nalt, natt, nlvl, ntask, ndesn, 
       nbeta, nZ, nTr, tscalecomp, tbetas, tbetabar, tSigmab, tSigmae, tOmega, tdelta, tDelta, sure=T)
  save.image("Simdata.RData")
}else{load("RealData.RData")}







