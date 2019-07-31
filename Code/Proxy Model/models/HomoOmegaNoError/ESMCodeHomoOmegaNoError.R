#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#   Task-level error scale modeling through tracking
#   Version 4 - Homogeneous Omega with no error in scale
#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


library(gtools)
library(bayesm)
library(gdata)

setwd("/Users/rogerbailey/Dropbox/Models/HomoOmegaNoError")

#Use simulated data
sim=1

if(sim==1){
  set.seed(77)
  nresp=500# number of respondents
  nalt=4# number of alternatives per task (inc outside)
  natt=3#number of discretely-leveled attributes (inc dummied price)
  nlvl=3#number of levels for the attributes(differenced)
  ntask=50#number of choice tasks
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
  #omega coeffcients
  tomega=matrix(runif(nTr),nr=nTr)#converts tracking into gamma
  #delta coefficients
  tdelta=matrix(runif(nTr),nr=nTr)*5#converts mean tracking into lambda
  #calculate lambda, gamma and total error for each respondent/task
  tscalecomp=matrix(double((ntask+1)*nresp),nc=ntask+1)
  tscalecomp[,1]=Tr1%*%tdelta
  tscalecomp[,2:(ntask+1)]=matrix(colSums(matrix(Tr2*matrix(rep(tomega,ntask*nresp),nc=ntask),nr=3)),nc=ntask)
  tscale=exp(tscalecomp[,2:(ntask+1)]+matrix(rep(tscalecomp[,1],ntask),nc=ntask))
  #Delta coefficients (upper level for betas)
  tDelta=matrix(c(runif(natt*nlvl,-2,3),.5),nc=1)#true betabar with out good
  tSigmab=diag(nbeta)*.25 +.25 #true cov
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
    ut=((desn%*%tbetas[resp,])/rep(tscale[resp,],each=nalt)+ matrix(-log(-log(runif(ntask*nalt))),nc=1))
    #Construct choices
    y=apply(matrix(ut,nr=nalt),2,which.max)
    ycomb=c(ycomb,y)
    #Store data
    data[[resp]]=list(X=desn,y=y,Z=Z,Tr1=Tr1,Tr2=Tr2)
  }
  keep(data, nresp, nalt, natt, nlvl, ntask, ndesn, 
       nbeta, nZ, nTr, tscalecomp, tbetas, tbetabar, tSigmab, 
       tomega, tdelta, tDelta, sim, sure=T)
  #save.image("Simdata.RData")
}else{load("RealData.RData")}


#Set priors
Deltabar=matrix(double(nZ*nbeta),nr=nZ)
ADelta=.01*diag(nZ)
omegabar=matrix(double(nTr),nr=nTr)
Aomega=.01*diag(nTr)
deltabar=matrix(double(nTr),nr=nTr)
Adelta=.01*diag(nTr)
nub=nbeta+3
Vb=diag(nbeta)*nub
R=25e3
accav=50
keep=25
space=5

Prior=list(Deltabar=Deltabar, ADelta=ADelta, omegabar=omegabar, 
           Aomega=Aomega, deltabar=deltabar, Adelta=Adelta, nub=nub, Vb=Vb)
MCMC=list(R=R, accav=accav, keep=keep, space=space)
Tvals=list(tscalecomp=tscalecomp, tbetas=tbetas, tDelta=tDelta, tSigmab=tSigmab, 
           tomega=tomega, tdelta=tdelta)
ycomb=NULL
for(i in 1:nresp){ycomb=c(ycomb, data[[i]]$y)}
tester=double(nresp*ntask)
tester[ycomb==4]=1
mean(tester)
#run the program
source("ESMEstHomoOmegaNoError.R")
output=ESMEstHomoOmegaNoError(data,Prior,MCMC,sim,Tvals)

plot.bayesm.mat(output$Deltadraw, burnin=0, tvalues=tDelta)
plot.bayesm.mat(output$omegadraw, burnin=0, tvalues=tomega)
plot.bayesm.mat(output$deltadraw, burnin=0, tvalues=tdelta)

for(resp in 1:5){
  plot.bayesm.mat(output$betadraw[,resp,], burnin=0, tvalues=tbetas[resp,])
}
  
  
  
}