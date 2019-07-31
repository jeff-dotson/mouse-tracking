#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#   Task-level error scale modeling through tracking
#   Version 3 - Heterogeneous Omega with error in scale
#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


library(gtools)
library(bayesm)
library(gdata)

setwd("/Users/rogerbailey/Dropbox/Models/HeteroOmegaWithError")

#Use simulated data
sim=1

if(sim==1){
  set.seed{22)
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
  
  #Delta coefficients (upper level for betas)
  tDelta=matrix(c(runif(natt*nlvl,-2,3),.5),nc=1)
  tSigmab=diag(nbeta)*.25 +.25 #true cov for betas
  #Omega coefficients (upper level for omegas)
  tOmega=matrix(runif(nTr,-1,2),nc=1)
  tSigmao=diag(nTr)*.25 +.25 #true cov for omegas
  # Generate individual level betas for each respondent.
  tbetas=matrix(double(nresp*nbeta),nc=nbeta)
  for (resp in 1:nresp){
    # Generate individual-level betas.
    tbetabar=tDelta%*%Z[resp,]
    tindbeta=tbetabar + t(chol(tSigmab))%*%rnorm(nbeta)
    tbetas[resp,]=tindbeta
  }
  # Generate individual level omegas for each respondent.
  tomegas=matrix(double(nresp*nTr),nc=nTr)
  for (resp in 1:nresp){
    # Generate individual-level omegas.
    tomegabar=tOmega%*%Z[resp,]
    tindomega=tomegabar + t(chol(tSigmao))%*%rnorm(nTr)
    tomegas[resp,]=tindomega#converts tracking into gamma
  }
  #draw distribution parameter error scale (error in lambda)
  tSigmae=.25
  #delta coefficients
  tdelta=matrix(runif(nTr),nr=nTr)*4-1#converts mean tracking into lambda
  #Determine errors
  terrs=rnorm(nresp)*tSigmae
  #calculate lambda, gamma and total error for each respondent/task
  tscalecomp=matrix(double((ntask+2)*nresp),nc=ntask+2)
  tscalecomp[,1]=Tr1%*%tdelta
  tscalecomp[,2]=terrs
  tscalecomp[,3:(ntask+2)]=matrix(colSums(matrix(matrix(rep(t(tomegas),ntask),nc=ntask)*Tr2,nr=nTr)),nc=ntask)
  tscale=exp(tscalecomp[,3:(ntask+2)]+matrix(rep(tscalecomp[,1],ntask),nc=ntask)+matrix(rep(tscalecomp[,2],ntask),nc=ntask))
  
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
       nbeta, nZ, nTr, tscalecomp, tbetas, tbetabar, tSigmab, tSigmae, 
       tomegas, tOmega, tdelta, tDelta, sim, sure=T)
  #save.image("Simdata.RData")
}else{load("RealData.RData")}


#Set priors
Deltabar=matrix(double(nZ*nbeta),nr=nZ)
ADelta=.01*diag(nZ)
deltabar=matrix(double(nTr),nr=nTr)
Adelta=.01*diag(nTr)
nub=nbeta+3
Vb=diag(nbeta)*nub
nul=4
Vl=matrix(nul)
Omegabar=matrix(double(nZ*nTr),nr=nZ)
AOmega=.01*diag(nZ)
nuo=nTr+3
Vo=diag(nTr)*nuo
R=20e4
accav=50
keep=100
space=5

Prior=list(Deltabar=Deltabar, ADelta=ADelta, Omegabar=Omegabar, 
           AOmega=AOmega, deltabar=deltabar, Adelta=Adelta, nub=nub, Vb=Vb, nul=nul, Vl=Vl, nuo=nuo, Vo=Vo)
MCMC=list(R=R, accav=accav, keep=keep, space=space)
Tvals=list(tscalecomp=tscalecomp, tbetas=tbetas, tDelta=tDelta, tSigmab=tSigmab, tSigmae=tSigmae, 
           tomegas=tomegas, tdelta=tdelta, tOmega=tOmega)
ycomb=NULL
for(i in 1:nresp){ycomb=c(ycomb, data[[i]]$y)}
tester=double(nresp*ntask)
tester[ycomb==4]=1
mean(tester)
#run the program
source("ESMEstHeteroOmegaWithError.R")
output=ESMEstHeteroOmegaWithError(data,Prior,MCMC,sim,Tvals)

plot.bayesm.mat(output$Deltadraw, burnin=0, tvalues=tDelta)
plot.bayesm.mat(output$Omegadraw, burnin=0, tvalues=tOmega)
plot.bayesm.mat(output$deltadraw, burnin=0, tvalues=tdelta)

for(resp in 1:5){
  plot.bayesm.mat(output$betadraw[,resp,], burnin=0, tvalues=tbetas[resp,])
}

  
View(output$Vlambdadraw)

  
}