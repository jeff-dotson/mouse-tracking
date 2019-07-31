#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#   Error Scaling with Mousetracking - Data Generation Code
#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


library(gtools)
library(bayesm)

setwd("/Users/rogerbailey/Desktop/Projects/ESM/Code")

set.seed(77)

nresp=200# number of respondents
nalt=4# number of alternatives per task (inc outside)
natt=5#number of discretely-leveled attributes (inc dummied price)
nlvl=3#number of levels for the attributes(differenced)
ntask=100#number of choice tasks
ndesn=100#number of different designs
nbeta=natt*nlvl+1
nZ=1#number of upper level covariates
nM=3#number of mouse-tracking variables


data=NULL
#*************************************************************
#draw covariates
Z=matrix(rep(1,nresp),nr=1)
if(nZ>1){
for(i in 1:(nZ-1)){
  Z=rbind(Z,runif(nresp)-.5) 
}
}

#*************************************************************
#draw mouse tracking variables
M=matrix(runif(nresp*nM),nc=nresp)


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
#Draw true values of parameters:

#upper level coefficients for betas
temp1=matrix(c(runif(natt*nlvl,-1.5,2),-1))#true betabar 
tDelta=matrix(temp1,nc=1)#Delta matrix that converts covariates to predicted betas
if(nZ>1){
  for(i in 1:(nZ-1)){
  tDelta=cbind(tDelta,runif(natt*nlvl+1))
  }
}
tSigma=diag(natt*nlvl+1)+.5#rwishart(natt*nlvl+3,(natt*nlvl+3)*diag(natt*nlvl+1))$IW #true cov 


#coeffcients for gamma
tOmega=matrix(runif(nM)*2/ntask,nr=1)#converts tracking into gamma

#draw true betas for each individual
tbetas=matrix(0,nr=nresp,nc=natt*nlvl+1)#true partworth values 
for(i in 1:nresp){
  tbetabar=tDelta%*%Z[,i]
  tbetas[i,]=tbetabar+t(chol(tSigma))%*%rnorm(natt*nlvl+1)
}

#compute true gamma for each individual
tgammas=matrix(0,nr=nresp,nc=1)#true partworth values
for(i in 1:nresp){
  tgammas[i]=tOmega%*%M[,i]
}

#construct utilities for all alternatives
utmat=matrix(0,nr=ntask*nalt,nc=nresp)
for(i in 1:nresp){
  constantut=desnarray[,,desnvec[i]]%*%matrix(tbetas[i,],nc=1)
  utmat[,i]=constantut*matrix((1+tgammas[i]*(c(0:(ntask-1))%x%rep(1,nalt)))^(-1))
}


#construct choices and save data in list
data=NULL
taskmat=matrix(utmat,nr=nalt)+matrix(-log(-log(runif(
  ntask*nalt*nresp))),nr=nalt)
ycomb=apply(taskmat,2,which.max)
for(i in 1:nresp){
  tempy=ycomb[((i-1)*ntask+1):(i*ntask)]
  data[[i]]=list(X=desnarray[,,desnvec[i]],y=tempy, Z=Z[,i],M=M[,i])
}
rm(altmat)
save.image("Simdata_BasicErrScaling.RData")


