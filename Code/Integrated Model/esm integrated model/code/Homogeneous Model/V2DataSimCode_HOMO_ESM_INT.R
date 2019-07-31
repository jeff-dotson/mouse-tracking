#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#   Error scaling trhough tracking model homogeneous
#                 Data Generation Code
#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


library(gtools)
library(bayesm)

setwd("/Users/rogerbailey/Dropbox/ESM Integrated Model/Code")

set.seed(77)

nresp=200# number of respondents
nalt=4# number of alternatives per task (inc outside)
natt=3#number of discretely-leveled attributes (inc dummied price)
nlvl=3#number of levels for the attributes(differenced)
ntask=3#number of choice tasks
ndesn=100#number of different designs
nbeta=natt*nlvl+1
nZ=1#number of upper level covariates
nTr=3#number of mouse-tracking variables


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
#Draw true values of parameters:

#draw homogeneous alpha for linear simple function form of engagement 
talpha=(runif(1)+1)/ntask

#construct engagement changes for each individual (homoegeneous)
teng=matrix(c(0:(ntask-1))*talpha,nc=1)

#draw homogeneous beta 
tbeta=matrix(c(runif(natt*nlvl,-2.5,3),-1))#true beta

#coeffcients for Lambda
tOmega=matrix(runif(nTr*2)/ntask,nc=nTr)#converts tracking into gamma

#construct covariance matrix for tracking variables
tSigmaTr=diag(runif(nTr,.5,1))*mean(cbind(double(ntask)+1,teng)%*%tOmega)*.01
#compute tracking variables for each individual

Tr=array(0,dim=c(nresp,ntask,nTr))#true partworth values
for(i in 1:nresp){
  Tr[i,,]=cbind(double(ntask)+1,teng)%*%tOmega+t(t(chol(tSigmaTr))%*%matrix(rnorm(nTr*ntask),nr=nTr))
}

#construct utilities for all alternatives
utmat=matrix(0,nr=ntask*nalt,nc=nresp)
for(i in 1:nresp){
  constantut=desnarray[,,desnvec[i]]%*%tbeta
  utmat[,i]=constantut*rep(matrix((1+teng)^(-1)), each=nalt)
}


#construct choices and save data in list
data=NULL
taskmat=matrix(utmat,nr=nalt)+matrix(-log(-log(runif(
  ntask*nalt*nresp))),nr=nalt)
ycomb=apply(taskmat,2,which.max)
for(i in 1:nresp){
  tempy=ycomb[((i-1)*ntask+1):(i*ntask)]
  data[[i]]=list(X=desnarray[,,desnvec[i]],y=tempy,Tr=Tr[i,,])
}
rm(altmat)
save.image("Simdata_HOMO_ESM_Int.RData")



