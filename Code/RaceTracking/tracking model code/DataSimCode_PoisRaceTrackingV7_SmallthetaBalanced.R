#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#   Poisson race tracking and response time- Estimation Code
#   - Data Gen Code
#   May 2016
#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

library(gtools)
library(bayesm)
library(gdata)

setwd("/Users/rogerbailey/Dropbox/Projects/RaceTracking/Tracking Model Code")

set.seed(22)

nresp=200# number of respondents
nalt=4# number of alternatives per task (no outside good)
natt=3#number of discretely-leveled attributes
nlvl=3#number of levels for the attributes(differenced)
ntask=30#number of choice tasks
ndesn=100#number of different designs
nbeta=natt*nlvl
nZ=1#number of upper level covariates
nTr=3#number of tracking variables

data=NULL

#*************************************************************
#construct list of all possible alternatives
altmat=rbind(diag(nlvl),double(nlvl))
for(i in 1:(natt-1)){
  ind=rep(1,nrow(altmat))
  altmat=cbind(rbind(diag(nlvl),double(nlvl))%x%ind,
               rep(1,nlvl+1)%x%altmat)
}

#*************************************************************
#Construct upper level for betas, deltas, lambdas, and thetas
#draw tracking covariates
Omega=cbind(matrix(rtrun(mu=rep(.5,nresp*ntask),sigma=rep(.5,nresp*ntask),a=double(nresp*ntask),b=rep(100,ntask*nresp)),nc=1), runif(nresp*ntask), sample(c(0,1,2),replace=T,nresp*ntask,prob=double(3)+1/3))
ttheta=-2


#true psi and thetas
while(min(ttheta)<.1){
tgamma=matrix(rnorm(nTr)+.25,nc=1)
tpsibar=runif(1)+.5
tSigmapsi=.5
tpsi=rnorm(nresp)*tSigmapsi+tpsibar
ttheta=Omega%*%tgamma + rep(tpsi,each=ntask)
}
ttheta=matrix(ttheta,byrow=T,nc=ntask)

#true upper level cov and mean for betas
tSigmabeta=diag(nbeta)*.1
tbetabar=runif(nbeta)*7-3
#true upper level cov and mean for deltas and lambda (coefficents determining 
#  effect of accessibility and attractiveness of alternatives) draw together
tSigmadeltalambda=diag(4)*.05
tdeltalambdabar=c(runif(1)-1,(runif(1)+.5)*10/ntask,runif(1)-2,(runif(1)-1.5)*5/ntask) 


#*************************************************************
#draw true values of paramters for each respondent
#draw true betas for each respondent
tbetas=t(tbetabar+t(chol(tSigmabeta))%*%matrix(rnorm(nbeta*nresp),nc=nresp))
#draw true K paramter for each respondent
tK=matrix(double(ntask*nresp),nc=ntask)
for(t in 1:ntask){
tK[,t]=matrix(rpois(nresp,ttheta[,t])+1)
}

#*************************************************************
#construct "balanced" choice tasks according to betabar

#create list of all possible tasks and evaluate logit
tasks=matrix(combn(c(1:length(altmat[,1])),nalt),nc=1)
npostasks=length(tasks)/nalt
acceptable=matrix(double(npostasks))
alts=altmat[tasks,]
uts=exp(matrix(alts%*%tbetabar, nr=nalt))
sums=matrix(rep(colSums(uts),nalt),byrow=T,nr=nalt)
uts= uts/sums 
for(t in 1:npostasks){
  order=sort(uts[,t])
  #determine acceptability 
  if((order[nalt]-order[nalt-1]<.2)&(order[nalt]-order[nalt-2]<.4)){
    acceptable[t]=1
  }
}
#list of acceptable tasks (via tbetabar) from all possible tasks
acctasks=matrix(tasks,nr=nalt)[,acceptable==1]

#construct array of designs that are balanced with respect to utility
desnarray=array(0,dim=c(ntask*nalt,ncol(altmat),ndesn))
altsample=matrix(double(ntask*ndesn),nc=ndesn)
for(d in 1:ndesn){
  tsample=sample(c(1:sum(acceptable)),ntask)
  altlist=NULL
  for(t in 1:ntask){
    taskalts=acctasks[,tsample[t]]
    altlist=c(altlist,sample(taskalts,nalt)) #permute tasks
  }
  desnarray[,,d]=altmat[altlist,]
}
#Draw design for each individual
desnvec=matrix(sample(c(1:ndesn),nresp,replace=T))

#*************************************************************
#draw total time for each alternative for each respodnent
#by constructing scaling parameters for each alt and
#task and using in a draw from a gamma distribution

#construct "base rates" for each alternative (anaogous to 
#alternative utilities in MNL model)
tbaserates=matrix(double(nresp*nalt*ntask),nr=nresp)
for(i in 1:nresp){
  tbaserates[i,]=exp(desnarray[,,desnvec[i]]%*%matrix(tbetas[i,],nc=1))
}

#construct adjusted rate for each respondent/task/alternative
trates=matrix(double(nresp*nalt*ntask),nr=nresp)
for(i in 1:nresp){
  trates[i,]=tbaserates[i,] 
}
#draw "response times" for each alternative and task
ttime=matrix(double(nresp*nalt*ntask),nr=nresp)
for(i in 1:nresp){
  shapemat=matrix(rep(tK[i,],each=nalt),nc=1)
  for(j in 1:(ntask*nalt)){
    ttime[i,j]=rgamma(1,shape=shapemat[j],scale=trates[i,j]) 
  }
}

#*************************************************************
#determine choices and save data

#determine choice for each task
taskmat=matrix(t(ttime),nr=nalt)
ycomb=apply(taskmat,2,which.min)
tcomb=apply(taskmat,2,min)
for(i in 1:nresp){
  tempy=ycomb[((i-1)*ntask+1):(i*ntask)]
  tempt=tcomb[((i-1)*ntask+1):(i*ntask)]
  data[[i]]=list(X=desnarray[,,desnvec[i]],y=tempy, time=tempt,Om=Omega[((i-1)*ntask+1):(i*ntask),])
}
keep(data, nresp, nalt, natt, nlvl, ntask, ndesn, 
     nbeta, nZ, nTr, tbetas, tbetabar, tSigmabeta, 
     trates, tK, tgamma, tpsi, tpsibar, sure=T)
save.image("Simdata_PoisRaceTrackingV7_Balanced.RData")
rm(list = ls())

















