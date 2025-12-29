##########################################################################
###### SECR function used for in "Number's count" for lion analysis ######
##########################################################################

##########################################################################
##### SCRi.fn.par1-lion.R is the master function that performs the spatial capture-recapture analysis #######
# scrobj = spatial capture-recapture object that is pulled in for analysis 
# modelno = a simple number so that the output is stored with this suffix
# nc = number of MCMC chains that are to be run for the analysis (about 3-4 should be sufficient)
# ni = total number of MCMC iterations
# burn = number of MCMC iterations to discard
# skip = thin rate (i.e., keep one in skip iterations unless there are memory issues)
# nz = number of "all zero" encounter histories to add
# theta = specify a value between 0.5 and 1. Otherwise, mentioning NA will estimate this parameter 
# Msigma = 0 fits non-spatial model and Msigma = 1 fits spatial model
# statespace (formerly GRID) is the discrete grid of the study area
# Mb = 0 fits model without behavioural response, Mb = 1 fits behavioural model
# Msex = 0 assumes that sex-specific detection probability at activity centre is the same. Msex =1 assumes there are sex-specific differences
# "effort" = Xeff = ntraps x ndays matrix or NULL. Assume log(eff). See "Number's count" paper for interpretation. 
# Xsex = NULL implies that sex ratio is known and can be taken from existing data. Otherwise, this is estimated. 
# ss.prob = NULL implies that pixel-specific covariates are not provided. Leave this as NULL as resource selection is not tested yet on this function. 
# coord.scale = this converts to the necessary unit. Default is in metres. So, 1000 means it is converted to kilometres. 
# area.per.pixel = area per pixel used. This is important for the final density estimation calculation. 
# thinstatespace = this is for thinning the statespace to a desired level. Leave it at 1 unless you have serious memory issues. 
# maxNN = the number of nearest neighbours to be considered for the MCMC jump to locate alternative activity centres for the next step of the MCMC iteration
# dumprate = NOT IMPLEMENTED YET but the plan is to store only specific results. 



SCRi.fn.par1 <-function(scrobj,modelno=1, nc=3,ni=1100,burn=100,skip=2,nz=200,theta=NA,Msigma=1,Mb=0,Msex=0,Msexsigma = 0,Xeff=NULL,Xsex=NULL, ss.prob=NULL,coord.scale=1000,area.per.pixel=1,thinstatespace=1,maxNN=20,dumprate=1000){



call <- match.call()

traps<-scrobj$traps
captures<-scrobj$captures
statespace<-scrobj$statespace
alive<-scrobj$alive
Xd<- scrobj$Xd

captures<- cbind(captures[,4],captures[,2],captures[,3])



if(   length(unique(captures[,2])) != length(min(captures[,2]):max(captures[,2])) ) {
 cat("Error: individuals not numbered sequentially, renumbering them now",fill=TRUE)
 captures[,2]<- as.numeric(factor(captures[,2]))
}

## Alternative e2dist function - utilized from program SPACECAP (Gopalaswamy et al. 2012)
	e2dist <- function(A, B)  {
		xdif <- outer(A[, 1], B[, 1], "-")
		ydif <- outer(A[, 2], B[, 2], "-")
		sqrt(xdif^2 + ydif^2)
	}	
	
 Y<-captures
 traplocs<-traps[,2:3]
 MASK<-as.matrix(traps[,4:ncol(traps)])
 nind<-max(Y[,2])
 T<-dim(MASK)[2]
 M<-nind+nz   # total size of data set
 ntraps<-nrow(traplocs)

totalarea<- nrow(statespace)*area.per.pixel

thinned<- seq(1,nrow(statespace),thinstatespace)
statespace<-statespace[thinned,]
Xd<-Xd[thinned]
goodbad<-statespace[,3]
G<-statespace[,1:2]
G<-G[goodbad==1,]
Gunscaled<-G
Xd<- Xd[goodbad==1]
nG<-nrow(G)

new.area.per.pixel<- totalarea/nG

# The following lines thin the statespace probability weights to match the thinned statespace
#########################################################
	if(!is.null(ss.prob)){
		ss.prob = ss.prob[thinned]
		ss.prob = ss.prob[goodbad==1]
	}
##########################################################


###
###
## Following lines scale the coordinate system based on input coord.scale 
###
###

mgx<-min(traplocs[,1])
mgy<-min(traplocs[,2])
traplocs[,1]<-(traplocs[,1]-min(traplocs[,1]))/coord.scale
traplocs[,2]<-(traplocs[,2]-min(traplocs[,2]))/coord.scale
G[,1]<-(G[,1]-mgx)/coord.scale
G[,2]<-(G[,2]-mgy)/coord.scale



###
# create "Data" vector but with trap mask information
###

msk2<-array(NA,c(nind+nz,T,ntraps))
for(i in 1:(nind+nz)){
msk2[i,1:T,1:ntraps]<-t(MASK[1:ntraps,1:T])
}
msk2<-as.vector(msk2)

###
# expand the data to a 3-d array
Ynew<-array(0,dim=c(nind,T,ntraps))
Ynew[cbind(Y[,2],Y[,3],Y[,1])]<-1
Y<-Ynew

###
### data augmentation
###
Yaug<-array(0,dim=c(nind+nz,T,ntraps))
for(j in 1:nind){
Yaug[j,1:T,1:ntraps]<-Y[j,1:T,1:ntraps]
}


if(!is.null(Xeff)){
Xeffnew<-array(0,dim=c(nind+nz,T,ntraps))
for(j in 1:M){
Xeffnew[j,1:T,1:ntraps]<-t(Xeff)
}
Xeff.tf<-TRUE
}
if(is.null(Xeff)){
Xeffnew<-array(0,dim=c(nind+nz,T,ntraps))
Xeff.tf<-FALSE
}
Xeff<-Xeffnew
if(!is.null(Xsex)){
Xsexnew<-c(Xsex,rep(NA,nz))
}
if(is.null(Xsex)){
Xsexnew<-rep(0,nind+nz)
}
Xsex<-Xsexnew
sex.naflag<-is.na(Xsex)


###
# create covariate of previous capture
###
prevcap<-array(0,c(nind+nz,T,ntraps))
for(i in 1:(nind)){
for(j in 1:ntraps){
tmp<-Yaug[i,1:T,j]
if(any(tmp==1)){
 fst<- min( (1:T)[tmp==1] )
 if(fst<T)
  prevcap[i,(fst+1):T,j]<-1
}
}
}
prevcap<-as.vector(prevcap)
###


### Here "alive" array is made into a full-dimensional array

alive.trues<-array(1, c(nind + nz, T, ntraps))
for(i in 1:nind)
{
for(t in 1:T)
{
alive.trues[i,T,1:ntraps]<-alive[i,t]
}
}
alive.trues<-as.vector(alive.trues)
aliveid<-alive.trues[msk2==1]


##
## vectorize all the data objects
##

arr.trues <- array(TRUE, c(nind+nz,T,ntraps))
idx<-which(arr.trues, arr.ind = TRUE)
y<-as.vector(Yaug)
y<-y[msk2==1 & alive.trues==1]
Xeff<-as.vector(Xeff)
Xeff<-Xeff[msk2==1 & alive.trues==1]
prevcap<-prevcap[msk2==1 & alive.trues==1]

indid.LM<- idx[msk2 == 1, 1]
indid <- idx[msk2 == 1 & alive.trues==1, 1]
repid<-idx[msk2==1 & alive.trues==1,2]
trapid<-idx[msk2==1 & alive.trues==1,3]


getNN<-function(maxNN,G){
###
### Data processing -- this block of code determines a neighborhood
### for every pixel. That information is used in the MCMC updating
### of activity centers. Some tuning may be required
# G = statespace grid

nG<-nrow(G)

NN<-matrix(NA,nrow=nG,ncol=maxNN)

# By having a constant maxNN then the neighborhoods are symmetric

for(i in 1:nG){
od<- sqrt( (G[i,1]-G[,1])^2  +  (G[i,2]-G[,2])^2  )
NN[i,1:maxNN]<-order(od)[1:maxNN]
}
numnn<-rep(0,nrow(NN))
for(i in 1:nG){

# set to NA any element that is not also a neighbor. i.e.,
# remove neighbors for which i is not in their neighborhood.

for(j in 1:ncol(NN)){
if(any(NN[NN[i,j],]==i,na.rm=TRUE)) next
else
 NN[i,j]<-NA
 }
numnn[i]<-sum(!is.na(NN[i,]))
NN[i,1:numnn[i]]<-NN[i,!is.na(NN[i,])]
if(maxNN>numnn[i]){NN[i,(numnn[i]+1):maxNN]<-NA} ## Cleaning up the duplication
}

if(min(numnn)<3)
 cat("State-space grid has isolated or nearly-isolated cells increase maxNN or modify state-space",fill=TRUE)

out<-list(NN=NN,numnn=numnn)
return(out)
}
hld<-getNN(maxNN,G)
NN<-hld$NN
numnn<-hld$numnn

###
###
## this block of code selects starting coordinates for each individual
###
###

centers1<-rep(NA,nind)
for(i in 1:nind){
tt<-t(as.matrix(Yaug[i,,]))
tt<-row(tt)[tt==1]   # which traps was animal captured in
xxx<-as.matrix(traplocs[tt,],ncol=2,byrow=FALSE)  ## coordinates of those traps
av.coord<-colSums(xxx)/nrow(xxx)

dvec<-as.vector(e2dist(matrix(av.coord,ncol=2),G))  # finds closest grid pt
centers1[i]<-(1:length(dvec))[dvec==min(dvec)][1]   # that is initial loc
}
# assigns uncaptured animals a center point;
centers2<-sample(1:nG,M-nind,replace=TRUE, prob=ss.prob)
centers<-c(centers1,centers2)
S<-G[centers,]   # initial locations for all M individuals



## Parallelize the operation (Cross-platform) ##
if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Package 'doParallel' is needed for this function to work. Please install it.")
}
require(doParallel)
require(foreach)

# Register the cores #
# Use makeCluster to be cross-platform (works on Windows and Linux)
cl <- makeCluster(nc)
registerDoParallel(cl)
on.exit(stopCluster(cl), add = TRUE)

# Parallelize the operations using "foreach" paradigm 

results <- foreach(ch=1:nc) %dopar% {

##### Initial values start here for other parameters #####
## Randomized here so that each chain has overdispersed starting values of parameters ## 

if(Msexsigma==0)
bsigma<- rbeta(1,1,40)
if(Msexsigma==1){
bsigmainit <- rbeta(1,1,50)
bsigma <-c(bsigmainit,bsigmainit)	
}

update.theta<-FALSE
if(is.na(theta)){
     theta<- runif(1,0.6,0.9)
     update.theta<-TRUE
 }

lam0<-rgamma(1,0.15,scale=1)
loglam0<-log(lam0)
beta.behave<-0
beta.sex<-0
# start beta1=0 so if there's no covariate this parameter is zeroed out
beta1<-0

psi<-runif(1,0.3,0.8)  # not a good starting values
psi.sex <- mean(Xsex,na.rm=TRUE)
z<-c(rep(1,nind),rbinom(nz,1,psi))
if(sum(sex.naflag)>0)
Xsex[sex.naflag]<-rbinom(sum(sex.naflag),1,psi.sex) 
if(is.null(Xd)){
      Xd<-rep(1,nG)
  }
  beta.den<-0

###### Initial value settings for parameters end here ######

###
###
##   initialization and constructing utility functions 
###
###

lik.fn<-function(lpc,y1){
llvector.new<- -1*exp(lpc)
part2<- exp(exp(lpc[y1])) - 1
part2[part2==0]<-.Machine$double.eps
llvector.new[y1]<- llvector.new[y1]  + log(part2)
llvector.new
}


trapgridbig<-traplocs[trapid,]   # stretches out the trap coord matrix
y1<-y==1
c1<- (S[indid,1]-trapgridbig[,1])^2
c2<- (S[indid,2]-trapgridbig[,2])^2

## print("Error trap 1.7")

gof.new<-gof.data<-rep(NA,(ni-burn)/skip)

out<-matrix(NA,nrow=(ni-burn)/skip,ncol=15)
dimnames(out)<-list(NULL,c("bsigma","sigma","bsigma2","sigma2","lam0","beta.behave","beta1(effort)","beta.sex","psi","psi.sex","Nsuper","theta","beta.density","D","D.adj"))
zout<-matrix(NA,nrow=(ni-burn)/skip,ncol=M)
Sout<-matrix(NA,nrow=(ni-burn)/skip,ncol=M)
LLout<-matrix(NA, nrow=(ni-burn)/skip,ncol=2)
Sexout <- matrix(NA, nrow=(ni-burn)/skip, ncol=M)


m<-1


LM1<-LM2<-matrix(0,nrow=M,ncol=length(indid.LM)/M)
ones<-rep(1,ncol(LM1))

if(Msexsigma==0)
lp.sigma<-Msigma*bsigma*(c1+c2)^theta
if(Msexsigma==1)
lp.sigma<-bsigma[Xsex[indid]+1]*(c1+c2)^theta

acc.count<-0
delta<-.05

start.time <- Sys.time()

	for(i in 1:ni){

########################
########################
########################
########################

# PART 1 OF THE MCMC ALGORITHM UPDATES THE REGRESSION PARAMETERS. FOR THIS MODEL
# THE REGRESSION PARAMETERS ARE (1) INTERCEPT (2) EFFECT OF PREVIOUS CAPTURE
# (BEHAVIORAL RESPONSE) (3) THE SPATIAL PARAMETER "sigma"
### Updating parameters here should only involve individuals with z = 1 (i.e., members of the population)

### update loglam0
repeat{
  lp<-   loglam0 + Mb*beta.behave*prevcap - lp.sigma + beta1*Xeff + Msex*beta.sex*Xsex[indid]
  loglam0c<-rnorm(1,loglam0,.1)
  lpc<-  loglam0c + Mb*beta.behave*prevcap - lp.sigma + beta1*Xeff  + Msex*beta.sex*Xsex[indid]
  llvector<-lik.fn(lp,y1)
  llvector.new<-lik.fn(lpc,y1)
  if(length(llvector.new[!is.finite(llvector.new)])<1) break
}

LM1[aliveid==1]<-llvector.new
LM2[aliveid==1]<- llvector

### Ensure that the likelihood does not become zero ###
mh.ratio <- expm1(sum(( (LM1[z==1,]-LM2[z==1,])%*%ones)))+1
if(mh.ratio==0){mh.ratio <- .Machine$double.xmin}
## print("Error trap 2")
if(runif(1)< mh.ratio){
 loglam0<-loglam0c
 lam0<-exp(loglam0)
 llvector<-llvector.new
 lp<-lpc
 LM2<-LM1
}


## update theta
if(update.theta){
lp<-   loglam0 + Mb*beta.behave*prevcap - lp.sigma + beta1*Xeff + Msex*beta.sex*Xsex[indid]
repeat{
## theta must be less than 0.5 
  thetac<-rnorm(1,theta,.02)
  while(thetac<0.5 | thetac>1) {thetac<-rnorm(1,theta,.02)} 	
      
    if(Msexsigma==0)
      {lp.sigmac<-Msigma*bsigma*(c1+c2)^thetac}
    if(Msexsigma==1)
      {lp.sigmac<-bsigma[Xsex[indid]+1]*(c1 + c2)^thetac}
    
    lpc<-  loglam0 + Mb*beta.behave*prevcap - lp.sigmac + beta1*Xeff  + Msex*beta.sex*Xsex[indid]
    llvector<-lik.fn(lp,y1)
    llvector.new<-lik.fn(lpc,y1)
    if(length(llvector.new[!is.finite(llvector.new)])<1) break
}

LM1[aliveid==1]<-llvector.new
LM2[aliveid==1]<- llvector

mh.ratio <- expm1(sum(( (LM1[z==1,]-LM2[z==1,])%*%ones)))+1
if(mh.ratio==0){mh.ratio <- .Machine$double.xmin}

if(runif(1)< mh.ratio){
 theta<-thetac
 lp.sigma<-lp.sigmac
 llvector<-llvector.new
 lp<-lpc
 LM2<-LM1
}


}


### Update bsigma
## do this differently depending on if sigma depends on sex
if(Msexsigma==0){
  repeat{
bsigmac<-exp(rnorm(1,log(bsigma),delta))
lp.sigmac<- Msigma*bsigmac*(c1+c2)^theta
lpc<-  loglam0 + Mb*beta.behave*prevcap - lp.sigmac + beta1*Xeff  + Msex*beta.sex*Xsex[indid]
llvector.new<- lik.fn(lpc,y1)
if(length(llvector.new[!is.finite(llvector.new)])<1) break
}

LM1[aliveid==1]<- llvector.new

mh.ratio <- expm1(sum(( (LM1[z==1,]-LM2[z==1,])%*%ones)))+1
if(mh.ratio==0){mh.ratio <- .Machine$double.xmin}

if(runif(1)< mh.ratio){
 lp.sigma<-lp.sigmac
 bsigma<-bsigmac
 llvector<-llvector.new
 lp<-lpc
 LM2<-LM1
 acc.count<-  acc.count - 1
 
}
else{
 acc.count<-acc.count+1
}
}

if(Msexsigma==1){
  repeat{
bsigmac<-c(exp(rnorm(1,log(bsigma[1]),2*delta)),bsigma[2])
lp.sigmac<- bsigmac[Xsex[indid]+1]*(c1+c2)^theta
lpc<-  loglam0 + Mb*beta.behave*prevcap - lp.sigmac + beta1*Xeff  + Msex*beta.sex*Xsex[indid]
llvector.new<- lik.fn(lpc,y1)
if(length(llvector.new[!is.finite(llvector.new)])<1) break
}

LM1[aliveid==1]<- llvector.new

mh.ratio <- expm1(sum(( (LM1[z==1,]-LM2[z==1,])%*%ones)))+1
if(mh.ratio==0){mh.ratio <- .Machine$double.xmin}

if(runif(1)< mh.ratio){
 lp.sigma<-lp.sigmac
 bsigma<-bsigmac
 llvector<-llvector.new
 lp<-lpc
 LM2<-LM1

# delta<-delta*1.2  ## This is meant to control the mh.ratio acceptance rate. Right now commented. 

}
else{

#   delta<- delta*.8  ## This is meant to control the mh.ratio acceptance rate. Right now commented.

}
repeat{
bsigmac<-c(bsigma[1],exp(rnorm(1,log(bsigma[2]),2*delta)))
lp.sigmac<- bsigmac[Xsex[indid]+1]*(c1+c2)^theta
lpc<-  loglam0 + Mb*beta.behave*prevcap - lp.sigmac + beta1*Xeff  + Msex*beta.sex*Xsex[indid]
llvector.new<- lik.fn(lpc,y1)
if(length(llvector.new[!is.finite(llvector.new)])<1) break
}

# July 2013
LM1[aliveid==1]<- llvector.new

mh.ratio <- expm1(sum(( (LM1[z==1,]-LM2[z==1,])%*%ones)))+1
if(mh.ratio==0){mh.ratio <- .Machine$double.xmin}


if(runif(1)< mh.ratio){
 lp.sigma<-lp.sigmac
 bsigma[2]<-bsigmac[2]
 llvector<-llvector.new
 lp<-lpc
 LM2<-LM1
}

}
## cat("DELTA: ",delta,fill=TRUE)

if(any(bsigma<0) ){
cat("negative bsigma....",fill=TRUE)
 return(0)
}

if(Msex==1){
  repeat{
beta.sexc<- rnorm(1,beta.sex,.1)
lpc<-  loglam0 + Mb*beta.behave*prevcap - lp.sigma + beta1*Xeff  + Msex*beta.sexc*Xsex[indid]
llvector.new<- lik.fn(lpc,y1)
if(length(llvector.new[!is.finite(llvector.new)])<1) break
}

LM1[aliveid==1]<- llvector.new

mh.ratio <- expm1(sum(( (LM1[z==1,]-LM2[z==1,])%*%ones)))+1
if(mh.ratio==0){mh.ratio <- .Machine$double.xmin}



if(runif(1)< mh.ratio){
 beta.sex<- beta.sexc
 llvector<-llvector.new
 lp<-lpc
 LM2<-LM1   # LM2 is current value
}
}

### This code allows a "local behavioral effect" as opposed to a conventional
###  "global effect" which would not be trap-specific.

if(Mb==1){
  repeat{
beta.behave.c<- rnorm(1,beta.behave,.1)
lpc<-  loglam0 + Mb*beta.behave.c*prevcap - lp.sigma + beta1*Xeff  + Msex*beta.sex*Xsex[indid]
llvector.new<- lik.fn(lpc,y1)
if(length(llvector.new[!is.finite(llvector.new)])<1) break
}

LM1[aliveid==1]<- llvector.new

mh.ratio <- expm1(sum(( (LM1[z==1,]-LM2[z==1,])%*%ones)))+1
if(mh.ratio==0){mh.ratio <- .Machine$double.xmin}



if(runif(1)< mh.ratio){
 beta.behave<- beta.behave.c
 llvector<-llvector.new
 lp<-lpc
 LM2<-LM1   # LM2 is current value
}
}



if(Xeff.tf){
### This block of code below deals with effort.
  repeat{
beta1c<- rnorm(1,beta1,.1)
lpc<-  loglam0 + Mb*beta.behave*prevcap - lp.sigma + beta1c*Xeff   + Msex*beta.sex*Xsex[indid]
llvector.new<- lik.fn(lpc,y1)
if(length(llvector.new[!is.finite(llvector.new)])<1) break
}

LM1[aliveid==1]<- llvector.new

mh.ratio <- expm1(sum(( (LM1[z==1,]-LM2[z==1,])%*%ones)))+1
if(mh.ratio==0){mh.ratio <- .Machine$double.xmin}


if(runif(1)< mh.ratio){
 beta1<- beta1c
 llvector<-llvector.new
 lp<-lpc
 LM2<-LM1   # LM2 is current value
}
}

##############################
##############################
########################
########################
########################
# PART 2 OF THE MCMC ALGORITHM UPDATES THE DATA AUGMENTATION PARAMETERS
# THIS INCLUDES THE LATENT "z" VARIABLES AS WELL AS THE CRITICAL
# PARAMETER "psi"
########################
########################

# This is the data augmentation part. This code updates each z[i]
# for all i=1,2,...M individuals. z[i]=1 is a "real" animal that lives
# in S, whereas z[i]=0 are excess zeros in the data set
# this is an application of Bayes rule to get Pr(z=1| y[i,,]=0)

probz<- exp( rowsum(llvector[indid>nind],indid[indid>nind]) ) # only for nind+1...M
probz<- (probz*psi )/(probz*psi + (1-psi))
z[(nind+1):M]<-rbinom(M-nind,1,probz)
psi<-rbeta(1,1+sum(z),1+M-sum(z))

###
###  This updates SEX identifier  and sex ratio
###

if(!is.null(Xsex)){
# here we only switch sex state if sex is missing for an individual
tmp.sex<-Xsex
tmp.sex[sex.naflag]<- 1-Xsex[sex.naflag]

if(Msexsigma==0)
lp.sigmac<-Msigma*bsigma*(c1+c2)^theta
if(Msexsigma==1)
lp.sigmac<-bsigma[tmp.sex[indid]+1]*(c1+c2)^theta

lpc<-  loglam0 + Mb*beta.behave*prevcap - lp.sigmac + beta1*Xeff  + Msex*beta.sex*tmp.sex[indid]     ## + sex contribution if sex is in the model!
llvector.new<-lik.fn(lpc,y1)
lik.othersex<- expm1(rowsum(llvector.new,indid))+1

lik.sex<- expm1( rowsum(llvector,indid) )+1 # need to do this only for MISSING SEX is done here for all animals
prior.curr<-(psi.sex^Xsex)*((1-psi.sex)^(1-Xsex))
prior.cand<-(psi.sex^tmp.sex)*((1-psi.sex)^(1-tmp.sex))

liksexcand <-lik.othersex*prior.cand
liksexcurr <- lik.sex*prior.curr

liksexcand[liksexcand==0] <- .Machine$double.xmin
liksexcurr[liksexcand==0] <- .Machine$double.xmin

## print("Error trap 5")
swtch<- sex.naflag & (runif(M,0,1)< (liksexcand/liksexcurr))
Xsex[swtch]<- 1-Xsex[swtch]
psi.sex<-rbeta(1,.1+sum(Xsex[z==1]),.1+sum(z)-sum(Xsex[z==1]))

if(Msexsigma==0)
lp.sigma<-Msigma*bsigma*(c1+c2)^theta
if(Msexsigma==1)
lp.sigma<-bsigma[Xsex[indid]+1]*(c1+c2)^theta

lp<-  loglam0 + Mb*beta.behave*prevcap - lp.sigma + beta1*Xeff  + Msex*beta.sex*Xsex[indid]
llvector.new<- lik.fn(lp,y1)
####LM1[1:length(LM1)]<- llvector.new
## July 2013
LM1[aliveid==1]<- llvector.new
llvector<-llvector.new
LM2<-LM1
}

########################
########################
## PART III OF THE ALGORITHM -- THIS BLOCK OF CODE UPDATES THE
## ACTIVITY CENTERS.
########################
########################


##################################################################
### Gets new centers with probability weightings - 10/9/13 JFG
#	if (!is.null(ss.prob)){
#		newcenters = rep(NA, M)
#		for (gc in 1:M){
#			newcenters[gc] <- sample(1:numnn[centers[gc]],size=1,replace=TRUE, prob=ss.prob[na.omit(unique(NN[centers[gc],]))])
#      	}
#	}
#	 if (is.null(ss.prob)){
#		 newcenters <- ceiling(runif(M,0,numnn[centers]))
#       }

## THIS IS CURRENTLY TURNED OFF SINCE WE ARE NOT FOCUSING ON THE RESOURCE SELECTION FUNCTION ## 

###################################################################	

repeat{
  newcenterschoice <- ceiling(runif(M,0,numnn[centers]))
  newcenters<- NN[cbind(centers,newcenterschoice)]

  qnew<- 1/numnn[centers]
  qold<- 1/numnn[newcenters]
  
  Sc<-G[newcenters,]
  c1c<- (Sc[indid,1]-trapgridbig[,1])^2
  c2c<- (Sc[indid,2]-trapgridbig[,2])^2
  
  if(Msexsigma==0)
    lp.sigmac<-Msigma*bsigma*(c1c+c2c)^theta
  if(Msexsigma==1)
    lp.sigmac<-bsigma[Xsex[indid]+1]*(c1c+c2c)^theta
  lpc<- loglam0+ Mb*beta.behave*prevcap - lp.sigmac + beta1*Xeff + Msex*beta.sex*Xsex[indid]
  
  llvector.new<- lik.fn(lpc,y1) 
  if(length(llvector.new[!is.finite(llvector.new)])<1) break
}

LM1[aliveid==1]<- llvector.new
likdiff<- (LM1-LM2)%*%ones

likdiff[z==0]<-0   # all other things equal, this lines sets acceptance prob to 1 for z=0 guys

logprior.new<-    Xd[newcenters]*beta.den
logprior.old<-    Xd[centers]*beta.den

likdiff<-likdiff + log(qold/qnew)  + (logprior.new-logprior.old)

accept.ratio <- expm1(likdiff)+1
accept.ratio[accept.ratio==0] <- .Machine$double.xmin

accept<- runif(M)< accept.ratio

S[accept,]<-Sc[accept,]
centers[accept]<-newcenters[accept]
c1<- (S[indid,1]-trapgridbig[,1])^2
c2<- (S[indid,2]-trapgridbig[,2])^2
LM2[accept,]<-LM1[accept,]


### This block of code updates a density covariate parameter
beta.den.c<- rnorm(1,beta.den,.25)
numerator.new <- exp(Xd*beta.den.c)
loglik.new <- Xd[centers]*beta.den.c - log(sum(numerator.new))
numerator.old <- exp(Xd*beta.den)
loglik.old <-  Xd[centers]*beta.den - log(sum(numerator.old))
if(  runif(1)< (expm1( sum(loglik.new) - sum(loglik.old) )+1)){
 beta.den<-beta.den.c
}


####
####
### update some objects that were left behind earlier
####
####
####

if(Msexsigma==0)
lp.sigma<-Msigma*bsigma*(c1+c2)^theta
if(Msexsigma==1)
lp.sigma<-bsigma[Xsex[indid]+1]*(c1+c2)^theta

lp<-  loglam0 + Mb*beta.behave*prevcap - lp.sigma + beta1*Xeff  + Msex*beta.sex*Xsex[indid]
llvector<- lik.fn(lp,y1)
LM2<-LM1


########################
########################
########################
##
## save output if iteration is past burn-in and using a thin-rate of "skip"
##
########################
########################
########################

if( (i>burn) & (i%%skip == 0) ){

sigma<- sqrt(1/(2*bsigma))
if(Msexsigma == 0){
sigmatmp<-c(sigma,sigma)
bsigmatmp<-c(bsigma,bsigma)
}
else{
sigmatmp<-sigma
bsigmatmp<-bsigma
}

####################################################################
## These are some updates for computing Goodness-Of-Fit statistics##
####################################################################
# mean of "real individuals" (not fixed zeros)

realguys<- z[indid]==1  # TRUE/FALSE vector


logmu<- loglam0 + Mb*beta.behave*prevcap - lp.sigma + beta1*Xeff  + Msex*beta.sex*Xsex[indid]
mu <- (1-exp(-exp(logmu[realguys])))

mu.array <- array(1-mu, dim=c(max(indid),max(repid), max(trapid)))
mu.array[as.matrix(captures[,c(2,3,1)])] <- 1 - mu.array[as.matrix(captures[,c(2,3,1)])]
logmu.array = log(mu.array)
n.nolik = sum(is.infinite(logmu.array))
logmu.array[is.infinite(logmu.array)]=.Machine$double.min.exp

newy <- rbinom(sum(realguys),1,mu)
gof.stats<-cbind(y[realguys],newy,mu)
gof.stats <- sqrt(rowsum(gof.stats, indid[realguys], reorder=FALSE))  
gof.data[m] <- sum((gof.stats[,1] - gof.stats[,3])^2)
gof.new[m]  <- sum((gof.stats[,2] - gof.stats[,3])^2)

####################################################################

density<- sum(z)/totalarea
densityADJ <- sum(z)/(nG*area.per.pixel)

zout[m,]<-z
Sout[m,]<- centers
Sexout[m,] <- Xsex

#### Compute the likelihood of the data given the model/parameter values

LLout[m,]<- c(sum(logmu.array), n.nolik) #prod(c(mu.array))
out[m,]<-c(bsigmatmp[1],sigmatmp[1],bsigmatmp[2],sigmatmp[2],
lam0, beta.behave, beta1,beta.sex,psi,psi.sex,sum(z),theta,beta.den,density,densityADJ)

if(m%%dumprate==0){
#write a file here - NOT IMPLEMENTED YET
}
m<-m+1
}


}

end.time <- Sys.time()
duration <- end.time-start.time

#### vector of model effects
parms.2.report<-
    c(TRUE,TRUE,Msexsigma==1,Msexsigma==1,
TRUE, Mb==1,
Xeff.tf,
 Msex==1, TRUE, (Msex==1|Msexsigma==1), TRUE, update.theta,
sum(Xd)>0,
TRUE ) 

out<- list(mcmchist=out,G=G,Gunscaled=Gunscaled,traplocs=traplocs,Sout=Sout,zout=zout,Sexout=Sexout, likelihood=LLout,statespace=statespace,gof.data=gof.data,gof.new=gof.new,call=call,parms2report=parms.2.report) 

### Extract output into files ###
ts <- format(Sys.time(), "%y%m%d_%H%M%S")
## Create a folder ##
folderName <- paste("NLDresultsModel_",modelno,"_", ts,"CH",ch,sep="")
dir.create(path=folderName)

## Store MCMC history ##
fname <- paste(folderName,"/MLDmcmchist_",ts,"CH",ch,".csv",sep="")
write.csv(file=fname,out$mcmchist) 
## Store activity centres ##
fname <- paste(folderName,"/AcCentres_",ts,"CH",ch,".csv",sep="")
write.csv(file=fname,out$Sout)
## Store real individuals ##
fname <- paste(folderName,"/RealIndividuals_",ts,"CH",ch,".csv",sep="")
write.csv(file=fname,out$zout)
## Store sex of all individuals ##
fname <- paste(folderName,"/SexIndividuals_",ts,"CH",ch,".csv",sep="")
write.csv(file=fname,out$Sexout)
## Store gof data ##
fname <- paste(folderName,"/gofdata_",ts,"CH",ch,".csv",sep="")
write.csv(file=fname,out$gof.data)
## Store gof new ##
fname <- paste(folderName,"/gofnew_",ts,"CH",ch,".csv",sep="")
write.csv(file=fname,out$gof.new)
## Store unscaled statespace for plotting ##
fname <- paste(folderName,"/SSunscaled_",ts,"CH",ch,".csv",sep="")
write.csv(file=fname,out$Gunscaled)
## Store likelihood values to calculated Deviance, BIC, Bayes Factors or other model selection criteria ##
fname <- paste(folderName,"/LogLikelihood_",ts,"CH",ch,".csv",sep="")
write.csv(file=fname,out$likelihood)
## Info file ##
fname <- paste(folderName,"/Info_",ts,"CH",ch,".txt",sep="")
sink(file=fname)
cat("Analysis Duration:",duration)
sink()

}


}
