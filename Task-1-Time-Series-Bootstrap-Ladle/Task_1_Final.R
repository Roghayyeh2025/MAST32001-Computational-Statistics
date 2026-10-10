###EXERCIZE 1###
library(readr)
timeseries_dat <- read_table("timeseries.dat.txt", 
                             col_names = FALSE)
View(timeseries_dat)
###PART A###
#centering
data=timeseries_dat
means=apply(data, 2, mean)
datac=sweep(data,2,means,'-') #centered data (I used the convention to center the data by the subtraction of every variable mean)
datac
apply(datac, 2, mean)
#covariance matrix
S=cov(datac) #for the normalization I chose n-1 (instead of n) because in estimating mu with \bar(x) I'm losing one 
#degree of freedom
S
#covariance matrix for the whitening
?svd
sing=svd(S)#with the svd decomposition S=UDV^T
#U and V are the same because S is a covariance matrix and so is symmetric.
D2=diag(sing$d^(-0.5)) #D^(-0.5)
D2
S2=sing$u%*%D2%*%t(sing$u)
S2 #S^(-0.5)
#whitening
datac2=as.matrix(datac)
dataw=datac2%*%S2 #my scaling convention is the multivariate whitening based on the sample covariance matrix.
cov(dataw)#We obtain almost I_10, so the whitening was done successfully
#Now we have whitened variables: they aren't correlated and they have variance=1 and mean=0. We still cannot
#say that they are indipendent between each other, and of course we didn't eliminate the temporal autocorrelation.
#why do we use whitening?
autoc=function(r){
  A=dataw[1:(1000-r),]
  B=dataw[(1+r):1000,] #lagging
  aut=(t(A)%*%B)/(1000-r) #it's the formula sum(t=1 to n-r)(x_t^w)(x_(t+r)^w)^T/(n-r) (we use n-r because we just have n-r available temporal couples)
  aut
}
autoc(1) #it gives us the relation between X_i,t and X_j,t+1. It doesn't have to be symmetric.
#The whitening didn't remove this correlation.
#I still didn't choose the relevant autocovariance matrices.

#
###PART B###
#
#Fixed block bootstrap
block=matrix(rep(0,10000), nrow =1000 )
l=20 #the block length doesn't have to be too small, because we would loose a lot of dependence, but l cannot be
#too big either, since the resampling variability would be too low.
datamat=as.matrix(data)
set.seed(123)
j=1
for (i in 1:50) {
  star=sample(1:981, size=1) #I didn't choose a moving block bootstrap
  block[j:(j+19),]=datamat[star:(star+19),]
  j=j+20
}
block
#stationary bootstrap
stationary=matrix(rep(0,10000), nrow =1000 )
datamat2=matrix(rep(0,20000), nrow=2000) #so that if we sample 1000 we can later add 1 to the block
datamat2[1:1000,]=datamat
datamat2[1001:2000,]=datamat
p=0.05 #to have an expected length of the blocks=20
set.seed(123)
star=sample(1:1000, size=1)
stationary[1,]=datamat2[star,]
set.seed(123)
for (i in 2:1000) {
  if (runif(1)>=p) {
    star=star+1
    stationary[i,]=datamat2[star,]
  } else {
    star=sample(1:1000, size=1)
    stationary[i,]=datamat2[star,]
  }
}
stationary
# The tuning parameters controlling the block length are l and p.
#In a stationary bootstrap the expected length of a block is 1/p since the block length has a geometric 
#distribution with parameter p. If the block doesn't continue we sample a new starting point, while for fixed
#blocks we have that the length is fixed and just have to sample the starting point of the block.
#iid resampling of individual observations is generally inappropriate for a dependent time series because it 
#doesn't preserves the indipendence of the observations.

#functions
#block function
funblock=function(dat,l){
  n=nrow(dat)
  nvar=ncol(dat)
  m=ceiling(n/l) #for example if we have n=105 and l=20 it gives m=6
  block=matrix(rep(0,m*l*nvar), nrow =m*l )#so that we can sample whole blocks
  j=1
  for (i in 1:m) {
    star=sample(1:(n-l+1), size=1) 
    block[j:(j+l-1),]=dat[star:(star+(l-1)),]
    j=j+l
  }
  block[1:n, ]#we just want the first n observations (that can be less than m*l)
}
funblock(datamat,20)
#stationary function
funstat=function(dat,p){
  n=nrow(dat)
  nvar=ncol(dat)
  stationary=matrix(rep(0,n*nvar), nrow =n )
  dat2=matrix(rep(0,n*nvar*2), nrow=n*2)
  dat2[1:n,]=dat
  dat2[(n+1):(2*n),]=dat
  star=sample(1:n, size=1)
  stationary[1,]=dat2[star,]
  for (i in 2:n) {
    if (runif(1)>=p) {
      star=star+1
      stationary[i,]=dat2[star,]
    } else {
      star=sample(1:n, size=1)
      stationary[i,]=dat2[star,]
    }
  }
  stationary
}
funstat(datamat, 0.05)


##graphical summaries
par(mfrow=c(1,2))#we work on two examples to reduce the computational work
acf(data[,4])#strong autocorrelation
acf(data[,3])#weak autocorrelation
datagrap=datamat[,3:4]
acf(datagrap)
#block bootstrap with l=20
block_20=funblock(datagrap,20)
par(mfrow=c(1,1))
acf(block_20)
#block bootstrap with l=5
block_5=funblock(datagrap,5)
acf(block_5)#the temporal dependence of X_4 is less preserved compared to l=20
#block bootstrap with l=50
block_50=funblock(datagrap,50)
acf(block_50) #the temporal dependence of X_4 is really well preserved compared to l=20
#stationary bootstrap with p=0.05
stat_005=funstat(datagrap,0.05)
acf(stat_005)#still good results
#stationary bootstrap with p=0.5
stat_05=funstat(datagrap,0.5)
acf(stat_05)#the expected length decreases and the results are bad: there isn't a good representation of the 
#autocorrelation of X_4.

#function for iid bootstrap
funiid=function(dat){
  n=nrow(dat)
  nvar=ncol(dat)
  iid=matrix(rep(0,n*nvar), nrow = n)
  for (i in 1:n) {
    star=sample(1:n,size = 1)
    iid[i,]=dat[star,]
  }
  iid
}
#graphic
iid=funiid(datagrap)
acf(iid)#our local dependence isn't preserved.
















