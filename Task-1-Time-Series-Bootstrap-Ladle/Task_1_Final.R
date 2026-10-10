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
#The stationary bootstrap gives a stationary series while the block bootstrap doesn't.

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

#I could also check some ACF for other bootstrap samples (repeating the same bootstrap twice for example)


#
###PART C###
#
install.packages('JADE')
library(JADE)
?rjd

#creation of the array 10x10x12 to use the rjd function. We are interested in different (12) autocorrelations 
#because SOBI uses them (while AMUSE uses just lag 1). The choice of 12 follows from the methodology in the paper.
#We work with the whitened data, since it is required in the preprocessing of SOBI.
arr=array(0,c(10,10,12))
for (i in 1:12) {
  arr[,,i]=autoc(i)
}
arr
dim(arr)
#joint diagonalization
arrdiag=rjd(arr)
str(arrdiag)
dim(arrdiag$V)#the orthogonal matrix that is used to diagonalize
arrdiag$V%*%t(arrdiag$V)#It's orthogonal: V*V^T=I. In the paper is called t(U_Tao)
dim(arrdiag$D)#the transformed matrices by the joint diagonalization
#rjd found the orthogonal matrix V so that every D_r=V^T*autoc_tao*V is almost diagonal (for every lag tao from 1 to 12)
#The different D_tao (tao:1,...,12) are not really diagonal, because normally doesn't exist a matrix that can
#diagonalize all 12 matrices at the same time. But in our case the diagonalization is quite bad.
toten=sum(arrdiag$D^2)#total energy
diagen=0
for (i in 1:12) {
  diagen=diagen+sum(diag(arrdiag$D[,,i])^2)
}
diagen#energy on every diagonal
offen=toten-diagen#energy outside the diagonals
offen/toten#almost half of the energy is stil outside the diagonals. The joint diagonalization tries to minimize
#this value, but some times the results will be bad just because of how the data are.
#Lets compare it to the result before the diagonalization
diagen=0
for (i in 1:12) {
  diagen=diagen+sum(diag(arr[,,i])^2)
}
diagen
(sum(arr^2)-diagen)/sum(arr^2)#before we had 82.8% of the energy outside the diagonals, so the diagonalization
#gave us a better situation
#The fact that the matrices autoc(tao) where not symmetric could have not let us have the best diagonalization

#Now I want to get the ten values that gives us the importance of each directions, to get them in SOBI we 
#combine the informations for each lag.
lambda=numeric(10)
for (i in 1:10) {
  lambda[i]=sum(arrdiag$D[i,i,]^2)
}
lambda
ord=order(lambda, decreasing = T)
lambdaord=numeric(10)
for (i in 1:10) {
  lambdaord[i]=lambda[ord[i]]
}
lambdaord #lambda ordered
#Now we can order the SOBI directions based on the associated values lambda
arrdiag$V
sobidir=matrix(rep(0,100), nrow = 10)
for (i in 1:10) {
  sobidir[,i]=arrdiag$V[,ord[i]]
}
sobidir #The SOBI directions are ordered
sobidir%*%t(sobidir) #still orthogonal

#definition of the function phi(k)
phi=numeric(10)
for (i in 1:10) {
  phi[i]=lambdaord[i]/(1+sum(lambdaord))
}
phi
plot(c(0,1,2,3,4,5,6,7,8,9),phi, type = 'l')
#we had from 1 to 10, but in the paper we have from 0 to 9: that's why I specified it.
#phi measurs how much second order serial dependence is associated to the sobi direction of rank k.
#An high phi(k) means that the directions presents a second order serial dependence relatively strong.
#We use it because the white noise variables doesn't have second order serial dependence (so we can recognize
#them also with phi).

#Now we need f(k), that evaluates the stability of the estimated subspaces.
#For the bootstrap I would choose a stattionary bootstrap with p=.005: in this way the resampled serie is 
#stationary, the temporal dependence is almost preserved (p isn't too big) and we will have also a good
#variability (p isn't too small): it's a good compromise between the two requests.
#I choose B=1000
B=1000
S2boot=array(0,c(10,10,1000))
sobidirboot=array(0,c(10,10,1000))
set.seed(123)
for (j in 1:B) {
  datab=funstat(datamat,0.05)
  means=apply(datab, 2, mean)
  datac=sweep(datab,2,means,'-') 
  datac
  apply(datac, 2, mean)
  S=cov(datac)
  sing=svd(S)
  D2=diag(sing$d^(-0.5)) 
  D2
  S2=sing$u%*%D2%*%t(sing$u)
  S2boot[,,j]=S2
  #whitening
  datac2=as.matrix(datac)
  dataw=datac2%*%S2 
  arr=array(0,c(10,10,12))
  for (i in 1:12) {
    arr[,,i]=autoc(i)
  }
  arrdiag=rjd(arr, maxiter = 300)
  lambda=numeric(10)
  for (i in 1:10) {
    lambda[i]=sum(arrdiag$D[i,i,]^2)
  }
  lambda
  ord=order(lambda, decreasing = T)
  lambdaord=numeric(10)
  for (i in 1:10) {
    lambdaord[i]=lambda[ord[i]]
  }
  lambdaord #lambda ordered
  #Now we can order the SOBI directions based on the associated values lambda
  arrdiag$V
  sobidir=matrix(rep(0,100), nrow = 10)
  for (i in 1:10) {
    sobidir[,i]=arrdiag$V[,ord[i]]
  }
  sobidirboot[,,j]=sobidir
}
S2boot
sobidirboot
#NOW RUN EVERYTHING ON TOP OF B=1000 JUST TO GET THE VALUES FOR THE ORIGINAL DATA
#Then
f_0=numeric(10)
detj=numeric(10)
for (k in 2:10) {
  B_k=sobidir[,1:(k-1), drop=F]
  for (j in 1:B) {
    B_kj=matrix(rep(0,10*(k-1)), nrow=10)
    B_kj=sobidirboot[,1:(k-1),j]
    detj[j]=1-abs(det(t(B_k)%*%B_kj))
  }
  f_0[k]=1/B*sum(detj)
}
f_0
plot(c(0,1,2,3,4,5,6,7,8,9),f_0, type = 'l')
#standardize to get f
f=numeric(10)
for (k in 1:10) {
  f[k]=f_0[k]/(1+sum(f_0))
}
f
plot(c(0,1,2,3,4,5,6,7,8,9),f, type = 'l')
#getting g
g=numeric(10)
for (k in 1:10) {
  g[k]=f[k]+phi[k]
}
g
plot(c(0,1,2,3,4,5,6,7,8,9),g, type = 'l')#we pick q=3

#point out the intermadiate checks.
#Then create simulate data to see if g works.



