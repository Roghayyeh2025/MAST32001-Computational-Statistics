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
autoc=function(dataw,r){
  n=nrow(dataw)
  A=dataw[1:(n-r),]
  B=dataw[(1+r):n,] #lagging
  aut=(t(A)%*%B)/(n-r) #it's the formula sum(t=1 to n-r)(x_t^w)(x_(t+r)^w)^T/(n-r) (we use n-r because we just have n-r available temporal couples)
  aut
}
autoc(dataw,1) #it gives us the relation between X_i,t and X_j,t+1. It doesn't have to be symmetric.
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
  star=sample(1:981, size=1) #I didn't choose a circular moving block bootstrap
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
  arr[,,i]=autoc(dataw,i)
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
#I choose B=500
B=500
sobidirboot=array(0,c(10,10,B))
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
  #whitening
  datac2=as.matrix(datac)
  dataw=datac2%*%S2 
  arr=array(0,c(10,10,12))
  for (i in 1:12) {
    arr[,,i]=autoc(dataw,i)
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

sobidirboot
#NOW RUN EVERYTHING ON TOP OF B=1000 JUST TO GET THE VALUES FOR THE ORIGINAL DATA
#Then
f_0=numeric(10)
detj=numeric(B)
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
#Notice that this code could be useful for part E, while for this part we should use probably simulated data!

#
###PART D###
#
nobs=1000#my choice
set.seed(123)
x1=arima.sim(model = list(ar=0.6, ma=c(0,-0.4)), n=nobs)
#but the paper tells us that in the SOS model E(z_t)=0 and Cov(z_t)=I. For every variable x1,...,x5 we already
#have E(x)=0 and cov between every x =0, but we need V(x)=1.We use the theoretical variance
var1=sum(c(1,ARMAtoMA(ar=0.6,ma=c(0,-0.4),lag.max = 1000))^2)
x1=x1/sqrt(var1)
x2=arima.sim(model = list(ar=c(0.4,0.1,0.3), ma=c(0.2,0.4)), n=nobs)
var2=sum(c(1,ARMAtoMA(ar=c(0.4,0.1,0.3),ma=c(0.2,0.4),lag.max = 1000))^2)
x2=x2/sqrt(var2)
x3=arima.sim(model = list(ar=0.4), n=nobs)
var3=1/(1-0.4^2)
x3=x3/sqrt(var3)
x4=arima.sim(model = list(ma=c(0.5,0.3,0.1,0.2)), n=nobs)
var4=1+sum(c(0.5,0.3,0.1,0.2)^2)
x4=x4/sqrt(var4)
x5=arima.sim(model = list(ar=c(0.1,0.1,0.3)), n=nobs)
var5=sum(c(1,ARMAtoMA(ar=c(0.1,0.1,0.3),lag.max = 1000))^2)
x5=x5/sqrt(var5)
x6=rnorm(nobs)
x7=rnorm(nobs)
x8=rnorm(nobs)
x9=rnorm(nobs)
x10=rnorm(nobs)
setting1=cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
setting1
#We obtained a dataset using the right theoretical hyphotesis. Now we need whitening to work with the sample 
#informations.
#whitening
#centering
means=apply(setting1, 2, mean)
settingc=sweep(setting1,2,means,'-') 
apply(settingc, 2, mean)
#covariance matrix
S=cov(settingc) 
S
sing=svd(S)
D2=diag(sing$d^(-0.5)) #D^(-0.5)
D2
S2=sing$u%*%D2%*%t(sing$u)
S2 #S^(-0.5)
#whitening
settingc2=as.matrix(settingc)
settingw=settingc2%*%S2 #my scaling convention is the multivariate whitening based on the sample covariance matrix.
cov(settingw)
#autocorrelation
autocset=function(settingw,r){
  n=nrow(settingw)
  A=settingw[1:(n-r),]
  B=settingw[(1+r):n,] #lagging
  aut=(t(A)%*%B)/(n-r) #it's the formula sum(t=1 to n-r)(x_t^w)(x_(t+r)^w)^T/(n-r) (we use n-r because we just have n-r available temporal couples)
  aut
}
#calculation of phi
arrset=array(0,c(10,10,12))
for (i in 1:12) {
  arrset[,,i]=autocset(settingw,i)
}
arrset
#joint diagonalization
arrsetdiag=rjd(arrset)
arrsetdiag$V%*%t(arrsetdiag$V)#It's orthogonal: V*V^T=I. In the paper is called t(U_Tao)
arrsetdiag$D
#get the lambda for phi
lambdaset=numeric(10)
for (i in 1:10) {
  lambdaset[i]=sum(arrsetdiag$D[i,i,]^2)
}
lambdaset
ordset=order(lambdaset, decreasing = T)
lambdaordset=numeric(10)
for (i in 1:10) {
  lambdaordset[i]=lambdaset[ordset[i]]
}
lambdaordset #lambda ordered
#Now we can order the SOBI directions based on the associated values lambda
arrsetdiag$V
sobidirset=matrix(rep(0,100), nrow = 10)
for (i in 1:10) {
  sobidirset[,i]=arrsetdiag$V[,ordset[i]]
}
sobidirset #The SOBI directions are ordered
phiset=numeric(10)
for (i in 1:10) {
  phiset[i]=lambdaordset[i]/(1+sum(lambdaordset))
}
phiset
plot(c(0,1,2,3,4,5,6,7,8,9),phiset, type = 'l')
#now we are getting different value of f (and so of g) by using different bootstraps
#we compare stationary bootstrap with p=0.02,0.05,0.10 and block bootstrap with l=50,20,10.
B=10
g=array(0,c(10,10,6))
sobidirboot=array(0,c(10,10,B))

set.seed(123)

for (i in 1:10) { #number of Monte Carlo repetitions
  
  nobs=1000
  
  x1=arima.sim(model=list(ar=0.6,ma=c(0,-0.4)),n=nobs)
  x1=x1/sqrt(var1)
  
  x2=arima.sim(model=list(ar=c(0.4,0.1,0.3),ma=c(0.2,0.4)),n=nobs)
  x2=x2/sqrt(var2)
  
  x3=arima.sim(model=list(ar=0.4),n=nobs)
  x3=x3/sqrt(var3)
  
  x4=arima.sim(model=list(ma=c(0.5,0.3,0.1,0.2)),n=nobs)
  x4=x4/sqrt(var4)
  
  x5=arima.sim(model=list(ar=c(0.1,0.1,0.3)),n=nobs)
  x5=x5/sqrt(var5)
  
  x6=rnorm(nobs)
  x7=rnorm(nobs)
  x8=rnorm(nobs)
  x9=rnorm(nobs)
  x10=rnorm(nobs)
  
  setting=cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
  
  
  # =====================================================
  # ORIGINAL DATASET: whitening, SOBI, phiset, sobidirset
  # =====================================================
  
  means=apply(setting,2,mean)
  settingc=sweep(setting,2,means,'-')
  
  S=cov(settingc)
  sing=svd(S)
  D2=diag(sing$d^(-0.5))
  S2=sing$u%*%D2%*%t(sing$u)
  
  settingw=as.matrix(settingc)%*%S2
  
  arrset=array(0,c(10,10,12))
  
  for (z in 1:12) {
    arrset[,,z]=autoc(settingw,z)
  }
  
  arrsetdiag=rjd(arrset,maxiter=300)
  
  lambdaset=numeric(10)
  
  for (z in 1:10) {
    lambdaset[z]=sum(arrsetdiag$D[z,z,]^2)
  }
  
  ordset=order(lambdaset,decreasing=TRUE)
  
  lambdaordset=numeric(10)
  
  for (z in 1:10) {
    lambdaordset[z]=lambdaset[ordset[z]]
  }
  
  # Original SOBI directions: never overwritten by bootstrap
  sobidirset=matrix(0,nrow=10,ncol=10)
  
  for (z in 1:10) {
    sobidirset[,z]=arrsetdiag$V[,ordset[z]]
  }
  
  # phiset is recalculated for every Monte Carlo repetition
  phiset=numeric(10)
  
  for (z in 1:10) {
    phiset[z]=lambdaordset[z]/(1+sum(lambdaordset))
  }
  
  
  # =====================================================
  # METHOD 1: stationary bootstrap, p=0.02
  # =====================================================
  
  for (j in 1:B) {
    
    datab=funstat(setting,0.02)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,'-')
    
    S=cov(datac)
    sing=svd(S)
    D2=diag(sing$d^(-0.5))
    S2=sing$u%*%D2%*%t(sing$u)
    
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for (z in 1:12) {
      arr[,,z]=autoc(dataw,z)
    }
    
    arrdiag=rjd(arr,maxiter=300)
    
    lambda=numeric(10)
    
    for (z in 1:10) {
      lambda[z]=sum(arrdiag$D[z,z,]^2)
    }
    
    ord=order(lambda,decreasing=TRUE)
    
    sobidir=matrix(0,nrow=10,ncol=10)
    
    for (z in 1:10) {
      sobidir[,z]=arrdiag$V[,ord[z]]
    }
    
    sobidirboot[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for (k in 2:10) {
    
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for (j in 1:B) {
      B_kj=sobidirboot[,1:(k-1),j,drop=FALSE]
      B_kj=matrix(B_kj,nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=numeric(10)
  
  for (k in 1:10) {
    f[k]=f_0[k]/(1+sum(f_0))
  }
  
  for (k in 1:10) {
    g[i,k,1]=f[k]+phiset[k]
  }
  
  
  # =====================================================
  # METHOD 2: stationary bootstrap, p=0.05
  # =====================================================
  
  for (j in 1:B) {
    
    datab=funstat(setting,0.05)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,'-')
    
    S=cov(datac)
    sing=svd(S)
    D2=diag(sing$d^(-0.5))
    S2=sing$u%*%D2%*%t(sing$u)
    
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for (z in 1:12) {
      arr[,,z]=autoc(dataw,z)
    }
    
    arrdiag=rjd(arr,maxiter=300)
    
    lambda=numeric(10)
    
    for (z in 1:10) {
      lambda[z]=sum(arrdiag$D[z,z,]^2)
    }
    
    ord=order(lambda,decreasing=TRUE)
    
    sobidir=matrix(0,nrow=10,ncol=10)
    
    for (z in 1:10) {
      sobidir[,z]=arrdiag$V[,ord[z]]
    }
    
    sobidirboot[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for (k in 2:10) {
    
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for (j in 1:B) {
      B_kj=sobidirboot[,1:(k-1),j,drop=FALSE]
      B_kj=matrix(B_kj,nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=numeric(10)
  
  for (k in 1:10) {
    f[k]=f_0[k]/(1+sum(f_0))
  }
  
  for (k in 1:10) {
    g[i,k,2]=f[k]+phiset[k]
  }
  
  
  # =====================================================
  # METHOD 3: stationary bootstrap, p=0.1
  # =====================================================
  
  for (j in 1:B) {
    
    datab=funstat(setting,0.1)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,'-')
    
    S=cov(datac)
    sing=svd(S)
    D2=diag(sing$d^(-0.5))
    S2=sing$u%*%D2%*%t(sing$u)
    
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for (z in 1:12) {
      arr[,,z]=autoc(dataw,z)
    }
    
    arrdiag=rjd(arr,maxiter=300)
    
    lambda=numeric(10)
    
    for (z in 1:10) {
      lambda[z]=sum(arrdiag$D[z,z,]^2)
    }
    
    ord=order(lambda,decreasing=TRUE)
    
    sobidir=matrix(0,nrow=10,ncol=10)
    
    for (z in 1:10) {
      sobidir[,z]=arrdiag$V[,ord[z]]
    }
    
    sobidirboot[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for (k in 2:10) {
    
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for (j in 1:B) {
      B_kj=sobidirboot[,1:(k-1),j,drop=FALSE]
      B_kj=matrix(B_kj,nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=numeric(10)
  
  for (k in 1:10) {
    f[k]=f_0[k]/(1+sum(f_0))
  }
  
  for (k in 1:10) {
    g[i,k,3]=f[k]+phiset[k]
  }
  
  
  # =====================================================
  # METHOD 4: fixed block bootstrap, l=50
  # =====================================================
  
  for (j in 1:B) {
    
    datab=funblock(setting,50)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,'-')
    
    S=cov(datac)
    sing=svd(S)
    D2=diag(sing$d^(-0.5))
    S2=sing$u%*%D2%*%t(sing$u)
    
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for (z in 1:12) {
      arr[,,z]=autoc(dataw,z)
    }
    
    arrdiag=rjd(arr,maxiter=300)
    
    lambda=numeric(10)
    
    for (z in 1:10) {
      lambda[z]=sum(arrdiag$D[z,z,]^2)
    }
    
    ord=order(lambda,decreasing=TRUE)
    
    sobidir=matrix(0,nrow=10,ncol=10)
    
    for (z in 1:10) {
      sobidir[,z]=arrdiag$V[,ord[z]]
    }
    
    sobidirboot[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for (k in 2:10) {
    
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for (j in 1:B) {
      B_kj=sobidirboot[,1:(k-1),j,drop=FALSE]
      B_kj=matrix(B_kj,nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=numeric(10)
  
  for (k in 1:10) {
    f[k]=f_0[k]/(1+sum(f_0))
  }
  
  for (k in 1:10) {
    g[i,k,4]=f[k]+phiset[k]
  }
  
  
  # =====================================================
  # METHOD 5: fixed block bootstrap, l=20
  # =====================================================
  
  for (j in 1:B) {
    
    datab=funblock(setting,20)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,'-')
    
    S=cov(datac)
    sing=svd(S)
    D2=diag(sing$d^(-0.5))
    S2=sing$u%*%D2%*%t(sing$u)
    
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for (z in 1:12) {
      arr[,,z]=autoc(dataw,z)
    }
    
    arrdiag=rjd(arr,maxiter=300)
    
    lambda=numeric(10)
    
    for (z in 1:10) {
      lambda[z]=sum(arrdiag$D[z,z,]^2)
    }
    
    ord=order(lambda,decreasing=TRUE)
    
    sobidir=matrix(0,nrow=10,ncol=10)
    
    for (z in 1:10) {
      sobidir[,z]=arrdiag$V[,ord[z]]
    }
    
    sobidirboot[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for (k in 2:10) {
    
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for (j in 1:B) {
      B_kj=sobidirboot[,1:(k-1),j,drop=FALSE]
      B_kj=matrix(B_kj,nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=numeric(10)
  
  for (k in 1:10) {
    f[k]=f_0[k]/(1+sum(f_0))
  }
  
  for (k in 1:10) {
    g[i,k,5]=f[k]+phiset[k]
  }
  
  
  # =====================================================
  # METHOD 6: fixed block bootstrap, l=10
  # =====================================================
  
  for (j in 1:B) {
    
    datab=funblock(setting,10)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,'-')
    
    S=cov(datac)
    sing=svd(S)
    D2=diag(sing$d^(-0.5))
    S2=sing$u%*%D2%*%t(sing$u)
    
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for (z in 1:12) {
      arr[,,z]=autoc(dataw,z)
    }
    
    arrdiag=rjd(arr,maxiter=300)
    
    lambda=numeric(10)
    
    for (z in 1:10) {
      lambda[z]=sum(arrdiag$D[z,z,]^2)
    }
    
    ord=order(lambda,decreasing=TRUE)
    
    sobidir=matrix(0,nrow=10,ncol=10)
    
    for (z in 1:10) {
      sobidir[,z]=arrdiag$V[,ord[z]]
    }
    
    sobidirboot[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for (k in 2:10) {
    
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for (j in 1:B) {
      B_kj=sobidirboot[,1:(k-1),j,drop=FALSE]
      B_kj=matrix(B_kj,nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=numeric(10)
  
  for (k in 1:10) {
    f[k]=f_0[k]/(1+sum(f_0))
  }
  
  for (k in 1:10) {
    g[i,k,6]=f[k]+phiset[k]
  }
  
}
g
result=matrix(rep(0,10*6),nrow=6)
for (j in 1:6) {
  for (i in 1:10) {
    result[j,i]=which.min(g[i,,j])-1
  }
}
result




####CHECKS###

####NEW TRY###
library(JADE)

B=200
nobs=1000

g=array(0,c(100,10,6))
sobidirboot=array(0,c(10,10,B))

set.seed(123)

for(i in 1:100){
  
  # Generazione dei segnali
  s1=as.numeric(arima.sim(model=list(ar=.6,ma=c(0,-.4)),n=nobs))/sqrt(var1)
  s2=as.numeric(arima.sim(model=list(ar=c(.4,.1,.3),ma=c(.2,.4)),n=nobs))/sqrt(var2)
  s3=as.numeric(arima.sim(model=list(ar=.4),n=nobs))/sqrt(var3)
  s4=as.numeric(arima.sim(model=list(ma=c(.5,.3,.1,.2)),n=nobs))/sqrt(var4)
  s5=as.numeric(arima.sim(model=list(ar=c(.1,.1,.3)),n=nobs))/sqrt(var5)
  
  setting=cbind(
    s1,s2,s3,s4,s5,
    rnorm(nobs),rnorm(nobs),rnorm(nobs),rnorm(nobs),rnorm(nobs)
  )
  
  # SOBI sul campione originale
  means=apply(setting,2,mean)
  settingc=sweep(setting,2,means,"-")
  S=cov(settingc)
  sing=svd(S)
  S2=sing$u%*%diag(sing$d^(-.5))%*%t(sing$u)
  settingw=as.matrix(settingc)%*%S2
  
  arrset=array(0,c(10,10,12))
  
  for(z in 1:12){
    A=autoc(settingw,z)
    arrset[,,z]=(A+t(A))/2
  }
  
  arrsetdiag=frjd(arrset,maxiter=300,eps=1e-3)
  
  lambdaset=numeric(10)
  
  for(z in 1:12){
    lambdaset=lambdaset+diag(t(arrsetdiag$V)%*%arrset[,,z]%*%arrsetdiag$V)^2
  }
  
  ordset=order(lambdaset,decreasing=TRUE)
  sobidirset=arrsetdiag$V[,ordset,drop=FALSE]
  lambdaset=lambdaset[ordset]
  phiset=lambdaset/(1+sum(lambdaset))
  
  
  # ========================================
  # 1. STATIONARY BOOTSTRAP p=0.02
  # ========================================
  
  for(j in 1:B){
    
    datab=funstat(setting,.02)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,"-")
    S=cov(datac)
    sing=svd(S)
    S2=sing$u%*%diag(sing$d^(-.5))%*%t(sing$u)
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for(z in 1:12){
      A=autoc(dataw,z)
      arr[,,z]=(A+t(A))/2
    }
    
    arrdiag=frjd(arr,maxiter=300,eps=1e-3)
    
    lambda=numeric(10)
    
    for(z in 1:12){
      lambda=lambda+diag(t(arrdiag$V)%*%arr[,,z]%*%arrdiag$V)^2
    }
    
    ord=order(lambda,decreasing=TRUE)
    sobidir=arrdiag$V[,ord,drop=FALSE]
    sobidirboot[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for(k in 2:10){
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for(j in 1:B){
      B_kj=matrix(sobidirboot[,1:(k-1),j],nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=f_0/(1+sum(f_0))
  
  for(k in 1:10){
    g[i,k,1]=f[k]+phiset[k]
  }
  
  
  # ========================================
  # 2. STATIONARY BOOTSTRAP p=0.05
  # ========================================
  
  for(j in 1:B){
    
    datab=funstat(setting,.05)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,"-")
    S=cov(datac)
    sing=svd(S)
    S2=sing$u%*%diag(sing$d^(-.5))%*%t(sing$u)
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for(z in 1:12){
      A=autoc(dataw,z)
      arr[,,z]=(A+t(A))/2
    }
    
    arrdiag=frjd(arr,maxiter=300,eps=1e-3)
    
    lambda=numeric(10)
    
    for(z in 1:12){
      lambda=lambda+diag(t(arrdiag$V)%*%arr[,,z]%*%arrdiag$V)^2
    }
    
    ord=order(lambda,decreasing=TRUE)
    sobidir=arrdiag$V[,ord,drop=FALSE]
    sobidirboot[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for(k in 2:10){
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for(j in 1:B){
      B_kj=matrix(sobidirboot[,1:(k-1),j],nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=f_0/(1+sum(f_0))
  
  for(k in 1:10){
    g[i,k,2]=f[k]+phiset[k]
  }
  
  
  # ========================================
  # 3. STATIONARY BOOTSTRAP p=0.1
  # ========================================
  
  for(j in 1:B){
    
    datab=funstat(setting,.1)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,"-")
    S=cov(datac)
    sing=svd(S)
    S2=sing$u%*%diag(sing$d^(-.5))%*%t(sing$u)
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for(z in 1:12){
      A=autoc(dataw,z)
      arr[,,z]=(A+t(A))/2
    }
    
    arrdiag=frjd(arr,maxiter=300,eps=1e-3)
    
    lambda=numeric(10)
    
    for(z in 1:12){
      lambda=lambda+diag(t(arrdiag$V)%*%arr[,,z]%*%arrdiag$V)^2
    }
    
    ord=order(lambda,decreasing=TRUE)
    sobidir=arrdiag$V[,ord,drop=FALSE]
    sobidirboot[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for(k in 2:10){
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for(j in 1:B){
      B_kj=matrix(sobidirboot[,1:(k-1),j],nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=f_0/(1+sum(f_0))
  
  for(k in 1:10){
    g[i,k,3]=f[k]+phiset[k]
  }
  
  
  # ========================================
  # 4. FIXED BLOCK BOOTSTRAP l=50
  # ========================================
  
  for(j in 1:B){
    
    datab=funblock(setting,50)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,"-")
    S=cov(datac)
    sing=svd(S)
    S2=sing$u%*%diag(sing$d^(-.5))%*%t(sing$u)
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for(z in 1:12){
      A=autoc(dataw,z)
      arr[,,z]=(A+t(A))/2
    }
    
    arrdiag=frjd(arr,maxiter=300,eps=1e-3)
    
    lambda=numeric(10)
    
    for(z in 1:12){
      lambda=lambda+diag(t(arrdiag$V)%*%arr[,,z]%*%arrdiag$V)^2
    }
    
    ord=order(lambda,decreasing=TRUE)
    sobidir=arrdiag$V[,ord,drop=FALSE]
    sobidirboot[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for(k in 2:10){
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for(j in 1:B){
      B_kj=matrix(sobidirboot[,1:(k-1),j],nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=f_0/(1+sum(f_0))
  
  for(k in 1:10){
    g[i,k,4]=f[k]+phiset[k]
  }
  
  
  # ========================================
  # 5. FIXED BLOCK BOOTSTRAP l=20
  # ========================================
  
  for(j in 1:B){
    
    datab=funblock(setting,20)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,"-")
    S=cov(datac)
    sing=svd(S)
    S2=sing$u%*%diag(sing$d^(-.5))%*%t(sing$u)
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for(z in 1:12){
      A=autoc(dataw,z)
      arr[,,z]=(A+t(A))/2
    }
    
    arrdiag=frjd(arr,maxiter=300,eps=1e-3)
    
    lambda=numeric(10)
    
    for(z in 1:12){
      lambda=lambda+diag(t(arrdiag$V)%*%arr[,,z]%*%arrdiag$V)^2
    }
    
    ord=order(lambda,decreasing=TRUE)
    sobidir=arrdiag$V[,ord,drop=FALSE]
    sobidirboot[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for(k in 2:10){
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for(j in 1:B){
      B_kj=matrix(sobidirboot[,1:(k-1),j],nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=f_0/(1+sum(f_0))
  
  for(k in 1:10){
    g[i,k,5]=f[k]+phiset[k]
  }
  
  
  # ========================================
  # 6. FIXED BLOCK BOOTSTRAP l=10
  # ========================================
  
  for(j in 1:B){
    
    datab=funblock(setting,10)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,"-")
    S=cov(datac)
    sing=svd(S)
    S2=sing$u%*%diag(sing$d^(-.5))%*%t(sing$u)
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for(z in 1:12){
      A=autoc(dataw,z)
      arr[,,z]=(A+t(A))/2
    }
    
    arrdiag=frjd(arr,maxiter=300,eps=1e-3)
    
    lambda=numeric(10)
    
    for(z in 1:12){
      lambda=lambda+diag(t(arrdiag$V)%*%arr[,,z]%*%arrdiag$V)^2
    }
    
    ord=order(lambda,decreasing=TRUE)
    sobidir=arrdiag$V[,ord,drop=FALSE]
    sobidirboot[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for(k in 2:10){
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for(j in 1:B){
      B_kj=matrix(sobidirboot[,1:(k-1),j],nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=f_0/(1+sum(f_0))
  
  for(k in 1:10){
    g[i,k,6]=f[k]+phiset[k]
  }
  
  cat("Monte Carlo:",i,"di 100 completato\n")
}

# Dimensioni stimate
qhat=matrix(0,100,6)

for(i in 1:100){
  for(j in 1:6){
    qhat[i,j]=which.min(g[i,,j])-1
  }
}

colnames(qhat)=c(
  "stat_0.02","stat_0.05","stat_0.1",
  "block_50","block_20","block_10"
)

print(qhat) #FINALLY RESULTS FOR N=1000, N OF MC=100 AND B=200


# Frequences
freq=apply(qhat,2,function(x){
  table(factor(x,levels=0:9))
})
print(freq) 
#we have 100/100 choice of 5 for p=0.02 and 0.05. We have 99/100 choice of 5 and 1/100 of 4 for p=0.1,l=50 and
#l=20. We have 97/100 choice of 5 and 3/100 choice of 4 for l=10
p=freq[6,]/100#estimator \hat(q) accuracy. Since we know that the real q is 5
#p is 1,1,0.99,0.99,0.99,0.97
#uncertainty
sqrt(p*(1-p)/100)# MC MCSE: 0 0 0.00990.0099 0.0099 0.0171
#For methods achieving 100% accuracy, the estimated Monte Carlo standard error is zero because no failures 
#were observed in the 100 replications. However, this does not imply the absence of Monte Carlo uncertainty, 
#as the true success probability may still be below one.
#variability
mean_q = colMeans(qhat)
sd_q = apply(qhat, 2, sd)
sd_q#0,0,0.1,0.1,0.171
#We have low sd, infact the different replies give often (almost always) the same result



#SAME FOR N=200
library(JADE)

B=200
nobs=500

g2=array(0,c(100,10,6))
sobidirboot2=array(0,c(10,10,B))

set.seed(123)

for(i in 1:100){
  
  # Generazione dei segnali
  s1=as.numeric(arima.sim(model=list(ar=.6,ma=c(0,-.4)),n=nobs))/sqrt(var1)
  s2=as.numeric(arima.sim(model=list(ar=c(.4,.1,.3),ma=c(.2,.4)),n=nobs))/sqrt(var2)
  s3=as.numeric(arima.sim(model=list(ar=.4),n=nobs))/sqrt(var3)
  s4=as.numeric(arima.sim(model=list(ma=c(.5,.3,.1,.2)),n=nobs))/sqrt(var4)
  s5=as.numeric(arima.sim(model=list(ar=c(.1,.1,.3)),n=nobs))/sqrt(var5)
  
  setting=cbind(
    s1,s2,s3,s4,s5,
    rnorm(nobs),rnorm(nobs),rnorm(nobs),rnorm(nobs),rnorm(nobs)
  )
  
  # SOBI sul campione originale
  means=apply(setting,2,mean)
  settingc=sweep(setting,2,means,"-")
  S=cov(settingc)
  sing=svd(S)
  S2=sing$u%*%diag(sing$d^(-.5))%*%t(sing$u)
  settingw=as.matrix(settingc)%*%S2
  
  arrset=array(0,c(10,10,12))
  
  for(z in 1:12){
    A=autoc(settingw,z)
    arrset[,,z]=(A+t(A))/2
  }
  
  arrsetdiag=frjd(arrset,maxiter=300,eps=1e-3)
  
  lambdaset=numeric(10)
  
  for(z in 1:12){
    lambdaset=lambdaset+
      diag(t(arrsetdiag$V)%*%arrset[,,z]%*%arrsetdiag$V)^2
  }
  
  ordset=order(lambdaset,decreasing=TRUE)
  sobidirset=arrsetdiag$V[,ordset,drop=FALSE]
  lambdaset=lambdaset[ordset]
  phiset=lambdaset/(1+sum(lambdaset))
  
  
  # ========================================
  # 1. STATIONARY BOOTSTRAP p=0.02
  # ========================================
  
  for(j in 1:B){
    
    datab=funstat(setting,.02)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,"-")
    S=cov(datac)
    sing=svd(S)
    S2=sing$u%*%diag(sing$d^(-.5))%*%t(sing$u)
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for(z in 1:12){
      A=autoc(dataw,z)
      arr[,,z]=(A+t(A))/2
    }
    
    arrdiag=frjd(arr,maxiter=300,eps=1e-3)
    
    lambda=numeric(10)
    
    for(z in 1:12){
      lambda=lambda+
        diag(t(arrdiag$V)%*%arr[,,z]%*%arrdiag$V)^2
    }
    
    ord=order(lambda,decreasing=TRUE)
    sobidir=arrdiag$V[,ord,drop=FALSE]
    sobidirboot2[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for(k in 2:10){
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for(j in 1:B){
      B_kj=matrix(sobidirboot2[,1:(k-1),j],nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=f_0/(1+sum(f_0))
  
  for(k in 1:10){
    g2[i,k,1]=f[k]+phiset[k]
  }
  
  
  # ========================================
  # 2. STATIONARY BOOTSTRAP p=0.05
  # ======================================== 
  
  for(j in 1:B){
    
    datab=funstat(setting,.05)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,"-")
    S=cov(datac)
    sing=svd(S)
    S2=sing$u%*%diag(sing$d^(-.5))%*%t(sing$u)
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for(z in 1:12){
      A=autoc(dataw,z)
      arr[,,z]=(A+t(A))/2
    }
    
    arrdiag=frjd(arr,maxiter=300,eps=1e-3)
    
    lambda=numeric(10)
    
    for(z in 1:12){
      lambda=lambda+
        diag(t(arrdiag$V)%*%arr[,,z]%*%arrdiag$V)^2
    }
    
    ord=order(lambda,decreasing=TRUE)
    sobidir=arrdiag$V[,ord,drop=FALSE]
    sobidirboot2[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for(k in 2:10){
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for(j in 1:B){
      B_kj=matrix(sobidirboot2[,1:(k-1),j],nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=f_0/(1+sum(f_0))
  
  for(k in 1:10){
    g2[i,k,2]=f[k]+phiset[k]
  }
  
  
  # ========================================
  # 3. STATIONARY BOOTSTRAP p=0.1
  # ========================================
  
  for(j in 1:B){
    
    datab=funstat(setting,.1)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,"-")
    S=cov(datac)
    sing=svd(S)
    S2=sing$u%*%diag(sing$d^(-.5))%*%t(sing$u)
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for(z in 1:12){
      A=autoc(dataw,z)
      arr[,,z]=(A+t(A))/2
    }
    
    arrdiag=frjd(arr,maxiter=300,eps=1e-3)
    
    lambda=numeric(10)
    
    for(z in 1:12){
      lambda=lambda+
        diag(t(arrdiag$V)%*%arr[,,z]%*%arrdiag$V)^2
    }
    
    ord=order(lambda,decreasing=TRUE)
    sobidir=arrdiag$V[,ord,drop=FALSE]
    sobidirboot2[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for(k in 2:10){
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for(j in 1:B){
      B_kj=matrix(sobidirboot2[,1:(k-1),j],nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=f_0/(1+sum(f_0))
  
  for(k in 1:10){
    g2[i,k,3]=f[k]+phiset[k]
  }
  
  
  # ========================================
  # 4. FIXED BLOCK BOOTSTRAP l=50
  # ========================================
  
  for(j in 1:B){
    
    datab=funblock(setting,50)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,"-")
    S=cov(datac)
    sing=svd(S)
    S2=sing$u%*%diag(sing$d^(-.5))%*%t(sing$u)
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for(z in 1:12){
      A=autoc(dataw,z)
      arr[,,z]=(A+t(A))/2
    }
    
    arrdiag=frjd(arr,maxiter=300,eps=1e-3)
    
    lambda=numeric(10)
    
    for(z in 1:12){
      lambda=lambda+
        diag(t(arrdiag$V)%*%arr[,,z]%*%arrdiag$V)^2
    }
    
    ord=order(lambda,decreasing=TRUE)
    sobidir=arrdiag$V[,ord,drop=FALSE]
    sobidirboot2[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for(k in 2:10){
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for(j in 1:B){
      B_kj=matrix(sobidirboot2[,1:(k-1),j],nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=f_0/(1+sum(f_0))
  
  for(k in 1:10){
    g2[i,k,4]=f[k]+phiset[k]
  }
  
  
  # ========================================
  # 5. FIXED BLOCK BOOTSTRAP l=20
  # ========================================
  
  for(j in 1:B){
    
    datab=funblock(setting,20)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,"-")
    S=cov(datac)
    sing=svd(S)
    S2=sing$u%*%diag(sing$d^(-.5))%*%t(sing$u)
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for(z in 1:12){
      A=autoc(dataw,z)
      arr[,,z]=(A+t(A))/2
    }
    
    arrdiag=frjd(arr,maxiter=300,eps=1e-3)
    
    lambda=numeric(10)
    
    for(z in 1:12){
      lambda=lambda+
        diag(t(arrdiag$V)%*%arr[,,z]%*%arrdiag$V)^2
    }
    
    ord=order(lambda,decreasing=TRUE)
    sobidir=arrdiag$V[,ord,drop=FALSE]
    sobidirboot2[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for(k in 2:10){
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for(j in 1:B){
      B_kj=matrix(sobidirboot2[,1:(k-1),j],nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=f_0/(1+sum(f_0))
  
  for(k in 1:10){
    g2[i,k,5]=f[k]+phiset[k]
  }
  
  
  # ========================================
  # 6. FIXED BLOCK BOOTSTRAP l=10
  # ========================================
  
  for(j in 1:B){
    
    datab=funblock(setting,10)
    
    means=apply(datab,2,mean)
    datac=sweep(datab,2,means,"-")
    S=cov(datac)
    sing=svd(S)
    S2=sing$u%*%diag(sing$d^(-.5))%*%t(sing$u)
    dataw=as.matrix(datac)%*%S2
    
    arr=array(0,c(10,10,12))
    
    for(z in 1:12){
      A=autoc(dataw,z)
      arr[,,z]=(A+t(A))/2
    }
    
    arrdiag=frjd(arr,maxiter=300,eps=1e-3)
    
    lambda=numeric(10)
    
    for(z in 1:12){
      lambda=lambda+
        diag(t(arrdiag$V)%*%arr[,,z]%*%arrdiag$V)^2
    }
    
    ord=order(lambda,decreasing=TRUE)
    sobidir=arrdiag$V[,ord,drop=FALSE]
    sobidirboot2[,,j]=sobidir
  }
  
  f_0=numeric(10)
  
  for(k in 2:10){
    B_k=sobidirset[,1:(k-1),drop=FALSE]
    detj=numeric(B)
    
    for(j in 1:B){
      B_kj=matrix(sobidirboot2[,1:(k-1),j],nrow=10,ncol=k-1)
      detj[j]=1-abs(det(t(B_k)%*%B_kj))
    }
    
    f_0[k]=sum(detj)/B
  }
  
  f=f_0/(1+sum(f_0))
  
  for(k in 1:10){
    g2[i,k,6]=f[k]+phiset[k]
  }
  
  cat("Monte Carlo:",i,"di 100 completato\n")
}

qhat2=matrix(0,100,6)

for(i in 1:100){
  for(j in 1:6){
    qhat2[i,j]=which.min(g2[i,,j])-1
  }
}

colnames(qhat2)=c(
  "stat_0.02","stat_0.05","stat_0.1",
  "block_50","block_20","block_10"
)

print(qhat2)
#analysis
freq2=apply(qhat2,2,function(x){
  table(factor(x,levels=0:9))
})
print(freq2)
p2=freq2[6,]/100#estimator \hat(q) accuracy. Since we know that the real q is 5
#p is 0.47-0.41-0.34-0.44-0.42-0.35
#uncertainty
sqrt(p2*(1-p2)/100)# MC MCSE: 0.0499 0.04918 0.0474 0.0496 0.0493 0.0477
#variability
mean_q2 = colMeans(qhat2)
sd_q2 = apply(qhat2, 2, sd)
sd_q2#1.851999  1.875170  1.873014  1.874092  1.876462  1.752891
#The joint diagonalization algorithm was run with a convergence tolerance of 10^(-3) and a maximum of 300 iterations
#This choice was made to have a lower computational cost.











