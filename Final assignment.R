###EXERCIZE 1###
library(readr)
timeseries_dat <- read_table("Downloads/timeseries.dat.txt", 
                             col_names = FALSE)
View(timeseries_dat)

###PART A###
#centering
data=timeseries_dat
means=apply(data, 2, mean)
datac=sweep(data,2,means,'-') #centered data
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
dataw=datac2%*%S2
cov(dataw)#We obtain almost I_10, so the whitening was done successfully
#Now we have whitened variables: they aren't correlated and they have variance=1 and mean=0. We still cannot
#say that they are indipendent between each other, and of course we didn't eliminate the temporal autocorrelation.
#why do we use whitening?
autoc=function(r){
  A=dataw[1:(1000-r),]
  B=dataw[(1+r):1000,]
  aut=(t(A)%*%B)/(1000-r) #it's the formula sum(t=1 to n-r)(x_t^w)(x_(t+r)^w)^T/(n-r)
  aut
}
autoc(1) #it gives us the relation between X_i,t and X_j,t+1. It doesn't have to be symmetric.
#The whitening didn't remove this correlation.
#I still didn't choose the relevant autocovariance matrices.





