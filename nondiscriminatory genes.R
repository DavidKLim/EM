source("C:/Users/David/Desktop/Research/Coding EM/EM.R")
library("MASS")

n=20
k=3
g=200
pi=c(.2,.3,.5)
b=matrix(rep(0,times=g*k),nrow=g)

sigma=diag(k)
b=matrix(rep(c(3,3.1,3.15),times=g),nrow=g,byrow=TRUE)

b[1:200,]<-mvrnorm(200,mu=c(3,3.1,3.15),sigma)

b[1:50,]<-mvrnorm(50,mu=c(3,3.5,4),sigma)
b[51:100,]<-mvrnorm(50,mu=c(1,1.2,6.5),sigma)
b[101:180,]<-mvrnorm(80,mu=c(2,4,5),sigma)
b[181:200,]<-mvrnorm(20,mu=c(1,1,1),sigma) # nondiscriminatory

X<-EM(n,k,g,pi,b)
