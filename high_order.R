rm(list=ls())
rm(list = setdiff(ls(), lsf.str()))
library(survival)
library(KMsurv)
library(MASS)
set.seed(1234)

####################Adaptive Assignment#################
assign.ps<-function(x,n,p){
  n.temp=n+rbind(x,x)
  g1=abs(n.temp[1,]-n[2,])%*%x
  g2=abs(n.temp[2,]-n[1,])%*%x
  if (g1<g2) t<-rbinom(1,1,p)
  else if (g1==g2) t<-rbinom(1,1,.5)
  else t<-1-rbinom(1,1,p)
  return(t)
}

generate.I<-function(n,treatment){
  l <- list(6) 
  l[[1]]<-c(1,1,0,0)
  l[[2]]<-c(1,0,1,0)
  l[[3]]<-c(1,0,0,1)
  l[[4]]<-c(0,1,1,0)
  l[[5]]<-c(0,1,0,1)
  l[[6]]<-c(0,0,1,1)
  block<-sample(1:6, 1)
  if (n%%4>0) {
    I1<-c(treatment[1:(n-n%%4)],l[[block]][1:(n%%4)])
  } else {
    I1<-c(treatment[1:(n-n%%4)])
  }
  return(I1)
}

assign.spb<-function(z1,z2){
  ID<-seq(1:length(z1))
  block1<-rep(c(1,1,0,0),50)
  n.strata<-rep(0,4)
  n<-length(z1)
  u<-runif(n,0,1)
  z<-data.frame(z1,z2)
  ones<-rep(1,n)
  n.strata<-aggregate( ones, by = as.list(z), FUN = sum)$x 
  z<-data.frame(ID,z1,z2,u)
  zz<-z[order(z2,z1,u),]
  I<-c(generate.I(n.strata[1],block1),generate.I(n.strata[2],block1),generate.I(n.strata[3],block1),generate.I(n.strata[4],block1))
  #I<-c(treatment[1:n.strata[1]],treatment[1:n.strata[2]],treatment[1:n.strata[3]],treatment[1:n.strata[4]])
  X<-data.frame(zz$ID,zz$z1,zz$z2,I)
  names(X)<-c("ID","z1","z2","I")
  Y<-X[order(X$ID),]
  return(Y)
}

l1<-c(1,-1,0,0)
l2<-c(1,-1)
l3<-c(1,-1,0)
l4<-c(1,-1,0)
gen.weights<-function(data){
  ####  n<-nrow(data)
  N <- nrow(data)
  ordereddata <- data[order(data$y),]
  ordereddata$j <- seq(1:N)
  ordereddata$p <- ((N-ordereddata$j)/(N-ordereddata$j+1))**ordereddata$delta
  ordereddata$cumprod<-cumprod(ordereddata$p)
  w<-rep(0,N)
  w[1] <- ordereddata$delta[1]/N
  for (k in 2:N){ 
    w[k]<-ordereddata$delta[k]/(N-k+1)*ordereddata$cumprod[k-1]
  }
  
  ordereddata$w <-w
  ordereddata$I1<-1-ordereddata$I
  return(ordereddata)
}



####################generate dataset####################
generate.cr<-function(n.sample,beta0,beta1,beta2,beta3,p1,p2,p3,t){
  err<-runif(n.sample,-0.5,0.5)
  c<-runif(n.sample,0,t)
  I<-rbinom(n.sample,1,p1)
  z1<-rbinom(n.sample,1,p2)-0.5
  #zz2<-rnorm(n.sample,0.25,1)
  zz2<-runif(n.sample,-p3,p3)
  z2<-sign(zz2>=0)
  log.t<-beta0*(1-I)+beta1*I+beta2*z1+beta3*(zz2^2)+err
  t<-exp(log.t)
  delta<-sign(t<=c)
  y=pmin(t,c)
  data<-data.frame(I,z1,z2,zz2,t,y,c,delta)
  return(data)
}
#data<-generate.cr(0.5,0.5,.5,.5,0.5,0.5,0.5)
#mean(data$delta)

generate.ps<-function(n.sample,beta0,beta1,beta2,beta3,p1,p2,t){
  err<-runif(n.sample,-0.5,0.5)
  c<-runif(n.sample,0,t)
  z1<-rbinom(n.sample,1,p1)
  zz2<-runif(n.sample,-p2,p2)
  z2<-sign(zz2>=0)
  n<-matrix(0,2,4)
  z<-cbind(z1,1-z1,z2,1-z2)
  I<-rep(0,n.sample)
  for (i in 1:n.sample){
    I[i]<-assign.ps(z[i,],n,0.75)
    if(I[i]) n[1,]<-n[1,]+z[i,]
    else n[2,]<-n[2,]+z[i,]
  }
  z1<-z1-0.5
  z2<-z2-0.5
  log.t<-beta0*(1-I)+beta1*I+beta2*z1+beta3*(zz2^2)+err
  t<-exp(log.t)
  delta<-sign(t<=c)
  y=pmin(t,c)
  data<-data.frame(I,z1,z2,zz2,t,y,c,delta)
  return(data)
}

generate.spb<-function(n.sample,beta0,beta1,beta2,beta3,p1,p2,t){
  err<-runif(n.sample,-0.5,0.5)
  c<-runif(n.sample,0,t)
  z1<-rbinom(n.sample,1,p1)
  zz2<-runif(n.sample,-p2,p2)
  z2<-sign(zz2>=0)
  data.temp<-assign.spb(z1,z2)
  z1<-data.temp$z1-0.5
  z2<-data.temp$z2-0.5
  I<-data.temp$I
  
  log.t<-beta0*(1-I)+beta1*I+beta2*z1+beta3*(zz2^2)+err
  t<-exp(log.t)
  delta<-sign(t<=c)
  y=pmin(t,c)
  data<-data.frame(I,z1,z2,zz2,t,y,c,delta)
  return(data)
}

###########################################################


power.cr<-function(n.sim,n.sample,tt,beta0,beta1){
  #xx.names <- c("beta.hat1","beta.hat2","beta.hat3","beta.hat4","type I Err","V","SD")
  xx<-vector("list", 7)
  #names(xx) <- xx.names
  reject<-rep(0,4)
  beta.hat1<-matrix(0,nr=n.sim,nc=4)
  beta.hat2<-matrix(0,nr=n.sim,nc=2)
  beta.hat3<-matrix(0,nr=n.sim,nc=3)
  beta.hat4<-matrix(0,nr=n.sim,nc=3)
  V1<-V2<-V3<-V4<-rep(0,n.sim)
  l.beta1<-l.beta2<-l.beta3<-l.beta4<-rep(0,n.sim)
  for (k in 1:n.sim){ 
    data<-generate.cr(n.sample,beta0,beta1,.5,.5,0.5,0.5,0.5,tt)
    ################Fit Models############################# 
    data1<-gen.weights(data)
    model1 <- lm(log(y) ~ I+I1+z1+zz2-1, weights=w,data=data1)
    model2 <- lm(log(y) ~ I+I1-1, weights=w,data=data1)
    model3 <- lm(log(y) ~ I+I1+z1-1, weights=w,data=data1)
    model4 <- lm(log(y) ~ I+I1+zz2-1, weights=w,data=data1)
    
    beta.hat1[k,] <- c(model1$coeff[1],model1$coeff[2],model1$coeff[3],model1$coeff[4])
    beta.hat2[k,] <- c(model2$coeff[1],model2$coeff[2])
    beta.hat3[k,] <- c(model3$coeff[1],model3$coeff[2],model3$coeff[3])
    beta.hat4[k,] <- c(model4$coeff[1],model4$coeff[2],model4$coeff[3])
    l.beta1[k]<-beta.hat1[k,]%*%l1
    l.beta2[k]<-beta.hat2[k,]%*%l2
    l.beta3[k]<-beta.hat3[k,]%*%l3
    l.beta4[k]<-beta.hat4[k,]%*%l4
    V1[k]<-sqrt(t(l1)%*%vcov(model1)%*%l1)
    V2[k]<-sqrt(t(l2)%*%vcov(model2)%*%l2)
    V3[k]<-sqrt(t(l3)%*%vcov(model3)%*%l3)
    V4[k]<-sqrt(t(l4)%*%vcov(model4)%*%l4)
    reject[1]<-reject[1]+sign(abs(l.beta1[k]/V1[k])>1.96)
    reject[2]<-reject[2]+sign(abs(l.beta2[k]/V2[k])>1.96)
    reject[3]<-reject[3]+sign(abs(l.beta3[k]/V3[k])>1.96)
    reject[4]<-reject[4]+sign(abs(l.beta4[k]/V4[k])>1.96)
  }
  V<-c(mean(V1),mean(V2),mean(V3),mean(V4))
  SD<-c(sd(l.beta1),sd(l.beta2),sd(l.beta3),sd(l.beta4))
  xx[[1]]<-colMeans(beta.hat1)
  xx[[2]]<-colMeans(beta.hat2)
  xx[[3]]<-colMeans(beta.hat3)
  xx[[4]]<-colMeans(beta.hat4)
  xx[[5]]<-reject/n.sim
  xx[[6]]<-V
  xx[[7]]<-SD
  return(xx)
}

power.ps<-function(n.sim,n.sample,tt,beta0,beta1){
  print(paste("running ps",n.sim,n.sample,tt,beta0,beta1, sep=" "))
  #xx.names <- c("beta.hat1","beta.hat2","beta.hat3","beta.hat4","type I Err","V","SD")
  xx<-vector("list", 7)
  #names(xx) <- xx.names
  reject<-rep(0,4)
  beta.hat1<-matrix(0,nr=n.sim,nc=4)
  beta.hat2<-matrix(0,nr=n.sim,nc=2)
  beta.hat3<-matrix(0,nr=n.sim,nc=3)
  beta.hat4<-matrix(0,nr=n.sim,nc=3)
  V1<-V2<-V3<-V4<-rep(0,n.sim)
  l.beta1<-l.beta2<-l.beta3<-l.beta4<-rep(0,n.sim)
  for (k in 1:n.sim){ 
    data<-generate.ps(n.sample,beta0,beta1,.5,.5,0.5,0.5,tt)
    ################Fit Models############################# 
    data1<-gen.weights(data)
    model1 <- lm(log(y) ~ I+I1+z1+zz2-1, weights=w,data=data1)
    model2 <- lm(log(y) ~ I+I1-1, weights=w,data=data1)
    model3 <- lm(log(y) ~ I+I1+z1-1, weights=w,data=data1)
    model4 <- lm(log(y) ~ I+I1+zz2-1, weights=w,data=data1)
    
    beta.hat1[k,] <- c(model1$coeff[1],model1$coeff[2],model1$coeff[3],model1$coeff[4])
    beta.hat2[k,] <- c(model2$coeff[1],model2$coeff[2])
    beta.hat3[k,] <- c(model3$coeff[1],model3$coeff[2],model3$coeff[3])
    beta.hat4[k,] <- c(model4$coeff[1],model4$coeff[2],model4$coeff[3])
    l.beta1[k]<-beta.hat1[k,]%*%l1
    l.beta2[k]<-beta.hat2[k,]%*%l2
    l.beta3[k]<-beta.hat3[k,]%*%l3
    l.beta4[k]<-beta.hat4[k,]%*%l4
    V1[k]<-sqrt(t(l1)%*%vcov(model1)%*%l1)
    V2[k]<-sqrt(t(l2)%*%vcov(model2)%*%l2)
    V3[k]<-sqrt(t(l3)%*%vcov(model3)%*%l3)
    V4[k]<-sqrt(t(l4)%*%vcov(model4)%*%l4)
    reject[1]<-reject[1]+sign(abs(l.beta1[k]/V1[k])>1.96)
    reject[2]<-reject[2]+sign(abs(l.beta2[k]/V2[k])>1.96)
    reject[3]<-reject[3]+sign(abs(l.beta3[k]/V3[k])>1.96)
    reject[4]<-reject[4]+sign(abs(l.beta4[k]/V4[k])>1.96)
  }
  V<-c(mean(V1),mean(V2),mean(V3),mean(V4))
  SD<-c(sd(l.beta1),sd(l.beta2),sd(l.beta3),sd(l.beta4))
  xx[[1]]<-colMeans(beta.hat1)
  xx[[2]]<-colMeans(beta.hat2)
  xx[[3]]<-colMeans(beta.hat3)
  xx[[4]]<-colMeans(beta.hat4)
  xx[[5]]<-reject/n.sim
  xx[[6]]<-V
  xx[[7]]<-SD
  return(xx)
}

power.spb<-function(n.sim,n.sample,tt,beta0,beta1){
  print(paste("running spb",n.sim,n.sample,tt,beta0,beta1, sep=" "))
  #xx.names <- c("beta.hat1","beta.hat2","beta.hat3","beta.hat4","type I Err","V","SD")
  xx<-vector("list", 7)
  #names(xx) <- xx.names
  reject<-rep(0,4)
  beta.hat1<-matrix(0,nr=n.sim,nc=4)
  beta.hat2<-matrix(0,nr=n.sim,nc=2)
  beta.hat3<-matrix(0,nr=n.sim,nc=3)
  beta.hat4<-matrix(0,nr=n.sim,nc=3)
  V1<-V2<-V3<-V4<-rep(0,n.sim)
  l.beta1<-l.beta2<-l.beta3<-l.beta4<-rep(0,n.sim)
  for (k in 1:n.sim){ 
    data<-generate.spb(n.sample,beta0,beta1,.5,.5,0.5,0.5,tt)
    ################Fit Models############################# 
    data1<-gen.weights(data)
    model1 <- lm(log(y) ~ I+I1+z1+zz2-1, weights=w,data=data1)
    model2 <- lm(log(y) ~ I+I1-1, weights=w,data=data1)
    model3 <- lm(log(y) ~ I+I1+z1-1, weights=w,data=data1)
    model4 <- lm(log(y) ~ I+I1+zz2-1, weights=w,data=data1)
    
    beta.hat1[k,] <- c(model1$coeff[1],model1$coeff[2],model1$coeff[3],model1$coeff[4])
    beta.hat2[k,] <- c(model2$coeff[1],model2$coeff[2])
    beta.hat3[k,] <- c(model3$coeff[1],model3$coeff[2],model3$coeff[3])
    beta.hat4[k,] <- c(model4$coeff[1],model4$coeff[2],model4$coeff[3])
    l.beta1[k]<-beta.hat1[k,]%*%l1
    l.beta2[k]<-beta.hat2[k,]%*%l2
    l.beta3[k]<-beta.hat3[k,]%*%l3
    l.beta4[k]<-beta.hat4[k,]%*%l4
    V1[k]<-sqrt(t(l1)%*%vcov(model1)%*%l1)
    V2[k]<-sqrt(t(l2)%*%vcov(model2)%*%l2)
    V3[k]<-sqrt(t(l3)%*%vcov(model3)%*%l3)
    V4[k]<-sqrt(t(l4)%*%vcov(model4)%*%l4)
    reject[1]<-reject[1]+sign(abs(l.beta1[k]/V1[k])>1.96)
    reject[2]<-reject[2]+sign(abs(l.beta2[k]/V2[k])>1.96)
    reject[3]<-reject[3]+sign(abs(l.beta3[k]/V3[k])>1.96)
    reject[4]<-reject[4]+sign(abs(l.beta4[k]/V4[k])>1.96)
  }
  V<-c(mean(V1),mean(V2),mean(V3),mean(V4))
  SD<-c(sd(l.beta1),sd(l.beta2),sd(l.beta3),sd(l.beta4))
  xx[[1]]<-colMeans(beta.hat1)
  xx[[2]]<-colMeans(beta.hat2)
  xx[[3]]<-colMeans(beta.hat3)
  xx[[4]]<-colMeans(beta.hat4)
  xx[[5]]<-reject/n.sim
  xx[[6]]<-V
  xx[[7]]<-SD
  return(xx)
}

#########################z2^2 continuous CR n=500 censor=20%######################
power.cr(5000,500,10,.5,.5)
#########################z2^2 continuous PS n=500 censor=20%######################
power.ps(5000,500,10,.5,.5)
#########################z2^2 continuous SPB n=500 censor=20%######################
power.spb(5000,500,10,.5,.5)
#########################zz2^2 continuous CR n=250 censor=20%######################
power.cr(5000,250,10,.5,.5)
#########################zz2^2 continuous PS n=250 censor=20%######################
power.ps(5000,250,10,.5,.5)
#########################z2^2 continuous SPB n=250 censor=20%######################
power.spb(5000,250,10,.5,.5)




