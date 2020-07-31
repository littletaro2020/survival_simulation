rm(list=ls())
library(survival)
library(KMsurv)
library(MASS)

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

reassign.ps<-function(data){
  n.sample<-nrow(data)
  z1<-data$z1
  z2<-data$z2
  n<-matrix(0,2,4)
  z<-cbind(z1,1-z1,z2,1-z2)
  I<-rep(0,n.sample)
  for (i in 1:n.sample){
    I[i]<-assign.ps(z[i,],n,0.75)
    if(I[i]) n[1,]<-n[1,]+z[i,]
    else n[2,]<-n[2,]+z[i,]
  }
  data$I<-I
  return(data)
}

reassign.spb<-function(data){
  n.sample<-nrow(data)
  z1<-data$z1
  z2<-data$z2
  data$I<-assign.spb(z1,z2)$I
  return(data)
}

########################################################

####################generate dataset####################
generate.ps<-function(n.sample,beta0,beta1,beta2,beta3,p1,p2,tt,rou){
  #err<-runif(n.sample,-0.5,0.5)
  err1<-rnorm(n.sample,0,0.25)
  err2 <-append(0,err1[1:n.sample-1])
  #  rou <- 0.1
  c<-runif(n.sample,0,tt)
  z1<-rbinom(n.sample,1,p1)
  z2<-rbinom(n.sample,1,p2)
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
  log.t<-beta0*(1-I)+beta1*I+beta2*z1+beta3*z2+err1+err2*rou
  t<-exp(log.t)
  delta<-sign(t<=c)
  y=pmin(t,c)
  data<-data.frame(I,z1,z2,t,c,y,delta)
  return(data)
}

generate.spb<-function(n.sample,beta0,beta1,beta2,beta3,p1,p2,tt,rou){
  #err<-runif(n.sample,-0.5,0.5)
  err1<-rnorm(n.sample,0,0.25)
  err2 <-append(0,err1[1:n.sample-1])
  #  rou <- 0.1
  c<-runif(n.sample,0,tt)
  z1<-rbinom(n.sample,1,p1)
  z2<-rbinom(n.sample,1,p2)
  data.temp<-assign.spb(z1,z2)
  z1<-data.temp$z1-0.5
  z2<-data.temp$z2-0.5
  I<-data.temp$I
  log.t<-beta0*(1-I)+beta1*I+beta2*z1+beta3*z2+err1+err2*rou
  t<-exp(log.t)
  delta<-sign(t<=c)
  y=pmin(t,c)
  data<-data.frame(I,z1,z2,t,c,y,delta)
  return(data)
}

###########################################################

fit.model <- function(data){
  retults<-vector("list", 3)
  data_w<-gen.weights(data)
  model2 <- lm(log(y) ~ I+I1-1, weights=w,data=data_w)
  model3 <- lm(log(y) ~ I+I1+z1-1, weights=w,data=data_w)
  model4 <- lm(log(y) ~ I+I1+z2-1, weights=w,data=data_w)
  beta.hat2<- c(model2$coeff[1],model2$coeff[2])
  beta.hat3<- c(model3$coeff[1],model3$coeff[2],model3$coeff[3])
  beta.hat4<- c(model4$coeff[1],model4$coeff[2],model4$coeff[3])
  retults[[1]]<-beta.hat2%*%l2
  retults[[2]]<-beta.hat3%*%l3
  retults[[3]]<-beta.hat4%*%l4
  return(retults)
}


args <- commandArgs(TRUE)
randomization_method <- args[1]
if (randomization_method == 'ps') {
  print("using ps method")
  generate <- generate.ps
  reassign <- reassign.ps
} else if (randomization_method == 'spb') {
  print("using spb method")
  generate <- generate.spb
  reassign <- reassign.spb
} else{
  print("method is not supported")
}
n.sim<-as.numeric(args[2])
n.sample<- as.numeric(args[3])
n.boot<-as.numeric(args[4])

tt<-10
beta0<-beta1<-0.5
reject<-0
##n.sim<-1000
##n.sample<- 500
##n.boot<-300

l.beta.hat2<-rep(0,n.sim)
l.beta.hat3<-rep(0,n.sim)
l.beta.hat4<-rep(0,n.sim)

l.beta2<-rep(0,n.boot)
l.beta3<-rep(0,n.boot)
l.beta4<-rep(0,n.boot)

reject<-rep(0,3)

V2<-V3<-V4<-rep(0,n.sim)
for (k in 1:n.sim){ 
  print(k)
  print(reject)
  data0<-generate(n.sample,beta0,beta1,.5,.5,0.5,0.5,tt,0.1)
  l.beta.hat2[k]<-fit.model(data0)[[1]][1,1]
  l.beta.hat3[k]<-fit.model(data0)[[2]][1,1]
  l.beta.hat4[k]<-fit.model(data0)[[3]][1,1]
  for (j in 1:n.boot){
    data<-reassign(data0[sample(nrow(data0), n.sample,replace=TRUE), ])
    ################Fit Models############################# 
    l.beta2[j]<-fit.model(data)[[1]][1,1]
    l.beta3[j]<-fit.model(data)[[2]][1,1]
    l.beta4[j]<-fit.model(data)[[3]][1,1]
  }
  V2[k]<-sqrt(var(l.beta2))
  V3[k]<-sqrt(var(l.beta3))
  V4[k]<-sqrt(var(l.beta4))
  reject[1]<-reject[1]+sign(abs(l.beta.hat2[k]/V2[k])>1.96)
  reject[2]<-reject[2]+sign(abs(l.beta.hat3[k]/V3[k])>1.96)
  reject[3]<-reject[3]+sign(abs(l.beta.hat4[k]/V4[k])>1.96)
}

type1error <- reject/n.sim
print(type1error)
#write.table(data.frame(l.beta.hat1,l.beta.hat2,V1,V2,beta.hat1,beta.hat2), "C:/Users/waykinglu/Desktop/Research/Zhu/survival/Simulations/summary.txt", sep="\t")
#save(reject,file="/Users/l0w013f/Desktop/SPH/New_simulations/summary.rda")
