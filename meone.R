library(R2OpenBUGS)
library(coda) 
library("BRugs")

load("E2.Rdata")
Y <- cbind(Y1,Y2,Y3)
n <- 365
p <- 6
dat <- list(Y1=Y1,Y2=Y2,Y3=Y3,Y=Y,p=6,n=25)

RE_model <- function(){
  
  Y1[1,1:p] ~ dmnorm(mu1[1:p],Omega1[1:p,1:p])
  Y2[1] ~ dnorm(mu1[1],tau)
  Y3[1] ~ dnorm(m3,tau)
  
  m3 <- mean(mu1[1:p])
  
  for(i in 2:n){
    Y1[i,1:p] ~ dmnorm(mpx[i-1,1:p],Omega2[1:p,1:p])
    Y2[i] ~ dnorm(mpxy[i-1],tau)
    Y3[i] ~ dnorm(mpxz[i-1],tau)
  }
  
  for(i in 1:n-1){
    for(j in 1:p){
      mpx[i,j] <- mu2[i,j] + rho*Y1[i,j]
    }
    mpxy[i] <- mu2[i,1] + rho*Y2[i]
    mpxz[i] <- mean(mu2[i,1:6]) + rho*Y3[i]
  }
  
  for(j in 1:p){
    mu1[j] ~ dnorm(0,0.01)
  }
  
  for(i in 1:n-1){
    for(j in 1:p){
      mu2[i,j] ~ dnorm(0,0.01)
    }
  }
  
  rho ~ dunif(0,1)
  tau ~ dgamma(0.01,0.01)
  
  k <- p+0.1
  Omega1[1:p,1:p] ~ dwish(R1[,],k)
  for(j1 in 1:p){for(j2 in 1:p){R1[j1,j2]<-0.1*equals(j1,j2)}} #R is diagonal
  
  Omega2[1:p,1:p] ~ dwish(R2[,],k)
  for(j1 in 1:p){for(j2 in 1:p){R2[j1,j2]<-0.1*equals(j1,j2)}} #R is diagonal
  
}



mlr_inits <- function() {
  list(mu1 = mean(dat$Y[1,1:7],na.rm=T),mu2=rowMeans(Y[2:n,1:7],na.rm=T),tau=0.01)
}


samps <- bugs(data = dat, 
              inits = mlr_inits,
              parameters.to.save = c("mu1","mpx","mpxy","mpxz"), 
              model.file = RE_model, 
              n.chains = 2, n.burnin=2e3, n.iter = 5e4,n.thin=10, DIC=F,debug = TRUE)

pdf(file="meone.pdf", width = 8, height = 12)
gelman.plot(samps)
dev.off()
