library(dclone)
library(rjags)
library(readr)

set.seed(1112)

T1 <- read_csv("T1.csv")

J <- 22  # number of time segments
a <- 0
for(k in 1:(J+1)){
  a[k] <- 6*(k-1)
}

#--------------------------------------------------
#      Outcome: DFS
#--------------------------------------------------

T1.os <- T1[T1$event!="OS",]
#T2.os <- T2[T2$Event=="OS",]
#T3.os <- T3[T3$Event=="OS",]

# Only one type of data will be used
#T2.os <- T2.os[!(T2.os$PMID %in% unique(T1.os$pmid)),]
#T3.os <- T3.os[!(T3.os$PMID %in% unique(T1.os$pmid)),]
#T3.os <- T3.os[!(T3.os$PMID %in% unique(T2.os$PMID)),]

idlist <- sort(unique(c(T1.os$pmid)))
T1.os$studyid <- sapply(T1.os$pmid, function(x) which(idlist==x))
#T2.os$studyid <- sapply(T2.os$PMID, function(x) which(idlist==x))
#T3.os$studyid <- sapply(T3.os$PMID, function(x) which(idlist==x))
n.Study <- length(idlist)


# Type 1
N1 <- dim(T1.os)[1]
study1 <- T1.os$studyid
x1 <- T1.os$group
censor <- T1.os$censor
t <- T1.os$survtime
d <- dt1 <- matrix(0,nrow = N1,ncol = J)
for(i in 1:N1){
  for(j in 1:J){
    d[i,j] <- censor[i]*(t[i]-a[j]>0)*(a[j+1]-t[i]>=0);
    dt1[i,j] <- (min(t[i],a[j+1])-a[j])*(t[i]-a[j]>0); 
  }
}


data <- list(n.Study=n.Study,J=J,
             N1=N1,study1=study1,x1=x1,d=d,dt1=dt1)


# Model
model <- function(){
  # type 1: reconstructed survival data
  for(i in 1:N1){  #N1 is the total number of patients with reconstructed survival data  
    for(j in 1:J){
      h[i,j] <- lambda[j]*exp(alpha[study1[i]]+(beta[study1[i]]+b[j])*x1[i]);  # study[] is a mapping from 1:N to 1:n.Study; n2 is the total number of type 2 studies
      
      d[i,j] ~ dpois(mu[i,j]);
      mu[i,j] <- dt1[i,j]*h[i,j];
    }
  }
  
  # priors for parameters
  hzd[1:J] ~ dmnorm.vcov(mu.h[1:J],Sigma.h[1:J,1:J])
  b[1:J] ~ dmnorm.vcov(mu.b[1:J],Sigma.b[1:J,1:J])
  for(i in 1:J){
    mu.h[i] <- m.hzd
    mu.b[i] <- m.b
    for(j in 1:J){
      covar.h[i,j] <- pow(rho.h,abs(j-i))   # autoregressive covariance structure
      covar.b[i,j] <- pow(rho.b,abs(j-i)) 
    }
  }
  Sigma.b[1:J,1:J] <- tau.b * covar.b[1:J,1:J]
  tau.b <- inverse(sigma.b)
  Sigma.h[1:J,1:J] <- tau.h * covar.h[1:J,1:J]
  tau.h <- inverse(sigma.hzd)
  
  for(i in 1:J){
    lambda[i] <- exp(hzd[i])
  }
 
  for(i in 1:n.Study){
    alpha[i] ~ dnorm(m.alpha,sigma.alpha)
    beta[i] ~ dnorm(m.beta,sigma.beta)
  }
  m.hzd ~ dnorm(0,0.001)
  m.alpha ~ dnorm(0,0.001)
  m.beta ~ dnorm(0,0.001)
  m.b ~ dnorm(0,0.001)
  sigma.hzd ~ dgamma(0.01,0.01)
  sigma.alpha ~ dgamma(0.01,0.01)
  sigma.beta ~ dgamma(0.01,0.01)
  sigma.b ~ dgamma(0.01,0.01)
  rho.h ~ dunif(-0.99, 0.99)
  rho.b ~ dunif(-0.99, 0.99)
  
  for (i in 1:n.Study){
    alpha.adj[i] <- alpha[i] - mean(alpha)
    beta.adj[i] <- beta[i] - mean(beta) 
  }
  
  beta.0 <- mean(beta) + mean(b)

  for (j in 1:J){
    lambda.adj[j] <- lambda[j]*exp(mean(alpha))
    b.adj[j] <- b[j] - mean(b) 
    b1[j] <- b.adj[j]+beta.0
    lambda.1[j] <- lambda.adj[j]*exp(b1[j])
  }
  for (j in 1:J){
    L.0[j] <- sum(lambda.adj[1:j])
    L.1[j] <- sum(lambda.1[1:j])
  }
}


params <- c("lambda.adj","lambda.1","L.0","L.1","alpha.adj","b.adj","b1","beta.adj","beta.0","rho.b","rho.h","sigma.hzd","sigma.b")

Start <- Sys.time()
cl <- makePSOCKcluster(3)
tmp <- clusterEvalQ(cl, library(dclone))
parLoadModule(cl, 'glm')
parLoadModule(cl, 'lecuyer')
parLoadModule(cl, 'dic')

model.t1 <- jags.parfit(cl = cl, data = data, params = params, model = model, 
                        n.chains = 3, 
                        n.update = 500,
                        n.iter = 300, 
                        thin = 1)

stopCluster(cl)
end <- Sys.time()

end-Start

save.image("T1.RData")
