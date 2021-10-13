library(dclone)
library(rjags)
library(readr)

T1 <- read_csv("T1.csv")
T2 <- read_csv("T2.csv")
T3 <- read_csv("T3.csv")

J <- 22  # number of time segments
a <- 0
for(k in 1:(J+1)){
  a[k] <- 6*(k-1)
}

#--------------------------------------------------
#      Outcome: DFS
#--------------------------------------------------


T1.os <- T1[T1$event!="OS",]
T2.os <- T2[T2$Event!="OS",]
T3.os <- T3[T3$Event!="OS",]

# Only one type of data will be used
T2.os <- T2.os[!(T2.os$PMID %in% unique(T1.os$pmid)),]
T3.os <- T3.os[!(T3.os$PMID %in% unique(T1.os$pmid)),]
T3.os <- T3.os[!(T3.os$PMID %in% unique(T2.os$PMID)),]

idlist <- sort(unique(c(T1.os$pmid,T2.os$PMID,T3.os$PMID)))
T1.os$studyid <- sapply(T1.os$pmid, function(x) which(idlist==x))
T2.os$studyid <- sapply(T2.os$PMID, function(x) which(idlist==x))
T3.os$studyid <- sapply(T3.os$PMID, function(x) which(idlist==x))
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

# Type 2
N2 <- dim(T2.os)[1]
study2 <- T2.os$studyid
y2 <- ifelse(T2.os$ref=="P",T2.os$LOG.HR,-T2.os$LOG.HR)
s2 <- 1/T2.os$LOG.HR.SE^2
t2 <- T2.os$time
dt2 <- matrix(0,nrow = N2,ncol = J)
for(i in 1:N2){
  for(j in 1:J){
    dt2[i,j] <- (min(t2[i],a[j+1])-a[j])*(t2[i]-a[j]>0)
  }
  dt2[i,] <- dt2[i,]/sum(dt2[i,])
}

logit <- function(x){log(x/(1-x))}
# Type 3
N3 <- dim(T3.os)[1]
study3 <- T3.os$studyid
x3 <- T3.os$MRD_group
y3 <- logit(T3.os$`survival rate`)
s3 <- (T3.os$`survival rate`*(1-T3.os$`survival rate`))^2/T3.os$std.err^2
t3 <- T3.os$`time (in months)`
dt3 <- matrix(0,nrow = N3,ncol = J)
for(i in 1:N3){
  for(j in 1:J){
    dt3[i,j] <- (min(t3[i],a[j+1])-a[j])*(t3[i]-a[j]>0)
  }
}


data <- list(n.Study=n.Study,J=J,
             N1=N1,study1=study1,x1=x1,d=d,dt1=dt1,
             N2=N2,study2=study2,y2=y2,s2=s2,dt2=dt2,
             N3=N3,study3=study3,x3=x3,y3=y3,s3=s3,dt3=dt3)


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
  
  # type 2: hazard ratios
  for (i in 1:N2){        #N1 is the total number of data points for type 1 studies 
      y2[i] ~ dnorm(m2[i],s2[i])
      m2[i] <- beta[study2[i]]+b[1]*dt2[i,1]+b[2]*dt2[i,2]+b[3]*dt2[i,3]+b[4]*dt2[i,4]+b[5]*dt2[i,5]+b[6]*dt2[i,6]+b[7]*dt2[i,7]+b[8]*dt2[i,8]+b[9]*dt2[i,9]+b[10]*dt2[i,10]+b[11]*dt2[i,11]+b[12]*dt2[i,12]+b[13]*dt2[i,13]+b[14]*dt2[i,14]+b[15]*dt2[i,15]+b[16]*dt2[i,16]+b[17]*dt2[i,17]+b[18]*dt2[i,18]+b[19]*dt2[i,19]+b[20]*dt2[i,20]+b[21]*dt2[i,21]+b[22]*dt2[i,22]
  }
  
  # type 3: survival rates
  for(i in 1:N3){
    y3[i] ~ dnorm(m3[i],s3[i])
    m3[i] <- log(u3[i]/(1-u3[i]))
    u3[i] <- exp(-sum(cum.hzd[i,])*exp(alpha[study3[i]]+beta[study3[i]]*x3[i]))
    for(j in 1:J){
      cum.hzd[i,j] <- lambda[j]*dt3[i,j]*exp(b[j]*x3[i])
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
                        n.update = 100000,
                        n.iter = 70000, 
                        thin = 1)

stopCluster(cl)
end <- Sys.time()

end-Start

save.image("All type.RData")
