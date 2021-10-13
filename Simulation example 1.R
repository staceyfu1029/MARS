## load R packages

library(simsurv)
library(survival)
library(purrr)
library(foreach)
library(doParallel)
library(parallel)
library(R2jags)


## Get job array id and setting from Batch
args <- commandArgs(trailingOnly = TRUE)

print(paste0("i = ", as.numeric(args[1])))
print(paste0("outFile = ", args[2]))

i <- as.numeric(args[1])
seed <- (i-1)*1000+1

## Simulation setting
n.Study <- 13
cen.rate <- 0.01
J <- 10
a <- 0:10*12 # 10 years of 12 months
n.sim <- 24
seed_list <- seq(seed, seed+n.sim, by = 1)


## Functions
simulated <- function(n.Study,seed){
  set.seed(seed)
  while(1){
    # Data generating -------------------------------
    alpha <- rnorm(n.Study,0,0.01)
    beta <- rnorm(n.Study,0,0.1)
    N <- sample(80:150,size=n.Study,replace = T)  # Sample size for each individual study
    Dur <- sample(16:24*5,size=n.Study,replace = T)  # Duration for each individual study
    
    data <- data.frame(x=numeric(),eventtime=numeric(),status=logical(),study=integer())
    
    for(i in 1:n.Study){
      p <- runif(1,min=0.45,max=0.55)  # Proportion of MRD negative in each study
      cov <- data.frame(id=1:N[i],x=as.numeric(rbernoulli(n=N[i],p=p))) 
      dt <- simsurv(lambdas=0.1,gammas=0.5,betas=c(x=-0.6+beta[i]),x=cov,maxt = Dur[i])
      dt <- merge(cov,dt,by="id")[,-1]
      cen <- rexp(N[i],cen.rate)
      dt$status <- (dt$eventtime<cen) & dt$status
      dt$eventtime <- round(pmin(dt$eventtime,cen),digits = 1)
      
      study <- rep(i,N[i])
      temp <- cbind(dt,study)
      data <- rbind(data,temp)
    }
    
    id <- 1:n.Study
    id.t1 <- sort(sample(1:n.Study,6))
    id.t2 <- sort(c(sample(id.t1,3),sample(id[!(id %in% id.t1)],3)))
    id.t3 <- sort(id[(!(id %in% id.t1)) & (!(id %in% id.t2))])
    
    T1 <- data[data$study %in% id.t1,]
    
    T2 <- data.frame(study=integer(),duration=integer(),LOG.HR=double(),LOG.HR.SE=double())
    for(i in id.t2)
    {
      d <- data[data$study==i,]
      fit <- coxph(Surv(eventtime,status)~x,data = d)
      temp <- data.frame(study=i,duration=Dur[i],LOG.HR=fit$coefficients,LOG.HR.SE=sqrt(fit$var[1,1]))
      T2 <- rbind(T2,temp)
    }
    
    T3 <- data.frame(study=integer(),time=integer(),SR=double(),SE=double(),group=integer())
    for(i in id.t3)
    {
      d <- data[data$study==i,]
      fit <- survfit(Surv(eventtime,status)~x,data = d)
      time <- sample(6:10*5,1)
      pred <- summary(fit,times = time)
      temp <- data.frame(study=rep(i,2),time=rep(time,2),SR=pred$surv,SE=pred$std.err,group=c(0,1))
      T3 <- rbind(T3,temp)
    }
    if(sum(T3$SR==0 | T3$SR==1)==0){
      break
    }
  }
  return(list(T1,T2,T3,data))
}

get_stats <- function(model){
  
  summ <- model$BUGSoutput$summary
  
  ## HR
  HR.hat <- summ[grep("beta.0",rownames(summ)),1]
  HR.se <- summ[grep("beta.0",rownames(summ)),2]
  HR.lower <- summ[grep("beta.0",rownames(summ)),3]
  HR.upper <- summ[grep("beta.0",rownames(summ)),7]                        
  
  b1 <- summ[grep("b1",rownames(summ)),1]
  w <- 0
  for(i in 1:J){
    w[i] <- (-0.6-b1[i])^2
  }
  HR.mse <- sum(w)/J
  
  #c.hr <- c("HR.hat","HR.se","HR.lower","HR.upper","HR.mse")
  hr <- c(HR.hat,HR.se,HR.lower,HR.upper,HR.mse)
  
  ## piecewise HR
  b <- c(summ[grep("b1",rownames(summ)),c(1,3,7,2)])
  
  ## SR
  sr0 <- exp(-c(summ[grep("L.0",rownames(summ)),c(1,7,3)])*12)
  sr1 <- exp(-c(summ[grep("L.1",rownames(summ)),c(1,7,3)])*12)
  
  
  dat <- model$BUGSoutput$sims.matrix
  l.0 <- dat[,grep("lambda.adj",colnames(dat))]
  l.1 <- dat[,grep("lambda.1",colnames(dat))]
  
  ## MS
  ms0 <- ms1 <- 0
  for(i in 1:dim(l.0)[1]){
    if(sum(cumsum(l.0[i,]*12)<log(2))==0){
      ms0[i] <- NA
    }else{
      index0 <- which.min(cumsum(l.0[i,]*12)<log(2))
      ms0[i] <- (log(2)-sum(l.0[i,1:(index0-1)]*12))/l.0[i,index0]+(index0-1)*12
    }
    if(sum(cumsum(l.1[i,]*12)<log(2))==0){
      ms1[i] <- NA
    }else{
      index1 <- which.min(cumsum(l.1[i,])*12<log(2))
      ms1[i] <- (log(2)-sum(l.1[i,1:(index1-1)]*12))/l.1[i,index1]+(index1-1)*12
    }
  }
  
  MS.0.hat <- mean(ms0,na.rm = T)
  MS.0.lower <- quantile(ms0,probs = 0.025,names = FALSE,na.rm = T)
  MS.0.upper <- quantile(ms0,probs = 0.975,names = FALSE,na.rm = T)
  MS.1.hat<- mean(ms1,na.rm = T)
  MS.1.lower <- quantile(ms1,probs = 0.025,names = FALSE,na.rm = T)
  MS.1.upper <- quantile(ms1,probs = 0.975,names = FALSE,na.rm = T)
  MS.D.hat <- mean(ms1-ms0,na.rm = T)
  MS.D.lower <- quantile(ms1-ms0,probs = 0.025,names = FALSE,na.rm = T)
  MS.D.upper <- quantile(ms1-ms0,probs = 0.975,names = FALSE,na.rm = T)
  
  ms <- c(MS.0.hat,MS.0.lower,MS.0.upper,MS.1.hat,MS.1.lower,MS.1.upper,MS.D.hat,MS.D.lower,MS.D.upper)
  
  ## RMST
  rmst <- c(t(sapply(c(1,3,5,10),RMST,l.0=l.0,l.1=l.1)))
  
  ## Correlation
  rho.h <- c(summ[grep("rho.h",rownames(summ)),c(5,3,7)])
  rho.b <- c(summ[grep("rho.b",rownames(summ)),c(5,3,7)])
  sigma.h <- c(summ[grep("sigma.hzd",rownames(summ)),c(5,3,7)])
  sigma.b <- c(summ[grep("sigma.b",rownames(summ)),c(5,3,7)])
  
  return(c(hr,b,sr0,sr1,ms,rmst,rho.h,rho.b,sigma.h,sigma.b))
}

RMST <- function(time,l.0,l.1){
  # X=0
  rmst.0 <- rmst.1 <- 0
  for(i in 1:dim(l.0)[1]){
    for(t in 1:time){
      if(t==1){
        rmst.0[i] <- (1-exp(-l.0[i,1]*12))/l.0[i,1]
      }else{
        rmst.0[i] <- rmst.0[i]+(exp(-sum(l.0[i,1:(t-1)])*12)-exp(-sum(l.0[i,1:t])*12))/l.0[i,t]
      }
    }
  }
  # X=1
  for(i in 1:dim(l.1)[1]){
    for(t in 1:time){
      if(t==1){
        rmst.1[i] <- (1-exp(-l.1[i,1]*12))/l.1[i,1]
      }else{
        rmst.1[i] <- rmst.1[i]+(exp(-sum(l.1[i,1:(t-1)])*12)-exp(-sum(l.1[i,1:t])*12))/l.1[i,t]
      }
    }
  }
  
  RMST0.hat <- mean(rmst.0)
  RMST0.lower <- quantile(rmst.0,probs = 0.025,names = FALSE)
  RMST0.upper <- quantile(rmst.0,probs = 0.975,names = FALSE)
  RMST1.hat<- mean(rmst.1)
  RMST1.lower <- quantile(rmst.1,probs = 0.025,names = FALSE)
  RMST1.upper <- quantile(rmst.1,probs = 0.975,names = FALSE)
  RMSTD.hat <- mean(rmst.1-rmst.0)
  RMSTD.lower <- quantile(rmst.1-rmst.0,probs = 0.025,names = FALSE)
  RMSTD.upper <- quantile(rmst.1-rmst.0,probs = 0.975,names = FALSE)
  RMSTR.hat <- median(rmst.1/rmst.0)
  RMSTR.lower <- quantile(rmst.1/rmst.0,probs = 0.025,names = FALSE)
  RMSTR.upper <- quantile(rmst.1/rmst.0,probs = 0.975,names = FALSE)
  
  c(RMST0.hat,RMST0.lower,RMST0.upper,RMST1.hat,RMST1.lower,RMST1.upper,
    RMSTD.hat,RMSTD.lower,RMSTD.upper,RMSTR.hat,RMSTR.lower,RMSTR.upper)
}

## Model 
model.full <- "
model{
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
    m2[i] <- beta[study2[i]]+b[1]*dt2[i,1]+b[2]*dt2[i,2]+b[3]*dt2[i,3]+b[4]*dt2[i,4]+b[5]*dt2[i,5]+b[6]*dt2[i,6]+b[7]*dt2[i,7]+b[8]*dt2[i,8]+b[9]*dt2[i,9]+b[10]*dt2[i,10]
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
"

params <- c("L.0","L.1","beta.0","lambda.adj","beta.adj","lambda.1","b1","rho.b","rho.h","sigma.hzd","sigma.b")

model.ph <- "
model{
  # type 2: hazard ratios
  for (i in 1:N2){        #N1 is the total number of data points for type 1 studies 
    y2[i] ~ dnorm(m2[i],s2[i])
    m2[i] <- beta[study2[i]]+b
  }
  
  # priors for parameters
  b ~ dnorm(0,0.001)
  for(i in 1:n.Study){
    beta[i] ~ dnorm(m.beta,sigma.beta)
  }
  m.beta ~ dnorm(0,0.001)
  sigma.beta ~ dgamma(0.01,0.01)
  
  for (i in 1:n.Study){
    beta.adj[i] <- beta[i] - mean(beta) 
  }
  
  b.adj <- mean(beta) + b
}
"

params2 <- c("b.adj","beta.adj")

model.type1 <- "
model{
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
"

params3 <- c("b.adj","beta.adj")


## Run simulation
start <- Sys.time()
registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))

rec.out <- foreach(count=1:n.sim, .combine = rbind, 
                   .packages = c("coxed","purrr","R2jags","simsurv","survival")) %dopar% {
                     seed.i <- seed_list[count]
                     print(paste("-------------------", count, "of ", n.sim, "simulated dataset with seed =", seed.i))
                     set.seed(seed.i)
                     
                     # Data generating -------------------------------
                     dt <- simulated(n.Study = n.Study,seed = seed.i)
                     T1 <- dt[[1]]
                     T2 <- dt[[2]]
                     T3 <- dt[[3]]
                     data <- dt[[4]]
                     
                     # IPD MA
                     # Type 1
                     N1 <- dim(data)[1]
                     study1 <- data$study
                     x1 <- data$x
                     censor <- as.numeric(data$status) # censor: 1=event,0=censored
                     t <- data$eventtime
                     d <- dt1 <- matrix(0,nrow = N1,ncol = J)
                     for(i in 1:N1){
                       for(j in 1:J){
                         d[i,j] <- censor[i]*(t[i]-a[j]>0)*(a[j+1]-t[i]>=0)
                         dt1[i,j] <- (min(t[i],a[j+1])-a[j])*(t[i]-a[j]>0) 
                       }
                     }
                     
                     data.list.t1 <- list(n.Study=n.Study,J=J,
                                          N1=N1,study1=study1,x1=x1,d=d,dt1=dt1)
                     
                     model.t1 <- jags(data = data.list.t1, 
                                      parameters.to.save = params, 
                                      model.file = textConnection(model.type1), 
                                      n.chains = 3, 
                                      n.burnin = 20000,
                                      n.iter = 30000, 
                                      n.thin = 1,
                                      progress.bar = "none")
                     
                     temp.t1 <- get_stats(model.t1)
                     
                     # AD analysis --------------------------------------
                     id <- sort(unique(T2$study))
                     N2 <- dim(T2)[1]
                     study2 <- sapply(T2$study,function(x) which(x==id))
                     y2 <- T2$LOG.HR
                     s2 <- 1/T2$LOG.HR.SE^2
                     
                     n.Study <- length(unique(study2))
                     
                     data.list2 <- list(n.Study=n.Study,
                                        N2=N2,study2=study2,y2=y2,s2=s2)
                     
                     model2 <- jags(data = data.list2, 
                                    parameters.to.save = params2, 
                                    model.file = textConnection(model.ph), 
                                    n.chains = 3, 
                                    n.burnin = 5000,
                                    n.iter = 10000, 
                                    n.thin = 1,
                                    progress.bar = "none")
                     
                     summ <- model2$BUGSoutput$summary
                     HR.ph.hat <- summ[1,1]
                     HR.ph.se <- summ[1,2]
                     HR.ph.lower <- summ[1,3]
                     HR.ph.upper <- summ[1,7]
                     HR.ph.mse <- (-0.6-HR.ph.hat)^2
                     
                     hr.ph <- c(HR.ph.hat,HR.ph.se,HR.ph.lower,HR.ph.upper,HR.ph.mse)
                     
                     # Res MA-------------------
                     T2 <- T2[!(T2$study %in% T1$study),]
                     
                     # Type 1
                     N1 <- dim(T1)[1]
                     study1 <- T1$study
                     x1 <- T1$x
                     censor <- as.numeric(T1$status) # censor: 1=event,0=censored
                     t <- T1$eventtime
                     d <- dt1 <- matrix(0,nrow = N1,ncol = J)
                     for(i in 1:N1){
                       for(j in 1:J){
                         d[i,j] <- censor[i]*(t[i]-a[j]>0)*(a[j+1]-t[i]>=0)
                         dt1[i,j] <- (min(t[i],a[j+1])-a[j])*(t[i]-a[j]>0) 
                       }
                     }
                     
                     # Type 2
                     N2 <- dim(T2)[1]
                     study2 <- T2$study
                     y2 <- T2$LOG.HR
                     s2 <- 1/T2$LOG.HR.SE^2
                     t2 <- T2$duration
                     dt2 <- matrix(0,nrow = N2,ncol = J)
                     for(i in 1:N2){
                       for(j in 1:J){
                         dt2[i,j] <- (min(t2[i],a[j+1])-a[j])*(t2[i]-a[j]>0)
                       }
                       dt2[i,] <- dt2[i,]/sum(dt2[i,])
                     }
                     
                     logit <- function(x){log(x/(1-x))}
                     # Type 3
                     N3 <- dim(T3)[1]
                     study3 <- T3$study
                     x3 <- T3$group
                     y3 <- logit(T3$SR)
                     s3 <- (T3$SR*(1-T3$SR))^2/T3$SE^2
                     t3 <- T3$time
                     dt3 <- matrix(0,nrow = N3,ncol = J)
                     for(i in 1:N3){
                       for(j in 1:J){
                         dt3[i,j] <- (min(t3[i],a[j+1])-a[j])*(t3[i]-a[j]>0)
                       }
                     }
                     
                     n.Study <- length(unique(c(study1,study2,study3)))
                     
                     data.list <- list(n.Study=n.Study,J=J,
                                       N1=N1,study1=study1,x1=x1,d=d,dt1=dt1,
                                       N2=N2,study2=study2,y2=y2,s2=s2,dt2=dt2,
                                       N3=N3,study3=study3,x3=x3,y3=y3,s3=s3,dt3=dt3)
                     
                     model <- jags(data = data.list, 
                                   parameters.to.save = params, 
                                   model.file = textConnection(model.full), 
                                   n.chains = 3, 
                                   n.burnin = 20000,
                                   n.iter = 30000, 
                                   n.thin = 1,
                                   progress.bar = "none")
                     
                     temp <- get_stats(model)
                     
                     c(temp.t1,temp,hr.ph)
                     
                   } 
end <- Sys.time()

end-start

## output data.frame
##  metrics to record: HR,b,SR,MS,RMST
stat <- c("hat","lower","upper")
group <- c("0","1")

nlist <- function(metric,list=NULL,side=T){
  name <- vector()
  if(length(list)==0){
    if(side==T){
      for(g in group){
        for(s in stat){
          name <- c(name,paste(metric,g,s,sep = "."))
        }
      }
    }else{
      for(s in stat){
        name <- c(name,paste(metric,s,sep = "."))
      }
    }
  }else{
    if(side==T){
      for(g in group){
        for(s in stat){ # order of loop corresponds to record
          for(i in list){
            name <- c(name,paste(metric,i,g,s,sep = "."))
          }
        }
      }
    }else{
      for(s in stat){ # order of loop corresponds to record
        for(i in list){
          name <- c(name,paste(metric,i,s,sep = "."))
        }
      }
    }
  }
  name
}

## HR
c.hr <- c("HR.hat","HR.se","HR.lower","HR.upper","HR.mse")

## piecewise HR
c.b <- nlist(metric="b",list=1:J,side=F)
c.b <- c(c.b,paste("b",1:J,"sd",sep="."))

## SR
c.sr <- nlist(metric="SR",list=1:J,side=T)

## MS
c.ms <- nlist(metric="MS",side=T)
c.msd <- nlist(metric = "MSD",side=F)

## RMST
rmst.list <- c(1,3,5,10)
c.rmst <- nlist(metric="RMST",list=rmst.list,side = T)
c.rmstd <- nlist(metric="RMSTD",list=rmst.list,side=F)
c.rmstr <- nlist(metric="RMSTR",list=rmst.list,side=F)

## correlation
c.rho.h <- nlist(metric="rho.h",side = F)
c.rho.b <- nlist(metric="rho.b",side = F)
c.sigma.h <- nlist(metric="sigma.h",side = F)
c.sigma.b <- nlist(metric="sigma.b",side = F)

## 3 methods
c.full <- c(c.hr,c.b,c.sr,c.ms,c.msd,c.rmst,c.rmstd,c.rmstr,c.rho.h,c.rho.b,c.sigma.h,c.sigma.b)
c.t1 <- unname(sapply(c.full, function(x) paste0(x,".t1")))
c.ph <- c("HR.ph.hat","HR.ph.se","HR.ph.lower","HR.ph.upper","HR.ph.mse")

colnames(rec.out) <- c(c.t1,c.full,c.ph)

save.image(args[2])





