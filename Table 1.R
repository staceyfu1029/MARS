rm(list = ls())
## Simulation result table 1

setwd("path")

## table
t <- data.frame(matrix(NA,nrow = 8,ncol = 6))
t[1,3] <- "log HR"
t[2,3:6] <- c("Average(SD)","MSE","CP","Length of CrIs")
t[c(3,6),1] <- c("Case 1","Case 2")
t[3:8,2] <- rep(c("MARS","IPD","AD"),2)

## Scenario 1
#load("Routput/sim1.RData")
load("Routput/sim1.RData")

out <- data.frame(out)
x <- cbind(out$HR.hat,out$HR.hat.t1,out$HR.ph.hat)
for(i in 1:3){
  temp <- x[,i]
  d <- paste0(format(round(mean(temp,na.rm = T),digits = 3),nsmall = 2)," (",
              format(round(sd(temp,na.rm = T),digits = 2),nsmall = 2),")")
  d[2] <- format(round(mean((temp-1)^2,na.rm = T),digits=2),nsmall = 2)
  if(i==1){
    t.sub <- d
  }else{
    t.sub <- rbind(t.sub,d)
  }
} 
t[3:5,3:4] <- t.sub
tv <- -0.6
t[3,5] <- round(mean(out$HR.lower<tv & out$HR.upper>tv,na.rm = T),digits = 2)
t[4,5] <- round(mean(out$HR.ph.lower<tv & out$HR.ph.upper>tv,na.rm = T),digits = 2)
t[5,5] <- round(mean(out$HR.lower.t1<tv & out$HR.upper.t1>tv,na.rm = T),digits = 2)

t[3,6] <- round(mean(out$HR.upper-out$HR.lower),digits = 2)
t[4,6] <- round(mean(out$HR.ph.upper-out$HR.ph.lower),digits = 2)
t[5,6] <- round(mean(out$HR.upper.t1-out$HR.lower.t1),digits = 2)

## Scenario 2
#load("Routput/sim2.RData")
load("Routput/sim2.RData")

out <- data.frame(out)
x <- cbind(out$HR.hat,out$HR.hat.t1,out$HR.ph.hat)
mse <- cbind(out$HR.mse,out$HR.mse.t1,out$HR.ph.mse)
for(i in 1:3){
  temp <- x[,i]
  d <- paste0(format(round(mean(temp,na.rm = T),digits = 3),nsmall = 2)," (",
              format(round(sd(temp,na.rm = T),digits = 2),nsmall = 2),")")
  d[2] <- format(round(mean(mse[,i],na.rm = T),digits=2),nsmall = 2)
  if(i==1){
    t.sub <- d
  }else{
    t.sub <- rbind(t.sub,d)
  }
} 
t[6:8,3:4] <- t.sub

writexl::write_xlsx(t,"Table/table 1.xlsx",col_names = F)
