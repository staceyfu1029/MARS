rm(list = ls())
setwd("Path")

load("Routput/sim2.RData")

library(ggpubr)

out <- data.frame(out)

## Panel A: HR
x <- as.matrix(out[,180:189])
x.i <- as.matrix(out[,6:15])
x.ph <- as.matrix(out$HR.ph.hat)

dt <- data.frame(logHR= c(x),
                 time=factor(rep(1:dim(x)[2],each=dim(x)[1])),
                 method=rep("MARS",each=dim(x)[1]*dim(x)[2]))
dt.i <- data.frame(logHR= c(x.i),
                   time=factor(rep(1:dim(x.i)[2],each=dim(x.i)[1])),
                   method=rep("IPD",each=dim(x.i)[1]*dim(x.i)[2]))
dt.ph <- data.frame(logHR = rep(x.ph,J),
                    time=factor(rep(1:dim(x)[2],each=dim(x)[1])),
                    method=rep("AD",each=dim(x)[1]*dim(x)[2]))
dt <- rbind(dt,dt.i,dt.ph)
dt$method <- factor(dt$method,levels=c("MARS","IPD","AD"))

f <- function(x){
  1-0.2^(x/3-1)
}

avg <- function(x){
  integrate(f,lower = x-1,upper = x)$value
}

y <- sapply(1:10, avg)
tv <- data.frame(x=1:10,y=y)
p2 <- ggboxplot(dt, y = "logHR",
                x = "time",
                color = "method",
                ylab = "HR",
                xlab = "Year",
                palette = c("#F8766DFF", "#00BFC4FF", "#00BA38FF")
)+
  scale_y_continuous(breaks=log(c(0.1,0.2,0.4,1,2)),
                     labels=c(0.1,0.2,0.4,1,2))

pA <- p2 + geom_segment(data = tv,aes(x=as.numeric(x)-0.5,xend=as.numeric(x)+0.5,
                                      y=y,yend=y),
                        linetype = "longdash",
                        size = 0.6)+
  theme_pubr(base_size=10)

## Panel B: RMSTD
hzd0 <- function(t){
  0.01
}
hzd1 <- function(t){
  0.01*exp(1-0.2^(t/30-1))
}
surv0 <- function(t){
  exp(-0.01*t)
}
surv1 <- function(t){
  sapply(t,function(x) {exp(-integrate(hzd1,0,x)$value)})
}
cumhzd <- function(x)
{
  integrate(hzd1,0,x)$value-log(2)
}

for(time in c(1,3,5,10)){
  eval(parse(text = paste0("x0 <- out$RMSTD.",time,".hat")))
  eval(parse(text = paste0("x0.i <- out$RMSTD.",time,".hat.t1")))
  
  dt <- data.frame(rmstd= c(x0),
                   time=time,
                   method=rep("MARS",length(x0)))
  dt.i <- data.frame(rmstd= c(x0.i),
                     time=time,
                     method=rep("IPD",length(x0.i)))
  dt <- rbind(dt,dt.i)
  dt$method <- factor(dt$method,levels=c("MARS","IPD"))
  
  ## Calculate true values of RMST
  y <- (integrate(surv1,0,time*12)$value-integrate(surv0,0,time*12)$value)/12
  temp <- data.frame(time=time,y=y)
  
  if(time==1){
    dt.overall <- dt
    tv <- temp
  }else{
    dt.overall <- rbind(dt.overall,dt)
    tv <- rbind(tv,temp)
  }
}

dt.overall$rmstd <- dt.overall$rmstd/12
tv$time <- factor(tv$time,levels = c(1,3,5,10))
p <- ggboxplot(dt.overall, y = "rmstd",
               x = "time",
               color = "method",
               ylab = "RMST difference (in years)",
               xlab = "Year",
               palette = c("#F8766DFF", "#00BFC4FF", "#00BA38FF")
)


pB <- p+geom_segment(data = tv,aes(x=as.numeric(time)-0.5,xend=as.numeric(time)+0.5,
                                   y=y,yend=y),
                     linetype = "longdash",
                     size = 0.6)+
  theme_pubr(base_size=10)

## Panel C,D: SR
x0 <- as.matrix(out[,220:229])
x0.i <- as.matrix(out[,46:55])

dt <- data.frame(sur= c(x0),
                 time=rep(1:10,each=dim(x0)[1]),
                 method=rep("MARS",each=dim(x0)[1]*dim(x0)[2]))
dt.i <- data.frame(sur= c(x0.i),
                   time=rep(1:10,each=dim(x0)[1]),
                   method=rep("IPD",each=dim(x0.i)[1]*dim(x0.i)[2]))
dt <- rbind(dt,dt.i)
dt$method <- factor(dt$method,levels=c("MARS","IPD","AD"))


## Calculate true values of SR
hzd0 <- function(t){
  0.01
}
surv0 <- function(t){
  exp(-0.01*t)
}

y <- surv0(1:10*12)
tv0 <- data.frame(sur=y,time=1:10)

p <- ggboxplot(dt, x="time", y = "sur",
               color = "method",
               ylab = "Survival rates",
               xlab = "Year",
               palette = c("#F8766DFF", "#00BFC4FF", "#00BA38FF")
)
p <- ggpar(p,ylim = c(0,1))

pC <- p+geom_segment(data = tv0,aes(x=as.numeric(time)-0.5,xend=as.numeric(time)+0.5,
                                    y=sur,yend=sur),
                     linetype = "longdash",
                     size = 0.6)+
  theme_pubr(base_size=10)

x1 <- as.matrix(out[,250:259])
x1.i <- as.matrix(out[,76:85])

dt <- data.frame(sur= c(x1),
                 time=rep(1:10,each=dim(x1)[1]),
                 method=rep("MARS",each=dim(x1)[1]*dim(x1)[2]))
dt.i <- data.frame(sur= c(x1.i),
                   time=rep(1:10,each=dim(x1)[1]),
                   method=rep("IPD",each=dim(x1.i)[1]*dim(x1.i)[2]))
dt <- rbind(dt,dt.i)
dt$method <- factor(dt$method,levels=c("MARS","IPD"))


## Calculate true values of SR
hzd1 <- function(t){
  0.01*exp(1-0.2^(t/30-1))
}
surv1 <- function(t){
  sapply(t,function(x) {exp(-integrate(hzd1,0,x)$value)})
}

y <- surv1(1:10*12)
tv <- data.frame(sur=y,time=1:10)

p <- ggboxplot(dt, x="time", y = "sur",
               color = "method",
               ylab = "",
               xlab = "Year",
               palette = c("#F8766DFF", "#00BFC4FF", "#00BA38FF"))
p <- ggpar(p,ylim = c(0,1))


pD <- p+geom_segment(data = tv,aes(x=as.numeric(time)-0.5,xend=as.numeric(time)+0.5,
                                   y=sur,yend=sur),
                     linetype = "longdash",
                     size = 0.6)+
  theme_pubr(base_size=10)

## MS
x0 <- out$MS.0.hat
x1 <- out$MS.1.hat
x0.i <- out$MS.0.hat.t1
x1.i <- out$MS.1.hat.t1

dt <- data.frame(y=c(x0,x1)/12,
                 time=factor(rep(1:2,each=length(x0))),
                 method=rep("MARS-MA",each=length(x0)))
dt.i <- data.frame(y=c(x0.i,x1.i)/12,
                   time=factor(rep(1:2,each=length(x0.i))),
                   method=rep("IPD",each=length(x0.i)))
dt <- rbind(dt,dt.i)
dt$method <- factor(dt$method,levels=c("MARS-MA","IPD"))

y <- c(log(2)/0.01,uniroot(cumhzd,interval = c(1,100))$root)/12
tv <- data.frame(time=factor(1:2),y=y)

p <- ggboxplot(dt, x="time", y = "y",
               color = "method",
               ylab = "MS",
               xlab = "Group"
)
pE <- p+geom_segment(data=tv, aes(x=as.numeric(time)-0.5,xend=as.numeric(time)+0.5,
                            y=y,yend=y),
               linetype = "longdash",
               size = 0.6)+
  theme_pubr(base_size = 10)+
  scale_x_discrete(labels = list("1"="x=0","2"="x=1"))


pdf("Figures/Figure2.pdf",width = 8,height = 6,onefile = FALSE)
ggarrange(pA,ggplot() + theme_void(),pC,pB,ggplot() + theme_void(),pD,
          labels = c("A","","B","C","","D"),widths = c(1, 0.08, 1),
          nrow = 2, ncol = 3, common.legend = T,
          font.label=list(color="black",size=8))
dev.off()
