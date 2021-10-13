rm(list=ls())

library(ggplot2)
library(ggpubr)

#setwd("path")
# load dataset
#load("All type.RData")


# plotting
record <- summ[[2]]

Lam.t.0 <- approx(0:J, c(0,(record[,3][grep("L.0",names(record[,1]))]*6)),n=J*6+1)  #Lambda(t) for MRD=0
S.t.0 <- exp(-Lam.t.0[[2]])   #s(t) = exp(-Lambda(t))
title(paste0("3 Types of Data"),cex.main=1)
# 95%
lcl.0 <- exp(-approx(0:J, c(0,(record[,1][grep("L.0",names(record[,1]))]*6)),n=J*6+1)[[2]])
ucl.0 <- exp(-approx(0:J, c(0,(record[,5][grep("L.0",names(record[,1]))]*6)),n=J*6+1)[[2]])
S.t.0.1 <- S.t.0

# MRD group
Lam.t.1 <- approx(0:J, c(0,(record[,3][grep("L.1",names(record[,1]))]*6)),n=J*6+1)  #Lambda(t) for MRD=0
S.t.1 <- exp(-Lam.t.1[[2]])   #s(t) = exp(-Lambda(t))
# 95%
lcl.1 <- exp(-approx(0:J, c(0,(record[,1][grep("L.1",names(record[,1]))]*6)),n=J*6+1)[[2]])
ucl.1 <- exp(-approx(0:J, c(0,(record[,5][grep("L.1",names(record[,1]))]*6)),n=J*6+1)[[2]])

S.t.1.1 <- S.t.1

dt <- data.frame(x=Lam.t.0[[1]]/2,s0=S.t.0,s0.lower=lcl.0,s0.upper=ucl.0)
dt$s1 <- S.t.1
dt$s1.lower <- lcl.1
dt$s1.upper <- ucl.1

## ggplot2
p <- ggplot(dt, aes(x=x)) + 
  ylim(0,1)+
  #xlim(0,12)+
  geom_ribbon(aes(ymin=s0.lower,ymax=s0.upper), fill="darksalmon",alpha=0.5) +
  geom_ribbon(aes(ymin=s1.lower,ymax=s1.upper), fill="deepskyblue3",alpha=0.5) +
  geom_line(aes(y=s0,col="MRD+"),size=1,show.legend = F) +
  geom_line(aes(y=s1,col="MRD-"),size=1,show.legend = F) +
  labs(x = "Year",
       y = "survival rate")+
  scale_color_manual(name = "group",
                     values = c("MRD+" = "coral1", "MRD-" = "blue2"))+ 
  theme_pubr(base_size = 16)+
  annotate("text", x = 7, y = 0.8, label = "MRD negative")+
  annotate("text", x = 5, y = 0.1, label = "MRD positive")+
  grids("y",linetype = "dashed")
p

# Figure 1b
# load dataset
load("T1.RData")


# plotting
record <- summ[[2]]

Lam.t.0 <- approx(0:J, c(0,(record[,3][grep("L.0",names(record[,1]))]*6)),n=J*6+1)  #Lambda(t) for MRD=0
S.t.0 <- exp(-Lam.t.0[[2]])   #s(t) = exp(-Lambda(t))
lcl.0 <- exp(-approx(0:J, c(0,(record[,1][grep("L.0",names(record[,1]))]*6)),n=J*6+1)[[2]])
ucl.0 <- exp(-approx(0:J, c(0,(record[,5][grep("L.0",names(record[,1]))]*6)),n=J*6+1)[[2]])
S.t.0.1 <- S.t.0

# MRD group
Lam.t.1 <- approx(0:J, c(0,(record[,3][grep("L.1",names(record[,1]))]*6)),n=J*6+1)  #Lambda(t) for MRD=0
S.t.1 <- exp(-Lam.t.1[[2]])   #s(t) = exp(-Lambda(t))
# 95%
lcl.1 <- exp(-approx(0:J, c(0,(record[,1][grep("L.1",names(record[,1]))]*6)),n=J*6+1)[[2]])
ucl.1 <- exp(-approx(0:J, c(0,(record[,5][grep("L.1",names(record[,1]))]*6)),n=J*6+1)[[2]]) 
S.t.1.1 <- S.t.1

dt <- data.frame(x=Lam.t.0[[1]]/2,s0=S.t.0,s0.lower=lcl.0,s0.upper=ucl.0)
dt$s1 <- S.t.1
dt$s1.lower <- lcl.1
dt$s1.upper <- ucl.1

## ggplot2
p1 <- ggplot(dt, aes(x=x)) + 
  ylim(0,1)+
  #xlim(0,12)+
  geom_ribbon(aes(ymin=s0.lower,ymax=s0.upper), fill="darksalmon",alpha=0.5) +
  geom_ribbon(aes(ymin=s1.lower,ymax=s1.upper), fill="deepskyblue3",alpha=0.5) +
  geom_line(aes(y=s0,col="MRD+"),size=1,show.legend = F) +
  geom_line(aes(y=s1,col="MRD-"),size=1,show.legend = F) +
  labs(x = "Year",
       y = "survival rate")+
  scale_color_manual(name = "group",
                     values = c("MRD+" = "coral1", "MRD-" = "blue2"))+ 
  theme_pubr(base_size = 16)+
  annotate("text", x = 8, y = 0.8, label = "MRD negative")+
  annotate("text", x = 5, y = 0.1, label = "MRD positive")+
  grids("y",linetype = "dashed")
p1

##RMST
sample <- do.call(rbind,model.t1)

l.0 <- sample[,grep("lambda.adj",colnames(sample))]
l.1 <- sample[,grep("lambda.1",colnames(sample))]

for(time in c(1,3,5,10)){
  rmst.p <- rmst.n <- 0
  for(i in 1:dim(l.0)[1]){
    rmst.p[i] <- (1-exp(-l.0[i,1]*6))/l.0[i,1]
    for(j in 2:(2*time)){
      rmst.p[i] <- rmst.p[i]+(exp(-sum(l.0[i,1:(j-1)])*6)-exp(-sum(l.0[i,1:j])*6))/l.0[i,j]
    }
  }
  
  for(i in 1:dim(l.0)[1]){
    rmst.n[i] <- (1-exp(-l.1[i,1]*6))/l.1[i,1]
    for(j in 2:(2*time)){
      rmst.n[i] <- rmst.n[i]+(exp(-sum(l.1[i,1:(j-1)])*6)-exp(-sum(l.1[i,1:j])*6))/l.1[i,j]
    }
  }
  temp <- matrix(c(time,mean(rmst.n)/12,quantile(rmst.n,probs=c(0.025,0.975))/12,1,
                   time,mean(rmst.p)/12,quantile(rmst.p,probs=c(0.025,0.975))/12,0),
                 ncol=5,byrow = T)
  if(time==1){
    dt <- temp
  }else{
    dt <- rbind(dt,temp)
  }
}

dt <- data.frame(dt)
colnames(dt) <- c("time","rmst","lower","upper","group")
for(i in 1:dim(dt)[1]){
  txt <- paste(format(round(dt$rmst[i],2),nsmall = 2)," (",format(round(dt$lower[i],2),nsmall = 2),"-",format(round(dt$upper[i],2),nsmall = 2),")",sep="")
  dt$text[i] <- txt
  dt$time[i] <- ifelse(dt$time[i]==1,"at 1 year",paste0("at ",dt$time[i]," years"))
}
dt$group <- factor(dt$group,levels = c(1,0))
dt$time <- factor(dt$time,levels = c("at 1 year","at 3 years","at 5 years","at 10 years"))

p3 <- ggplot(data=dt,
       aes(x = group,y = rmst, ymin = lower, ymax = upper))+
  geom_pointrange(aes(col=group),size=0.4)+
  xlab('')+ ylab("RMST (95% Credible Interval)")+
  ylim(0,10)+
  scale_x_discrete(labels=NULL,expand = c(0.2,0.2))+
  #scale_y_continuous(expand=c(0,0))+
  geom_errorbar(aes(ymin=lower, ymax=upper,col=group),width=0.3,cex=1)+ 
  facet_wrap(~time,strip.position="top",nrow=4,scales = "free_y") +
  scale_color_manual(values = c("0"="darksalmon","1"="deepskyblue3"),
                     labels = c("MRD positive","MRD negative"))+
  geom_text(aes(y=9,label=text),size=4)+
  theme(plot.title=element_text(size=16,face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray98"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=16),
        axis.title=element_text(size=16),
        strip.text = element_text(hjust = 0,size=16),
        strip.background = element_blank(),
        legend.position = "top")+
  coord_flip()

pdf("figure1.pdf",width = 16,height = 6)
ggarrange(p,p1,p3,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()

