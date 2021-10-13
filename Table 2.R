rm(list=ls())
## Simulation result table 2

setwd("path")

## table
t <- data.frame(matrix(NA,nrow = 17,ncol = 16))
t[1,c(3,10)] <- c("x=0","x=1")
t[2,c(3,7,10,14)] <- rep(c("MARS","IPD"),2)
t[3,3:16] <- rep(c("True","Average(SD)","MSE","CP","Average(SD)","MSE","CP"),2)
t[c(4,11),1] <- c("Case 1","Case 2")
t[4:17,2] <- rep(c("SR at 3 yrs","SR at 5 yrs","SR at 10 yrs",
                   "MS","RMST at 5 yrs","RMST Diff","RMST Ratio"),2)
#t[4:17,3] <- c("0.549","0.461","0.334","48.05","36.42","9.00","1.247",
#               "0.703","0.549","0.301","69.31","45.12","5.72","1.127")
#t[4:17,10] <- c("0.719","0.654","0.548","-","45.42","-","-",
#                "0.835","0.540","0.116","63.28","50.84","-","-")
t[4:17,3] <- c("0.55","0.46","0.33","4.00","3.04","0.75","1.25",
               "0.70","0.55","0.30","5.78","3.76","0.48","1.13")
t[4:17,10] <- c("0.72","0.65","0.55","-","3.76","-","-",
                "0.84","0.54","0.12","5.27","4.24","-","-")

tv <- c("0.549","0.461","0.334","48.05","36.42","9.00","1.247","0.703","0.549","0.301","69.31","45.12","5.72","1.127")
tv <- cbind(tv,c("0.719","0.654","0.548","-","45.42","-","-","0.835","0.540","0.116","63.28","50.84","-","-"))

## Scenario 1
load("Routput/sim1.RData")

m <- c("hat","lower","upper")
keyword.list <- c("SR.3.0","SR.5.0","SR.10.0","MS.0","RMST.5.0","RMSTD.5","RMSTR.5")

get_cell <- function(t,keyword.list,group,method="MARS",scene=1){
  if(group==1){
    cl <- 10
  }else{
    cl <- 3
  }
  if(method=="t1"){
    cl2 <- cl+3
  }else{
    cl2 <- cl
  }
  if(scene==1){
    nrow <- 3
  }else{
    nrow <- 10
  }
  for(i in 1:length(keyword.list)){
    if(method=="t1"){
      w <- paste(keyword.list[i],m,"t1",sep=".") 
    }else{
      w <- paste(keyword.list[i],m,sep=".")
    }
    temp <- out[,w]
    tv.x <- as.numeric(tv[i,1+(cl==10)])
    if(is.na(tv.x)){
      t[i+nrow,1:3+cl2] <- "-"
    }else{
      ave <- format(round(mean(temp[,1]),digits = 2),nsmall=2)
      var <- format(round(sd(temp[,1]),digits = 2),nsmall=2)
      t[i+nrow,1+cl2] <- paste0(ave," (",var,")")
      t[i+nrow,2+cl2] <- round(mean((temp[,1]-tv.x)^2),digits = 4)
      t[i+nrow,3+cl2] <- format(round(mean(temp[,2]<tv.x & temp[,3]>tv.x),digits = 2),nsmall=2)
    }
  }
  return(t)
}

t <- get_cell(t=t,keyword.list = keyword.list,group = 0)
t <- get_cell(t=t,keyword.list = keyword.list,group = 0,method = "t1")
keyword.list_2 <- c("SR.3.1","SR.5.1","SR.10.1","MS.1","RMST.5.1","RMSTD.5","RMSTR.5")
t <- get_cell(t=t,keyword.list = keyword.list_2,group = 1)
t <- get_cell(t=t,keyword.list = keyword.list_2,group = 1,method = "t1")

# Scenario 2
load("Routput/sim2.RData")
t <- get_cell(t=t,keyword.list = keyword.list,group = 0,scene = 2)
t <- get_cell(t=t,keyword.list = keyword.list,group = 0,method = "t1",scene = 2)
keyword.list_2 <- c("SR.3.1","SR.5.1","SR.10.1","MS.1","RMST.5.1")
t <- get_cell(t=t,keyword.list = keyword.list_2,group = 1,scene = 2)
t <- get_cell(t=t,keyword.list = keyword.list_2,group = 1,method = "t1",scene = 2)

writexl::write_xlsx(t,"Table/table 2.xlsx",col_names = F)
