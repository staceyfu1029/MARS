rm(list=ls())
setwd("path")

txt <- "Routput/sim2/"

rd.list <- list.files(txt)
rd.list <- rd.list[grep("RData",rd.list)]

for(i in 1:length(rd.list)){
  load(paste0(txt,rd.list[i]))
  temp <- rec.out
  if (i==1) {
    out <- rbind(temp)
  } else {
    out <- rbind(out, temp)
  }
}



  

