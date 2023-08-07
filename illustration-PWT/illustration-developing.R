###################################################################################
## Inference for DEA Estimators of Malmquist Productivity Indices: 
## An Overview, Further Improvements and a Guide for Practitioners
## Author: Valentin Zelenyuk, Shirong Zhao
## Date: June 27, 2023
## The programming codes used in this paper involve 
## some earlier codes from Paul Wilson
## All rights reserved. 
## It is free for academic use only with adequate citation and acknowledgments.
## For any other use, contact the authors.
##################################################################################
rm(list=ls())
require(readxl)
require(rDEA)
source("./Functions/coverage.simple.R")
source("./Functions/coverage.agg.R")
#### Set Seed ####
if (exists(".Random.seed")) {
  save.seed=.Random.seed
  flag.seed=TRUE
} else {
  flag.seed=FALSE
}
set.seed(88888888)
################################################
np=2
nq=1
######################################################
pwt100 <- read_excel("Data/pwt100.xlsx", sheet = "Data")
######################################################
bhz<-c("Albania","Argentina","Armenia","Australia", "Austria", "Azerbaijan","Belarus",
       "Belgium", "Bolivia (Plurinational State of)","Brazil","Bulgaria","Canada","Chile",
       "China","Colombia","Costa Rica","Croatia","Czech Republic","Denmark","Dominican Republic",
       "Ecuador","Estonia","Finland","France","Germany","Greece","Guatemala","Honduras","China, Hong Kong SAR",
       "Hungary","Iceland","India","Indonesia" ,"Ireland","Israel","Italy","Jamaica","Japan",
       "Kazakhstan","Kenya","Republic of Korea","Kyrgyzstan","Latvia", "Lithuania","North Macedonia",
       "Madagascar","Malawi","Malaysia",
       "Mauritius","Mexico","Republic of Moldova","Morocco","Netherlands","New Zealand","Nigeria","Norway","Panama",
       "Paraguay","Peru","Philippines","Poland","Portugal","Romania","Russian Federation","Sierra Leone" ,
       "Singapore","Slovakia","Slovenia","Spain","Sri Lanka","Sweden","Switzerland","Syrian Arab Republic",
       "Taiwan", "Tajikistan","Thailand","Turkey","Ukraine","United Kingdom","Uruguay","United States","Venezuela (Bolivarian Republic of)",
       "Zambia","Zimbabwe")

developed<-c("Australia","Austria","Belgium","Canada","Denmark","Finland","France","Germany","Greece","China, Hong Kong SAR",
             "Iceland","Ireland","Israel","Italy","Japan","Republic of Korea","Netherlands","New Zealand","Norway","Portugal","Singapore","Spain","Sweden",
             "Switzerland","Taiwan","United Kingdom","United States")


ii=which(pwt100$country %in% bhz)
df<-pwt100[ii,]

ii=which(df$year>=1990)
df<-df[ii,]

ii=which(df$country %in% developed)
dfa<-df[ii,]
dfb<-df[-ii,]


####### developing countries ########

year=c(1990,1995,2000,2005,2010,2015,2019)
res.agg=matrix(NA,nrow=length(year),ncol=28)
res.simple=matrix(NA,nrow=length(year),ncol=28)
res.agg.sd=matrix(NA,nrow=length(year),ncol=4)
res.simple.sd=matrix(NA,nrow=length(year),ncol=4)

for (i in 1:length(year)) {
  
  cat (i,"\n")
  
  if (i<length(year)){
    
    i1=which(dfb$year==year[i])
    df1b=dfb[i1,]
    x1=cbind(df1b$emp,df1b$cn)
    y1=matrix(df1b$rgdpo,ncol=1)
    i2=which(dfb$year==year[i+1])
    df2b=dfb[i2,]
    x2=cbind(df2b$emp,df2b$cn)
    y2=matrix(df2b$rgdpo,ncol=1)
    
  } else {
    
    i1=which(dfb$year==year[1])
    df1b=dfb[i1,]
    x1=cbind(df1b$emp,df1b$cn)
    y1=matrix(df1b$rgdpo,ncol=1)
    i2=which(dfb$year==year[i])
    df2b=dfb[i2,]
    x2=cbind(df2b$emp,df2b$cn)
    y2=matrix(df2b$rgdpo,ncol=1)
    
  }
  
  tt=coverage.agg(x1=x1,y1=y1,x2=x2,y2=y2)
  
  res.agg[i,1:2]=exp(tt$estimate[1:2])
  res.agg[i,3:4]=exp(tt$estimate[5:6])
  
  # 90%
  res.agg[i,5:6]=tt$bounds.psz[1,]
  res.agg[i,7:8]=tt$bounds.szz[1,]
  res.agg[i,9:10]=tt$bounds.nsz[1,]
  res.agg[i,11:12]=tt$bounds.nsz.szz[1,]
  
  # 95%
  res.agg[i,13:14]=tt$bounds.psz[2,]
  res.agg[i,15:16]=tt$bounds.szz[2,]
  res.agg[i,17:18]=tt$bounds.nsz[2,]
  res.agg[i,19:20]=tt$bounds.nsz.szz[2,]
  
  # 99%
  res.agg[i,21:22]=tt$bounds.psz[3,]
  res.agg[i,23:24]=tt$bounds.szz[3,]
  res.agg[i,25:26]=tt$bounds.nsz[3,]
  res.agg[i,27:28]=tt$bounds.nsz.szz[3,]
  
  #
  res.agg.sd[i,]=tt$sig
  
  ###############
  ss=coverage.simple(x1=x1,y1=y1,x2=x2,y2=y2)
  
  res.simple[i,1:2]=exp(ss$estimate[1:2])
  res.simple[i,3:4]=exp(ss$estimate[5:6])
  
  # 90%
  res.simple[i,5:6]=ss$bounds.ksw[1,]
  res.simple[i,7:8]=ss$bounds.szz[1,]
  res.simple[i,9:10]=ss$bounds.nsz[1,]
  res.simple[i,11:12]=ss$bounds.nsz.szz[1,]
  
  # 95%
  res.simple[i,13:14]=ss$bounds.ksw[2,]
  res.simple[i,15:16]=ss$bounds.szz[2,]
  res.simple[i,17:18]=ss$bounds.nsz[2,]
  res.simple[i,19:20]=ss$bounds.nsz.szz[2,]
  
  # 99%
  res.simple[i,21:22]=ss$bounds.ksw[3,]
  res.simple[i,23:24]=ss$bounds.szz[3,]
  res.simple[i,25:26]=ss$bounds.nsz[3,]
  res.simple[i,27:28]=ss$bounds.nsz.szz[3,]
  
  #
  res.simple.sd[i,]=ss$sig
}


res.agg.sig=matrix(NA,nrow=length(year),ncol=12)
res.simple.sig=matrix(NA,nrow=length(year),ncol=12)

res.agg.sig[,c(1,2)]=res.agg[,c(1,3)]
res.agg.sig[,c(3,4)]=log(res.agg[,c(1,3)])
res.agg.sig[,5:8]=res.agg.sd
res.agg.sig=formatC(res.agg.sig,width=6,digits = 4,format = "f")

res.simple.sig[,c(1,2)]=res.simple[,c(1,3)]
res.simple.sig[,c(3,4)]=log(res.simple[,c(1,3)])
res.simple.sig[,5:8]=res.simple.sd
res.simple.sig=formatC(res.simple.sig,width=6,digits = 4,format = "f")

for (i in 1:dim(res.agg)[1]) {
  
  row.i=res.agg[i,]
  
  if (row.i[21]>0 | 0>row.i[22]){
    res.agg.sig[i,9]="***"
  } else if (row.i[13]>0 | 0>row.i[14]){
    res.agg.sig[i,9]="** "
  } else if  (row.i[5]>0 | 0>row.i[6]){
    res.agg.sig[i,9]="*  "
  } else {
    res.agg.sig[i,9]="-- "
  }
  
  if (row.i[23]>0 | 0>row.i[24]){
    res.agg.sig[i,10]="***"
  } else if (row.i[15]>0 | 0>row.i[16]){
    res.agg.sig[i,10]="** "
  } else if  (row.i[7]>0 | 0>row.i[8]){
    res.agg.sig[i,10]="*  "
  } else {
    res.agg.sig[i,10]="-- "
  }
  
  
  if (row.i[25]>0 | 0>row.i[26]){
    res.agg.sig[i,11]="***"
  } else if (row.i[17]>0 | 0>row.i[18]){
    res.agg.sig[i,11]="** "
  } else if  (row.i[9]>0 | 0>row.i[10]){
    res.agg.sig[i,11]="*  "
  } else {
    res.agg.sig[i,11]="-- "
  }
  
  if (row.i[27]>0 | 0>row.i[28]){
    res.agg.sig[i,12]="***"
  } else if (row.i[19]>0 | 0>row.i[20]){
    res.agg.sig[i,12]="** "
  } else if  (row.i[11]>0 | 0>row.i[12]){
    res.agg.sig[i,12]="*  "
  } else {
    res.agg.sig[i,12]="-- "
  }
  
}

# for simple mean MPI
for (i in 1:dim(res.simple)[1]) {
  
  row.i=res.simple[i,]
  
  if (row.i[21]>0 | 0>row.i[22]){
    res.simple.sig[i,9]="***"
  } else if (row.i[13]>0 | 0>row.i[14]){
    res.simple.sig[i,9]="** "
  } else if  (row.i[5]>0 | 0>row.i[6]){
    res.simple.sig[i,9]="*  "
  } else {
    res.simple.sig[i,9]="-- "
  }
  
  if (row.i[23]>0  | 0>row.i[24]){
    res.simple.sig[i,10]="***"
  } else if (row.i[15]>0 | 0>row.i[16]){
    res.simple.sig[i,10]="** "
  } else if  (row.i[7]>0 | 0>row.i[8]){
    res.simple.sig[i,10]="*  "
  } else {
    res.simple.sig[i,10]="-- "
  }
  
  
  if (row.i[25]>0  | 0>row.i[26]){
    res.simple.sig[i,11]="***"
  } else if (row.i[17]>0 | 0>row.i[18]){
    res.simple.sig[i,11]="** "
  } else if  (row.i[9]>0 | 0>row.i[10]){
    res.simple.sig[i,11]="*  "
  } else {
    res.simple.sig[i,11]="-- "
  }
  
  if (row.i[27]>0  | 0>row.i[28]){
    res.simple.sig[i,12]="***"
  } else if (row.i[19]>0 | 0>row.i[20]){
    res.simple.sig[i,12]="** "
  } else if  (row.i[11]>0 | 0>row.i[12]){
    res.simple.sig[i,12]="*  "
  } else {
    res.simple.sig[i,12]="-- "
  }
  
}

res.simple.sig
res.agg.sig

