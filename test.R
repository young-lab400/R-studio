library("cluster")
library("clusterSim")
library("MASS")
library("modeest")
library("nlme")
library("analysis")
library(stats)
setwd("C:/Users/LAB401/Documents/data/Disease data")
FolderA = "~/data/Disease data/CD"
filesA = list.files(FolderA,pattern="*.txt")
CD_read = NULL
CD_normalized = NULL

for (x in 1:length(filesA))
{
  temppath = paste(FolderA,filesA[x],sep = "/")
  temppath
  CD_read <-read.table(temppath, header=F,sep="\t",stringsAsFactors=FALSE)
  #刪除末3年資料
  CD_read <- CD_read[,-14]
  CD_read <- CD_read[,-13]
  CD_read <- CD_read[,-12]
  CD_read
  CD_read_numeric <- as.numeric(unlist(CD_read))
  CD_normalized[[x]] <- data.Normalization(CD_read_numeric,type="n1")
  CD_normalized[[x]]
}

FolderB = "~/data/Disease data/DD"
filesB = list.files(FolderB,pattern="*.txt")
DD_read = NULL
DD_normalized = NULL
for(x in 1:length(filesB))
{
  temppath = paste(FolderB,filesB[x],sep="/")
  DD_read <- read.table(temppath, header=F,sep="\t",stringsAsFactors=FALSE)
  #刪除末3年資料
  DD_read <- DD_read[,-14]
  DD_read <- DD_read[,-13]
  DD_read <- DD_read[,-12]
  DD_read_numeric <- as.numeric(unlist(DD_read))
  DD_normalized[[x]] <- data.Normalization(DD_read_numeric,type = "n1")
  DD_normalized[[x]]
}

FolderC = "~/data/Disease data/feature"
filesC = list.files(FolderC,pattern="*.txt")
f_read = NULL
f_normalized = NULL
for (x in 1:length(filesC))
{
  temppath = paste(FolderC,filesC[x],sep="/")
  temppath
  f_read <-read.table(temppath, header=T,sep="\t",row.names=1,stringsAsFactors=FALSE,fileEncoding = "UCS-2LE")
  f_read[f_read=="."]<-NA
  f_read <- as.data.frame(t(f_read), stringsAsFactors=FALSE)
  f_read
  #刪除方法去除缺失
  ppp <- complete.cases(f_read)
  ppp
  f_read <- f_read[complete.cases(f_read),]
  f_read
  
  f_read_numeric <- as.numeric(unlist(f_read))
  f_read_numeric
  f_normalized[[x]] <- data.Normalization(f_read_numeric,type="n1")
  f_normalized[[x]]
}

#CD
CD_result <- data.frame()
for (x in 1:length(filesC))
{
  for (y in 1:length(filesA))
  {
    f_normalized[[x]]
    CD_normalized[[y]]
    CD_result[x,y]<-ks.test(unique(f_normalized[[x]]),unique(CD_normalized[[y]]),alternative = "two.sided")$p.value
  }
}
names(CD_result) <- filesA
row.names(CD_result) <- filesC
CD_result
#DD
DD_result <- data.frame()
for (x in 1:length(filesC))
{
  for (y in 1:length(filesB))
  {
    f_normalized[[x]]
    DD_normalized[[y]]
    DD_result[x,y]<-ks.test(unique(f_normalized[[x]]),unique(DD_normalized[[y]]),alternative = "two.sided")$p.value
  }
}
names(DD_result) <- filesB
row.names(DD_result) <- filesC
DD_result
