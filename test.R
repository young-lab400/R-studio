library("cluster")
library("MASS")
library("modeest")
library("clusterSim")
library("nlme")
library("analysis")
library(stats)

#自定義函數 
data_ks_test <- function(data)
{
    test_data <- rnorm(length(data), mean = 0, sd = 1)
    #檢定測試資料是否來自normal(Yes : p.value>0.05)
    data_result <- ks.test(unique(test_data),"pnorm")$p.value
    #檢定data是否與測試資料同分布(Yes : p.value>0.05)
    test_result <- ks.test(unique(data),unique(test_data))$p.value
    result <- c(data_result,test_result)
    return (result)
}

setwd("C:/Users/LAB401/Documents/data/Disease data")
#FolderA = "~/data/Disease data/CD"
FolderA = "./CD"
filesA = list.files(FolderA,pattern="*.txt")
CD_read = NULL
CD_normalized = NULL
for (x in 1:length(filesA))
{
  temppath = paste(FolderA,filesA[x],sep = "/")
  temppath
  CD_read <-read.table(temppath, header=F,sep="\t",stringsAsFactors=FALSE)
  #刪除末3年資料
  CD_read <- CD_read[,-3]
  CD_read <- CD_read[,-2]
  CD_read <- CD_read[,-1]
  #CD_read
  CD_read_numeric <- as.numeric(unlist(CD_read))
  CD_normalized[[x]] <- data.Normalization(CD_read_numeric,type="n1")
  #CD_normalized[[x]]
}


FolderB = "./DD"
filesB = list.files(FolderB,pattern="*.txt")
DD_read = NULL
DD_normalized = NULL
for(x in 1:length(filesB))
{
  temppath = paste(FolderB,filesB[x],sep="/")
  DD_read <- read.table(temppath, header=F,sep="\t",stringsAsFactors=FALSE)
  #刪除末3年資料
  DD_read <- DD_read[,-3]
  DD_read <- DD_read[,-2]
  DD_read <- DD_read[,-1]
  DD_read_numeric <- as.numeric(unlist(DD_read))
  DD_normalized[[x]] <- data.Normalization(DD_read_numeric,type = "n1")
  #DD_normalized[[x]]
}

FolderC = "./feature"
filesC = list.files(FolderC,pattern="*.txt")
f_read = NULL
f_normalized = NULL
for (x in 1:length(filesC))
{
  temppath = paste(FolderC,filesC[x],sep="/")
  #temppath
  f_read <-read.table(temppath, header=T,sep="\t",row.names=1,stringsAsFactors=FALSE,fileEncoding = "UCS-2LE")
  f_read[f_read=="."]<-NA
  f_read <- as.data.frame(t(f_read), stringsAsFactors=FALSE)
  #f_read
  #刪除方法去除缺失
  ppp <- complete.cases(f_read)
  ppp
  #月份正規化
  f_read <- f_read[complete.cases(f_read),]
  f_read[1,] <-f_read[1,]/31*30
  f_read[3,] <-f_read[3,]/31*30
  f_read[5,] <-f_read[5,]/31*30
  f_read[7,] <-f_read[7,]/31*30
  f_read[8,] <-f_read[8,]/31*30
  f_read[10,] <-f_read[10,]/31*30
  f_read[12,] <-f_read[12,]/31*30
  
  f_read_numeric <- as.numeric(unlist(f_read))
  #f_read_numeric
  f_normalized[[x]] <- data.Normalization(f_read_numeric,type="n1")
  #f_normalized[[x]]
}

#BL_test for CD 
for(x in 1:length(filesA))
{
  CD_BL_result = Box.test(CD_normalized[[x]],lag = 12,type = "Ljung-Box")
  CD_BL_result
  CD_BL_result_2 = data.frame(t(matrix(unlist(CD_BL_result))))
  names(CD_BL_result_2) <- names(CD_BL_result)
}
CD_BL_result
CD_BL_result[3]
CD_BL_result_2

#ks.test for CD
CD_result <- NULL
for (y in 1:length(filesA))
{
 #CD_result[x,y]<-ks.test(unique(f_normalized[[x]]),unique(CD_normalized[[y]]),alternative = "two.sided")$p.value
  CD_result[y]<-ks.test(unique(CD_normalized[[y]]),"pnorm")$p.value
}
names(CD_result) <- filesA

#ks.test for DD
DD_result <- NULL
for (y in 1:length(filesB))
{
    #DD_result[x,y]<-ks.test(unique(f_normalized[[x]]),unique(DD_normalized[[y]]),alternative = "two.sided")$p.value
    DD_result[y]<-ks.test(unique(DD_normalized[[y]]),"pnorm")$p.value
}
names(DD_result) <- filesB

#ks.test for feature
f_result <- NULL
for (x in 1:length(filesC))
{
  f_result[x]<-ks.test(unique(f_normalized[[x]]),"pnorm")$p.value
}
names(f_result) <- filesC
#顯示結果and <0.05的值
#CD_result
#length(CD_result[CD_result<0.05])
#DD_result
#length(DD_result[DD_result<0.05])
#f_result
#length(f_result[f_result<0.05])

#ks.test for data and test_data(One time)
rnorm_testresult <- NULL
for (CD_index in 1:length(CD_normalized))
{
  rnorm_testresult[[CD_index]] <- data_ks_test(unlist(CD_normalized[CD_index]))
}
names(rnorm_testresult) <- filesA
rnorm_testresult
