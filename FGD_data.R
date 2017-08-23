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
  #test_data <- rnorm(length(data), mean = 0, sd = 1)
  test_data <-NULL
  for (x in 1:length(data))
  {
    test_data <- c(test_data, rnorm(1, mean = 0, sd = 1))
  }
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
  f_read[2,]
  #二月份正規化 2000~2011
  for (year in 1:length(f_read[2,]))
  {
    year1 <- year+1999
    if (year1 %% 400==0)
    {
      f_read[2,year] <- f_read[2,year]/29*30
    }
    else if(year1 %% 4==0 & year1 %% 100!=0)
    {
      f_read[2,year] <- f_read[2,year]/29*30
    }
    else
    {
      f_read[2,year] <- f_read[2,year]/28*30
    }
  }
  
  #f_read_numeric
  f_read_numeric <- as.numeric(unlist(f_read))
  f_normalized[[x]] <- data.Normalization(f_read_numeric,type="n1")
  #f_normalized[[x]] <- f_read_numeric
  #f_normalized[[x]]
}

setwd("C:/Users/LAB401/Documents/data/Disease data/Normalization/CD")
for (x in 1:length(filesA))
{
  file_result <- data.frame(CD_normalized[x])
  row.names(file_result)<-NULL
  names(file_result)<-NULL
  write.table(file_result,filesA[x],row.names = FALSE,sep = "\n")
}

#DD_normalized
setwd("C:/Users/LAB401/Documents/data/Disease data/Normalization/DD")
for (x in 1:length(filesB))
{
  file_result <- data.frame(DD_normalized[x])
  row.names(file_result)<-NULL
  names(file_result)<-NULL
  write.table(file_result,filesB[x],row.names = FALSE,sep = "\n")
}

#f_normalized
setwd("C:/Users/LAB401/Documents/data/Disease data/Normalization/F")
for (x in 1:length(filesC))
{
  file_result <- data.frame(f_normalized[x])
  row.names(file_result)<-NULL
  names(file_result)<-NULL
  write.table(file_result,filesC[x],row.names = FALSE,sep = "\n")
}


#顯示結果and <0.05的值
#CD_result
#length(CD_result[CD_result<0.05])
#DD_result
#length(DD_result[DD_result<0.05])
#f_result
#length(f_result[f_result<0.05])



