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
setwd("C:/Users/LAB401/Documents/data/Disease data/test_output")
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
CD_ks_times <- NULL
for(smile in 1:100)
{
  rnorm_testresult <- NULL
  for (CD_index in 1:length(CD_normalized))
  {
    rnorm_testresult[[CD_index]] <- data_ks_test(unlist(CD_normalized[CD_index]))
  }
  names(rnorm_testresult) <- filesA
  rnorm_testresult <- data.frame(rnorm_testresult)
  CD_ks_times[smile] <- length(rnorm_testresult[2,rnorm_testresult[2,]<0.05])
}

DD_ks_times <- NULL
for(smile in 1:100)
{
  rnorm_testresult <- NULL
  for (DD_index in 1:length(DD_normalized))
  {
    rnorm_testresult[[DD_index]] <- data_ks_test(unlist(DD_normalized[DD_index]))
  }
  names(rnorm_testresult) <- filesB
  rnorm_testresult <- data.frame(rnorm_testresult)
  DD_ks_times[smile] <- length(rnorm_testresult[2,rnorm_testresult[2,]<0.05])
}

f_ks_times <- NULL
for(smile in 1:100)
{
  rnorm_testresult <- NULL
  for (f_index in 1:length(f_normalized))
  {
    rnorm_testresult[[f_index]] <- data_ks_test(unlist(f_normalized[f_index]))
  }
  names(rnorm_testresult) <- filesC
  rnorm_testresult <- data.frame(rnorm_testresult)
  f_ks_times[smile] <- length(rnorm_testresult[2,rnorm_testresult[2,]<0.05])
}

#CD_ks_times
#DD_ks_times
#f_ks_times
#ending test

#BL_test for CD 
CD_BL_result <- NULL
for(x in 1:length(filesA))
{
  CD_BL_result[[x]] <- Box.test(CD_normalized[[x]],lag = 12,type = "Ljung-Box")
  CD_BL_result[[x]]
  #CD_BL_result_2 = data.frame(t(matrix(unlist(CD_BL_result))))
  #names(CD_BL_result_2) <- names(CD_BL_result)
}
sink("CD_BL_test.txt",append = FALSE)
CD_BL_result
sink()

#BL_test for DD 
DD_BL_result <- NULL
for(x in 1:length(filesB))
{
  DD_BL_result[[x]] <- Box.test(DD_normalized[[x]],lag = 12,type = "Ljung-Box")
  DD_BL_result[[x]]
}
DD_BL_result
sink("DD_BL_test.txt",append = FALSE)
DD_BL_result
sink()

#BL_test for f 
f_BL_result <- NULL
for(x in 1:length(filesC))
{
  f_BL_result[[x]] <- Box.test(f_normalized[[x]],lag = 12,type = "Ljung-Box")
  f_BL_result[[x]]
}
sink("f_BL_test.txt",append = FALSE)
f_BL_result
sink()

#shapiro.test()
CD_SP_result <- NULL
for (x in 1:length(filesA))
{
  CD_SP_result[[x]] <-shapiro.test(CD_normalized[[x]])$p.value
  #CD_SP_result[[x]] <- unlist(CD_SP_result[[x]])
}

CD_SP_result <- data.frame(CD_SP_result)
row.names(CD_SP_result) <-filesA
sink("CD_SP_test.txt",append = FALSE)
CD_SP_result
CD_SP_result[CD_SP_result<0.05]
sink()

DD_SP_result <- NULL
for (x in 1:length(filesB))
{
  DD_SP_result[[x]] <-shapiro.test(DD_normalized[[x]])$p.value
  #CD_SP_result[[x]] <- unlist(CD_SP_result[[x]])
}

DD_SP_result <- data.frame(DD_SP_result)
row.names(DD_SP_result) <-filesB
sink("DD_SP_test.txt",append = FALSE)
DD_SP_result
DD_SP_result[DD_SP_result<0.05]
sink()

f_SP_result <- NULL
for (x in 1:length(filesC))
{
  f_SP_result[[x]] <-shapiro.test(f_normalized[[x]])$p.value
  #CD_SP_result[[x]] <- unlist(CD_SP_result[[x]])
}
f_SP_result <- data.frame(f_SP_result)
row.names(f_SP_result) <-filesC
f_SP_result
row.names(f_SP_result)
sink("f_SP_test.txt",append = FALSE)
f_SP_result
f_SP_result[f_SP_result<0.05]
sink()
