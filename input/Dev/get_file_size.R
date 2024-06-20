library(dplyr)
library(tibble)

a1 <- read.table('/projectnb/planet/PLSP/runLogs/get_file_size.sh.o2615996')
a2 <- read.table('/projectnb/planet/PLSP/runLogs/get_file_size.sh.o2615997')

aa <- rbind(a1,a2)

bb <- matrix(NA,dim(aa)[1],4)
for(i in 1:dim(aa)[1]){
  temp1 <- unlist(strsplit(aa$V1[i],'kB'))
  temp2 <- unlist(strsplit(aa$V2[i],'/'))
  
  bb[i,1:length(temp2)] <- temp2
  bb[i,1             ] <- temp1
}
b1 <- which(is.na(bb[,4]))

# all files stoared
cc <- bb[-b1,]
cc <- as.data.frame(cc)



dat <- cc %>% 
  rownames_to_column() %>% 
  filter(!duplicated(V4))

sites <- unique(dat$V2)

dat1 <- dat[dat$V3==2016,]
dat2 <- dat[dat$V3==2017,]
dat3 <- dat[dat$V3==2018,]
dat4 <- dat[dat$V3==2019,]
dat5 <- dat[dat$V3==2020,]


sT1 <- sum(as.numeric(cc$V1))/(10^9)
sT2 <- sum(as.numeric(dat$V1))/(10^9)

sT3 <- sum(as.numeric(dat1$V1))/(10^9)
sT4 <- sum(as.numeric(dat2$V1))/(10^9)
sT5 <- sum(as.numeric(dat3$V1))/(10^9)
sT6 <- sum(as.numeric(dat4$V1))/(10^9)
sT7 <- sum(as.numeric(dat5$V1))/(10^9)


