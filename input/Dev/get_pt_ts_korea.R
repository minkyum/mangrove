library(sp)
library(raster)
library(terra)

library(rjson)
library(geojsonR)

library(doMC)
library(doParallel)

###############################
args <- commandArgs()
print(args)

vv <- args[3]
# vv <- 'dbf_class5'
numSite <- 37

print(vv)


########################################
params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/PLSP_Parameters.json')
source(params$setup$rFunctions)


########################################
## Get site name and image directory
geojsonDir <- params$setup$geojsonDir

strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

# ckDir <- paste0(params$setup$outDir,strSite,'/chunk')
ckDir <- paste0('/projectnb/modislc/users/mkmoon/Planet/rawImage/chunks/',strSite)
print(ckDir)



########################################
##
imgBase <- raster(paste0(params$setup$outDir,strSite,'/base_image.tif'))
numPix <- length(imgBase)
imgNum <- setValues(imgBase,1:numPix)
numChunks <- params$setup$numChunks
chunk <- numPix%/%numChunks

# nn <- sprintf('%03d',numSite)
# ptShp <- shapefile(paste0(params$setup$workDir,'Shape_file/pt_',nn,'.shp'))
ptShp <- shapefile(paste0('/projectnb/modislc/users/mkmoon/KoreaProject/shpfile/points/d2/',vv,'.shp'))
ptShp <- spTransform(ptShp,crs(imgBase))

pixNums <- extract(imgNum,ptShp)


##
pixNum  <- pixNums[1]
ckNum <- sprintf('%03d',(pixNum%/%chunk+1))
file <- list.files(path=ckDir,pattern=glob2rx(paste0('*',ckNum,'.rda')),full.names=T)
load(file)

dat <- matrix(NA,length(dates),1000)
# for(i in 1:length(pixNums)){
registerDoMC()
dat <- foreach(i=1:1000, .combine=cbind) %dopar%{
  pixNum  <- pixNums[i]

  ckNum <- sprintf('%03d',(pixNum%/%chunk+1))
  file <- list.files(path=ckDir,pattern=glob2rx(paste0('*',ckNum,'.rda')),full.names=T)

  load(file)

  blue  <- band1[pixNum%%chunk,]
  green <- band2[pixNum%%chunk,]
  red   <- band3[pixNum%%chunk,]
  nir   <- band4[pixNum%%chunk,]
  phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr

  blue <- blue/10000; green <- green/10000; red <- red/10000; nir <- nir/10000

  pheno_pars <- params$phenology_parameters
  qa_pars    <- params$qa_parameters

  vi     <- 2.5*(nir - red) / (nir + 2.4*red + 1)
  # dat[,i] <- vi
  #
  # print(i)
}
print(dim(dat))

dat <- data.frame(dates,dat)

write.csv(dat,file=paste0('/projectnb/modislc/users/mkmoon/KoreaProject/data/pt_ts/d2/',vv,'.csv'))

# ########################################
# path <- '/projectnb/modislc/users/mkmoon/KoreaProject/data/pt_ts/'
# files <- list.files(path,pattern=glob2rx('*csv'),full.names=T)
# 
# dat <- read.csv(files[1])
# d1 <- rowMeans(dat[,3:1002])
# dat <- read.csv(files[2])
# d2 <- rowMeans(dat[,3:1002])
# dat <- read.csv(files[3])
# d3 <- rowMeans(dat[,3:1002])
# dat <- read.csv(files[4])
# d4 <- rowMeans(dat[,3:1002])
# dat <- read.csv(files[5])
# d5 <- rowMeans(dat[,3:1002])
# 
# dates <- dat[,2]
# 
# 
# ##
# dat <- read.csv(files[1+5])
# e1 <- rowMeans(dat[,3:1002])
# dat <- read.csv(files[2+5])
# e2 <- rowMeans(dat[,3:1002])
# dat <- read.csv(files[3+5])
# e3 <- rowMeans(dat[,3:1002])
# dat <- read.csv(files[4+5])
# e4 <- rowMeans(dat[,3:1002])
# dat <- read.csv(files[5+5])
# e5 <- rowMeans(dat[,3:1002])
# 
# dates <- dat[,2]
# 
# 
# setwd('/projectnb/modislc/users/mkmoon/KoreaProject/figure/')
# png(filename='pt_ts_d2.png',width=14,height=8,unit='in',res=300)
# 
# par(mfrow=c(2,1),oma=c(1,1,2,2),mar=c(4,4,1,1),mgp=c(2.7,0.8,0))
# plot(as.Date(dates),d1,pch=19,ylim=c(0,0.7),
#      xlab='Dates',ylab='EVI2')
# points(as.Date(dates),d2,col=2,pch=19)
# points(as.Date(dates),d3,col=3,pch=19)
# points(as.Date(dates),d4,col=4,pch=19)
# points(as.Date(dates),d5,col=5,pch=19)
# legend('topleft',c('C1','C2','C3','C4','C5'),
#        col=1:5,pch=19,cex=1.2,bty='n',title='DBF')
# 
# plot(as.Date(dates),e1,pch=19,ylim=c(0,0.7),
#      xlab='Dates',ylab='EVI2')
# points(as.Date(dates),e2,col=2,pch=19)
# points(as.Date(dates),e3,col=3,pch=19)
# points(as.Date(dates),e4,col=4,pch=19)
# points(as.Date(dates),e5,col=5,pch=19)
# legend('topleft',c('C1','C2','C3','C4','C5'),
#        col=1:5,pch=19,cex=1.2,bty='n',title='ENF')
# 
# dev.off()