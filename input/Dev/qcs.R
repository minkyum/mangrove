library(raster)
library(rgdal)
library(gdalUtils)
library(ncdf4)

library(rjson)
library(geojsonR)

library(RColorBrewer)


###############################
args <- commandArgs()
print(args)

numSite <- as.numeric(substr(args[3],1,3))
yy      <- as.numeric(substr(args[3],4,7))
# numSite <- 1; yy <- 2017


########################################
workPath  <- '/projectnb/planet/PLSP/'
sites     <- list.files(paste0(workPath,'Product_netCDF_fillVal'))
dirs      <- list.files(paste0(workPath,'Product_netCDF_fillVal'),full.names=T)

ptFiles <- list.files(paste0(dirs[numSite]),pattern=glob2rx('*.nc'),full.names=T)
nc       <- nc_open(ptFiles[1])
varNames <- names(nc$var)


## All layers as raster
qcOutDir <- paste0(workPath,'Output_images/fillVal')
if (!dir.exists(qcOutDir)) {dir.create(qcOutDir)}

  png(filename=paste0(qcOutDir,'/1_AllRast_',sites[numSite],'_',yy,'.png'),width=18,height=8.5,units='in',res=100)
  par(mfrow=c(4,6),oma=c(0,0,0,2),mar=c(0,1,2,4))

  for(i in 2:25){

    rast <- raster(ptFiles[(yy-2016)],varname=varNames[i])
    pcFill <- round(sum(values(rast)==32767)/length(rast)*100,4)
    rast[rast==32767] <- NA

    if(i==2){
      plot(rast,axes=F,box=F,
           main=paste0(varNames[i],'_',pcFill,'%'),
           colNA='grey45',cex.main=1.3)
    }else{
      plot(rast,axes=F,box=F,
           main=paste0(varNames[i],'_',pcFill,'%'),
           colNA='grey45',cex.main=1.3)
    }
    print(i)
  }

  dev.off()






########################################
#
spec <- rev(brewer.pal(11,'Spectral'))
mycolRamp = colorRampPalette(c('White',spec))


# for(yy in 2017:2020){
  ptDir1 <- paste0(workPath,'Archive/Product_010/',sites[numSite],'/',yy)
  ptFiles1 <- list.files(ptDir1,pattern='*.tif',full.names=T)
  ptFile <- ptFiles[(yy-2016)]
  
  metTable <- read.csv('/projectnb/modislc/users/mkmoon/Planet/data_paper/data/site_full_w_pc_hls_PCrevised_ext.csv')

  png(filename=paste0(qcOutDir,'/2_Comp_',sites[numSite],'_',yy,'.png'),width=10,height=6,units='in',res=150)
  for(i in c(4,8)){
    
    rast1 <- raster(ptFiles1[i-1])
    rast2 <- raster(ptFile,varname=varNames[i])
    pcFill <- round(sum(values(rast2)==32767)/length(rast2)*100,4)
    rast1[rast1==32767] <- NA
    rast2[rast2==32767] <- NA
    
    if(i==4){
      x1 <- values(rast1)
      y1 <- values(rast1)
      par(fig=c(0,0.3,0.5,1),oma=c(2,2,2,2),mar=c(2,2,1,2),mgp=c(2,1,0))
      smoothScatter(x1,y1,xlim=c(-100,500),ylim=c(-100,500),
                    colramp=mycolRamp,nbin=850,nrpoints=0,
                    transformation=function(x)x^.6,axe=F,
                    xlab='Proto',ylab='3-day',main=compareRaster(rast1,rast2,res=T,orig=T),
                    cex.lab=1)
      axis(1,at=seq(-100,500,100),cex.axis=1)
      axis(2,at=seq(-100,500,100),cex.axis=1)
      abline(0,1,lty=5)

      qt02  <- quantile(values(rast1),0.05,na.rm=T)
      qt98  <- quantile(values(rast1),0.95,na.rm=T)

      rast1[rast1<qt02] <- qt02
      rast1[rast1>qt98] <- qt98

      par(fig=c(0.3,0.65,0.5,1),oma=c(0,0,0,2),mar=c(0,1,2,4),new=T)
      plot(rast1,axes=F,box=F,
           main=paste0(varNames[i],'_',yy,'_20'),
           colNA='grey45',cex.main=1)

      rast2[rast2<qt02] <- qt02
      rast2[rast2>qt98] <- qt98

      par(fig=c(0.65,1,0.5,1),oma=c(0,0,0,2),mar=c(0,1,2,4),new=T)
      plot(rast2,axes=F,box=F,
           main=paste0(varNames[i],'_',pcFill,'%'),
           colNA='grey45',cex.main=1)
    }else{
      # x1 <- values(rast1)
      # y1 <- values(rast1)
      # par(fig=c(0,0.3,0,0.5),oma=c(2,2,2,2),mar=c(2,2,1,2),mgp=c(2,1,0),new=T)
      # smoothScatter(x1,y1,xlim=c(-100,500),ylim=c(-100,500),
      #               colramp=mycolRamp,nbin=850,nrpoints=0,
      #               transformation=function(x)x^.6,axe=F,
      #               xlab='Proto',ylab='3-day',
      #               cex.lab=1)
      # axis(1,at=seq(-100,500,100),cex.axis=1)
      # axis(2,at=seq(-100,500,100),cex.axis=1)
      # abline(0,1,lty=5)
      
      par(fig=c(0,0.3,0,0.5),oma=c(0,0,0,2),mar=c(0,1,2,4),new=T)
      waterDir <- paste0(workPath,'Img_cliped/',metTable$Site[numSite],'/water_mask_30_1.tif')

      rastWater <- raster(waterDir)
      plot(rastWater,axes=F,box=F,zlim=c(0,100),
           # main=paste0(productTable$short_name[i],'_',yy,'_gap'),
           main='water',
           colNA='grey45',cex.main=1)

      qt02  <- quantile(values(rast1),0.05,na.rm=T)
      qt98  <- quantile(values(rast1),0.95,na.rm=T)

      rast1[rast1<qt02] <- qt02
      rast1[rast1>qt98] <- qt98

      par(fig=c(0.3,0.65,0,0.5),oma=c(0,0,0,2),mar=c(0,1,2,4),new=T)
      plot(rast1,axes=F,box=F,
           main=paste0(varNames[i],'_',yy,'_20'),
           colNA='grey45',cex.main=1)

      rast2[rast2<qt02] <- qt02
      rast2[rast2>qt98] <- qt98

      par(fig=c(0.65,1,0,0.5),oma=c(0,0,0,2),mar=c(0,1,2,4),new=T)
      plot(rast2,axes=F,box=F,
           main=paste0(varNames[i],'_',pcFill,'%'),
           colNA='grey45',cex.main=1)
    }
  }
  dev.off()
# }



# ########################################
# # All images versus 3-day
# qcOutDir <- paste0(params$setup$workDir,'Output_images/gap')
# if (!dir.exists(qcOutDir)) {dir.create(qcOutDir)}
# 
# # for(yy in 2017:2020){
# 
# # yy <- 2020
# 
#   ptDir2 <- paste0(params$setup$workDir,'Product/gap_',strSite,'/',yy)
#   print(ptDir2)
# 
#   ptFiles2 <- list.files(ptDir2,pattern='*.tif',full.names=T)
#   print(length(ptFiles2))
# 
#   png(filename=paste0(qcOutDir,'/3_gap_',strSite,'_',yy,'.png'),width=12,height=6,units='in',res=150)
#   par(mfrow=c(2,3),oma=c(0,0,0,2),mar=c(0,1,2,4))
#   for(i in c(1,12,10,24,3,7)){
#     if(i==1|i==12){
#       rast2 <- raster(ptFiles2[i])
# 
#       plot(rast2,axes=F,box=F,
#            main=paste0(productTable$short_name[i],'_',yy,'_gap'),
#            colNA='grey45',cex.main=1)
#     }else if(i==24){
#       waterDir <- paste0(params$setup$workDir,'Img_cliped/',strSite,'/water_mask.tif')
#       
#       rast2 <- raster(waterDir)
#       plot(rast2,axes=F,box=F,zlim=c(0,100),
#            # main=paste0(productTable$short_name[i],'_',yy,'_gap'),
#            main='water',
#            colNA='grey45',cex.main=1)
#     }else{
#       rast2 <- raster(ptFiles2[i])
# 
#       qt02  <- quantile(values(rast2),0.02,na.rm=T)
#       qt98  <- quantile(values(rast2),0.98,na.rm=T)
# 
#       rast2[rast2<qt02] <- qt02
#       rast2[rast2>qt98] <- qt98
# 
#       plot(rast2,axes=F,box=F,
#            main=paste0(productTable$short_name[i],'_',yy,'_gap'),
#            colNA='grey45',cex.main=1)
#     }
#   }
#   dev.off()
#   print(yy)
# # }
  
 
  
  