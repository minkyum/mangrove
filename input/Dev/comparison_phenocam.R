library(raster)
library(rgdal)
library(gdalUtils)
library(rgeos)

library(rjson)

library(RColorBrewer)
library(viridis)

########################################
args <- commandArgs()
print(args)

ss <- as.numeric(substr(args[3],1,2))
# ss <- 3

vari <- c('50PCGI','50PCGD')


########################################
Rast <- vector('list',8); ind <- 1
for(vv in 1:2){
  for(yy in 2017:2020){
    imgs <- list.files('/projectnb/planet/PLSP/Product',pattern=glob2rx(paste0('*',yy,'_',vari[vv],'.tif')),full.names=T,recursive=T)  
    
    Rast[[ind]] <- raster(imgs[ss])    
    
    print(ind)
    ind  <- ind + 1 
  }
}

siteStr <- unlist(strsplit(imgs[ss],'/'))[6]
print(siteStr)


siteFull <- read.csv('/projectnb/modislc/users/mkmoon/Planet/data_paper/data/site_full_w_pc_hls.csv')

afcd <- siteFull$Site_amerifulx[which(siteFull$Site==substr(siteStr,5,80))]
vgty <- siteFull$Veg[which(siteFull$Site==substr(siteStr,5,80))]

print(afcd)
print(vgty)




########################################
# Get site name, image directory and coordinate
rast <- Rast[[1]]
cx <- (rast@extent@xmax+rast@extent@xmin)/2
cy <- (rast@extent@ymax+rast@extent@ymin)/2
  
# Point Shape file
utm_crs = rast@crs
site <- data.frame(1,cx,cy)
colnames(site) <- c('id','xcrd','ycrd')
xy   <- site[,c(2,3)]
bb   <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=utm_crs)

geog_crs = CRS("+proj=longlat +datum=WGS84")
bb   <- spTransform(bb,geog_crs)
  
  
# Find associate PhenoCam site
pcMeta <- list.files('/projectnb/planet/PhenoCam/PhenoCam_V3/data_record_1',pattern='*_meta.json',full.names=T)

siteWPC <- c()
for(i in 1:length(pcMeta)){
  meta <- fromJSON(file=pcMeta[i])
  
  # Point Shape file
  geog_crs = CRS("+proj=longlat +datum=WGS84")
  site <- data.frame(1,meta$phenocam_site$lon,meta$phenocam_site$lat)
  colnames(site) <- c('id','xcrd','ycrd')
  xy   <- site[,c(2,3)]
  bb   <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=geog_crs)

  utm_crs = rast@crs
  bb   <- spTransform(bb,utm_crs)
  
  intPt <- intersect(bb,rast)  

  if(length(intPt)>0){siteWPC <- c(siteWPC,meta$phenocam_site$sitename)} 
  
  if(i%%100==0) print(i)
}


print(length(siteWPC))


# Get phenometric
if(length(siteWPC)>0){
 
  
  for(i in 1:length(siteWPC)){
    datMetrics <- matrix(NA,14,6)
    
    fJson <- list.files('/projectnb/planet/PhenoCam/PhenoCam_V3/data_record_1',pattern=glob2rx(paste0(siteWPC[i],'_*json')),full.names=T)
    fMetr <- list.files('/projectnb/planet/PhenoCam/PhenoCam_V3/data_record_5',pattern=glob2rx(paste0(siteWPC[i],'_*3day_transition_dates.csv')),full.names=T)
    
    # if(length())
    
    meta <- fromJSON(file=fJson)
    datP <- c()
    for(j in 1:length(fMetr)){
      tempP <- read.csv(fMetr[j],skip=16)  
      datP  <- rbind(datP,tempP)
    }
    
    
    ## Get PhenoCam data
    dat1 <- datP[datP$direction=='rising'&datP$gcc_value=='gcc_90',8]
    dat2 <- datP[datP$direction=='falling'&datP$gcc_value=='gcc_90',8]
    
    if(length(which(substr(dat1,1,4)==2016))>0){datMetrics[1,1] <- as.numeric(strftime(dat1[which(substr(dat1,1,4)==2016)][1], format = "%j"))}
    if(length(which(substr(dat1,1,4)==2017))>0){datMetrics[1,2] <- as.numeric(strftime(dat1[which(substr(dat1,1,4)==2017)][1], format = "%j"))}
    if(length(which(substr(dat1,1,4)==2018))>0){datMetrics[1,3] <- as.numeric(strftime(dat1[which(substr(dat1,1,4)==2018)][1], format = "%j"))}
    if(length(which(substr(dat1,1,4)==2019))>0){datMetrics[1,4] <- as.numeric(strftime(dat1[which(substr(dat1,1,4)==2019)][1], format = "%j"))}
    if(length(which(substr(dat1,1,4)==2020))>0){datMetrics[1,5] <- as.numeric(strftime(dat1[which(substr(dat1,1,4)==2020)][1], format = "%j"))}
    if(length(which(substr(dat1,1,4)==2021))>0){datMetrics[1,6] <- as.numeric(strftime(dat1[which(substr(dat1,1,4)==2021)][1], format = "%j"))}
  
    if(length(which(substr(dat2,1,4)==2016))>0){datMetrics[2,1] <- as.numeric(strftime(dat2[which(substr(dat2,1,4)==2016)][1], format = "%j"))}
    if(length(which(substr(dat2,1,4)==2017))>0){datMetrics[2,2] <- as.numeric(strftime(dat2[which(substr(dat2,1,4)==2017)][1], format = "%j"))}
    if(length(which(substr(dat2,1,4)==2018))>0){datMetrics[2,3] <- as.numeric(strftime(dat2[which(substr(dat2,1,4)==2018)][1], format = "%j"))}
    if(length(which(substr(dat2,1,4)==2019))>0){datMetrics[2,4] <- as.numeric(strftime(dat2[which(substr(dat2,1,4)==2019)][1], format = "%j"))}
    if(length(which(substr(dat2,1,4)==2020))>0){datMetrics[2,5] <- as.numeric(strftime(dat2[which(substr(dat2,1,4)==2020)][1], format = "%j"))}
    if(length(which(substr(dat2,1,4)==2021))>0){datMetrics[2,6] <- as.numeric(strftime(dat2[which(substr(dat2,1,4)==2021)][1], format = "%j"))}
  
    
    # Get Planet data
    # Point Shape file
    geog_crs = CRS("+proj=longlat +datum=WGS84")
    site <- data.frame(1,meta$phenocam_site$lon,meta$phenocam_site$lat)
    colnames(site) <- c('id','xcrd','ycrd')
    xy   <- site[,c(2,3)]
    bb   <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=geog_crs)
    
    utm_crs = rast@crs
    bb   <- spTransform(bb,utm_crs)
    
    site <- data.frame(1,bb@bbox[1],(bb@bbox[2]+30))
    colnames(site) <- c('id','xcrd','ycrd')
    xy   <- site[,c(2,3)]
    bb1  <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=utm_crs)
    
    # Get values
    for(j in 1:8){
      rast <- Rast[[j]]
      
      pheP1 <- extract(rast,bb)
      pheP2 <- unlist(extract(rast,bb,buffer=15))
      pheP3 <- extract(rast,bb1)
      pheP4 <- unlist(extract(rast,bb1,buffer=15))  
      
      if(j<=4){
        datMetrics[ 3,j] <- pheP1
        datMetrics[ 4,j] <- round(mean(pheP2,na.rm=T))
        datMetrics[ 5,j] <- round(median(pheP2,na.rm=T))
        datMetrics[ 6,j] <- pheP3
        datMetrics[ 7,j] <- round(mean(pheP4,na.rm=T))
        datMetrics[ 8,j] <- round(median(pheP4,na.rm=T))
      }else{
        datMetrics[ 9,(j-4)] <- pheP1
        datMetrics[10,(j-4)] <- round(mean(pheP2,na.rm=T))
        datMetrics[11,(j-4)] <- round(median(pheP2,na.rm=T))
        datMetrics[12,(j-4)] <- pheP3
        datMetrics[13,(j-4)] <- round(mean(pheP4,na.rm=T))
        datMetrics[14,(j-4)] <- round(median(pheP4,na.rm=T))
      }
    }
   
    rownames(datMetrics) <- c('PhenoCam_gup',
                              'PhenoCam_gdn',
                              'Planet_gup_center_single',
                              'Planet_gup_center_30m_mean',
                              'Planet_gup_center_30m_median',
                              'Planet_gup_offset_single',
                              'Planet_gup_offset_30m_mean',
                              'Planet_gup_offset_30m_median',
                              'Planet_gdn_center_single',
                              'Planet_gdn_center_30m_mean',
                              'Planet_gdn_center_30m_median',
                              'Planet_gdn_offset_single',
                              'Planet_gdn_offset_30m_mean',
                              'Planet_gdn_offset_30m_median')
    
    
    
    
    
    
    
    
    
    if(sum(!is.na(datMetrics[1:8]))>0 & sum(!is.na(datMetrics[9:56]))>0){
      save(afcd,vgty,meta, datMetrics,
           file=paste0('/projectnb/modislc/users/mkmoon/Planet/data_paper/data/comp_phenocam/V3_006/pc_',siteStr,'_',siteWPC[i],'.rda'))  
    } 
  
    print(i)
  }
  

}



# ########################################
# #
# ########################################
# files <- list.files('/projectnb/modislc/users/mkmoon/Planet/data_paper/data/comp_phenocam/V3',pattern=glob2rx('*.rda'),full.names=T)
# 
# dat <- matrix(NA,length(files),2)
# for(i in 1:length(files)){
#   load(files[i])
# 
#   dat[i,] <- c(afcd,vgty)
#   if(vgty=='WSA') print(i)
# }
# datMeta <- as.data.frame(dat)
# datAF   <- unique(dat[,1])
# unique(dat[,2])
# datAF[63]
# 
# par(mfrow=c(2,3))
# for(j in 3:8){
#   dat <- c()
#   for(i in 1:length(files)){
#     load(files[i])
# 
#     temp <- cbind(datMetrics[1,],datMetrics[j,])
#     dat <- rbind(dat,temp)
#   }
#   plot(dat,xlim=c(-100,500),ylim=c(-100,500))
#   abline(0,1)
# 
#   x1 <- dat[,1]
#   y1 <- dat[,2]
# 
#   reg <- formatC(round(cor(x1,y1,use='na.or.complete'),2), digits=2,format="fg", flag="#")
#   rmse <- round(sqrt(mean((x1-y1)^2,na.rm=T)))
#   bias <- round(mean(x1-y1,na.rm=T),1)
# 
#   print(rmse)
# }
# 
# par(mfrow=c(2,3))
# for(j in 9:14){
#   dat <- c()
#   for(i in 1:length(files)){
#     load(files[i])
# 
#     temp <- cbind(datMetrics[2,],datMetrics[j,])
#     dat <- rbind(dat,temp)
#   }
#   plot(dat,xlim=c(-100,500),ylim=c(-100,500))
#   abline(0,1)
# 
#   x1 <- dat[,1]
#   y1 <- dat[,2]
# 
#   reg <- formatC(round(cor(x1,y1,use='na.or.complete'),2), digits=2,format="fg", flag="#")
#   rmse <- round(sqrt(mean((x1-y1)^2,na.rm=T)))
#   bias <- round(mean(x1-y1,na.rm=T),1)
# 
#   print(rmse)
# }
# 
# 
# 
# ##
# setwd('/projectnb/modislc/users/mkmoon/Planet/data_paper/figure/')
# png(filename='1to1_pc_.png',width=11,height=5.5,unit='in',res=600)
# 
# mycol <- brewer.pal(7,'Set1')
# mycol <- c(mycol[1:5],mycol[7])
# mycol <- brewer.pal(4,'Greens')
# 
# 
# for(j in 1:2){
#   dat <- c()
#   if(j==1){vv <- c(1,6)
#   }else{vv <- c(2,12)}
# 
#   for(i in 1:length(files)){
#     load(files[i])
#     stcd <- which(datAF==afcd)
# 
#     if(vgty=='CRO'){       cc=1
#     }else if(vgty=='DBF'|vgty=='MF'){  cc=2
#     }else if(vgty=='ENF'){             cc=3
#     }else if(vgty=='GRA'){             cc=4
#     }else if(vgty=='OSH'|vgty=='CSH'){ cc=5
#     }else if(vgty=='WET'){             cc=6
#     }else if(vgty=='SAV'|vgty=='WSA'){ cc=7
#     }else if(vgty=='CVM'){             cc=8
#     }else if(vgty=='EBF'){             cc=9}
# 
# 
#     temp <- c(datMetrics[vv[1],],datMetrics[vv[2],],cc,stcd)
#     dat <- rbind(dat,temp)
#   }
# 
#   par(mfrow=c(3,3),oma=c(1,1,1,1),mar=c(4,4,1,1),mgp=c(2.5,1,0))
#   for(i in 1:9){
#     for(yy in 1:4){
#       if(yy==1){
#         plot(dat[dat[,9]==i,c(1,5)],xlim=c(-200,600),ylim=c(-200,600),pch=21,bg=mycol[yy])
#         text(dat[dat[,9]==i,1],dat[dat[,9]==i,5],dat[dat[,9]==i,10],pos=4)
#         abline(0,1)
#       }else{
#         points(dat[dat[,9]==i,c(yy,(4+yy))],pch=21,bg=mycol[yy])
#         text(dat[dat[,9]==i,yy],dat[dat[,9]==i,yy+4],dat[dat[,9]==i,10],pos=4)
#       }
#     }
#   }
# 
#   x1 <- dat[,1]
#   y1 <- dat[,2]
# 
#   reg <- formatC(round(cor(x1,y1,use='na.or.complete'),2), digits=2,format="fg", flag="#")
#   rmse <- round(sqrt(mean((x1-y1)^2,na.rm=T)))
#   bias <- round(mean(x1-y1,na.rm=T),1)
# }


