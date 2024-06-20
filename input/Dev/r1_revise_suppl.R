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

ss <- as.numeric(args[3])
# ss <- 24

vari <- c('50PCGI','50PCGD')


########################################
Rast <- vector('list',10); ind <- 1
for(vv in 1:2){
  for(yy in 2017:2021){
    imgs <- list.files('/projectnb/planet/PLSP/Product_010',pattern=glob2rx(paste0('*',yy,'_',vari[vv],'.tif')),full.names=T,recursive=T)  
    
    Rast[[ind]] <- raster(imgs[ss])    
    
    print(ind)
    ind  <- ind + 1 
  }
}

siteStr <- unlist(strsplit(imgs[ss],'/'))[6]
print(siteStr)


siteFull <- read.csv('/projectnb/modislc/users/mkmoon/Planet/data_paper/data/site_full_w_pc_hls_PCrevised_ext.csv')

# afcd <- siteFull$Site_amerifulx[which(siteFull$Site==substr(siteStr,5,80))]
# vgty <- siteFull$Veg[which(siteFull$Site==substr(siteStr,5,80))]
afcd <- siteFull$Site_amerifulx[which(siteFull$Site==siteStr)]
vgty <- siteFull$Veg[which(siteFull$Site==siteStr)]

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
  
  if(length(intPt)>0){
    siteWPC <- c(siteWPC,meta$phenocam_site$sitename)
    
    siteMeta <- c(afcd,
                  meta$phenocam_site$sitename,
                  meta$phenocam_site$long_name,
                  meta$phenocam_site$camera_orientation,
                  meta$phenocam_site$lat,
                  meta$phenocam_site$lon,
                  meta$phenocam_site$elevation,
                  meta$phenocam_site$primary_veg_type,
                  meta$phenocam_site$secondary_veg_type,
                  meta$phenocam_site$date_start,
                  meta$phenocam_site$date_end,
                  meta$phenocam_site$landcover_igbp,
                  meta$phenocam_site$site_acknowledgements)
    save(siteMeta,file=paste0('/projectnb/modislc/users/mkmoon/Planet/data_paper/data/comp_phenocam/site_pc_list_2/pc_',afcd,'_',meta$phenocam_site$sitename,'.rda'))
  } 
  
  if(i%%100==0) print(i)
}


print(length(siteWPC))





########################################
# PhenoCam meta
library(dplyr)
library(tibble)

files <- list.files('/projectnb/modislc/users/mkmoon/Planet/data_paper/data/comp_phenocam/site_pc_list_2',full.names = T)

df <- as.data.frame(matrix(NA,length(files),12))
for(i in 1:length(files)){
  load(files[i])
  df[i,1:length(siteMeta)] <- siteMeta
}
length(unique(df$V1))
length(unique(df$V2))

dat <- df %>%
  rownames_to_column() %>%
  filter(!duplicated(V2))
# 
# write.csv(dat,file='/projectnb/modislc/users/mkmoon/Planet/data_paper/data/comp_phenocam/pc_site_meta_ext.csv')

