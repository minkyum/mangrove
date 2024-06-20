library(raster)
library(rgdal)
library(gdalUtils)
library(rgeos)

library(rjson)
library(geojsonR)


########################################
imgBase <- list.files('/projectnb/planet/PLSP/Img_cliped',pattern='base_i*',full.names=T,recursive=T)

# Get site name, image directory and coordinate
sitePlanet <- matrix(NA,length(imgBase),3)
for(i in 1:length(imgBase)){
  sitePlanet[i,1] <- unlist(strsplit(imgBase[i],'/'))[6]

  rast <- raster(imgBase[i])

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


  sitePlanet[i,2] <- bb@coords[1]
  sitePlanet[i,3] <- bb@coords[2]
}
colnames(sitePlanet) <- c('Site_Planet','Long','Lat')
sitePlanet <- as.data.frame(sitePlanet)


# PhenoCam site
pcMeta <- list.files('/projectnb/planet/PhenoCam/PhenoCam_V3/data_record_1',pattern='*_meta.json',full.names=T)

sitePC <- matrix(NA,length(pcMeta),4)
for(i in 1:length(pcMeta)){
  meta <- fromJSON(file=pcMeta[i])
  
  sitePC[i,1] <- meta$phenocam_site$sitename
  if(length(meta$phenocam_site$flux_sitenames)>0){
    sitePC[i,2] <- meta$phenocam_site$flux_sitenames 
  }
  sitePC[i,3] <- meta$phenocam_site$lon
  sitePC[i,4] <- meta$phenocam_site$lat
}
colnames(sitePC) <- c('Site_PC','Site_AF','Long','Lat')
sitePC <- as.data.frame(sitePC)

# List site having PhenoCam site
siteWPC <- c()
for(i in 1:dim(sitePlanet)[1]){
  rast <- raster(imgBase[i])
  
  for(j in 1:dim(sitePC)[1]){
    # Point Shape file
    geog_crs = CRS("+proj=longlat +datum=WGS84")
    site <- data.frame(1,as.numeric(sitePC$Long[j]),as.numeric(sitePC$Lat[j]))
    colnames(site) <- c('id','xcrd','ycrd')
    xy   <- site[,c(2,3)]
    bb   <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=geog_crs)
    
    utm_crs = rast@crs
    bb   <- spTransform(bb,utm_crs)
    
    log <- try({
      intPt <- intersect(bb,rast)  
      
      if(length(intPt)>0){
        temp <- c(sitePlanet[i,],sitePC[j,])
        siteWPC <- rbind(siteWPC,temp)
      } 
    },silent=T)
  }
  print(i)
}
siteFull <- as.data.frame(siteWPC)



########################################
# HLS tile
shpHLS <- shapefile('/projectnb/modislc/projects/landsat_sentinel/shapefiles/sentinel2_tiles_world/sentinel2_tiles_world/sentinel2_tiles_north_america_Albers.shp')

tileHLS <- matrix(NA,dim(siteFull)[1],3)
for(i in 1:dim(siteFull)[1]){
  # Point Shape file
  geog_crs = CRS("+proj=longlat +datum=WGS84")
  site <- data.frame(1,siteFull$Long[i],siteFull$Lat[i])
  colnames(site) <- c('id','xcrd','ycrd')
  xy   <- site[,c(2,3)]
  bb   <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=geog_crs)
  bb   <- spTransform(bb,shpHLS@proj4string)

  aa <- intersect(shpHLS,bb)

  if(length(aa$Name)>0){
    tileHLS[i,1] <- length(aa$Name)
    tileHLS[i,2:(length(aa$Name)+1)] <- aa$Name
  }else{
    tileHLS[i,1] <- 0
  }
}
# 
# 
# siteFull <- data.frame(siteFull,tileHLS)
# colnames(siteFull) <- c('Site_amerifulx','Site_full','Lat','Long','Elev','Veg','Clim','MAT','MAP','Site Start','Site End','BASE Start','BASE End','','Site','Lat','Long',
#                         '# of PhenoCam','PC1','PC2','PC3','PC4','# of HLS','HLS1','HLS2')
# 
# 
# # write.csv(siteFull,'/projectnb/modislc/users/mkmoon/Planet/data_paper/data/site_full_w_pc_hls.csv')




########################################
# Rename directory
siteFull <- read.csv('/projectnb/modislc/users/mkmoon/Planet/data_paper/data/site_full_w_pc_hls_PCrevised_ext.csv')

dirAll <- list.dirs('/projectnb/planet/PLSP/Product_netCDF_fillVal',full.names=T,recursive=F)

for(i in 1:dim(siteFull)[1]){
  dirOld <- dirAll[i]

  str1 <- unlist(strsplit(dirOld,'/'))[6]
  str2 <- siteFull$Site_full[which(siteFull$Site==str1)]
  afcd <- siteFull$Site_amerifulx[which(siteFull$Site==str1)]
  str2 <- unlist(strsplit(str2,' '))
  # if(str2[1]=='NEON'){
  #   str2 <- str2[2:length(str2)]
  # }
  str2 <- paste(str2,collapse='_')
  strF <- paste0(afcd,'__',str2)

  dirNew <- paste0('/projectnb/planet/PLSP/Product_netCDF_fillVal/',strF)

  system(paste0('mv ',dirOld,' ',dirNew))

  # print(paste(str1,'     ',strF))
}



########################################
siteFull <- read.csv('/projectnb/modislc/users/mkmoon/Planet/data_paper/data/site_full_w_pc_hls_PCrevised_ext.csv')

shpNA <- shapefile('/projectnb/planet/PLSP/Shape_file/map_na.shp')

plot(shpNA)

# siteYpc <- siteFull[siteFull$X..of.PhenoCam>0,]
# siteNpc <- siteFull[siteFull$X..of.PhenoCam==0,]

# Point Shape file
geog_crs = CRS("+proj=longlat +datum=WGS84")
# site <- data.frame(siteYpc$Site_amerifulx,siteYpc$Long,siteYpc$Lat)
site <- data.frame(siteFull$Site_amerifulx,siteFull$Long,siteFull$Lat)
colnames(site) <- c('id','xcrd','ycrd')
xy   <- site[,c(2,3)]
bb   <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=geog_crs)
bb   <- spTransform(bb,shpNA@proj4string)

setwd('/projectnb/modislc/users/mkmoon/Planet/data_paper/shp/')
writeOGR(bb,".","site_point_ext",driver="ESRI Shapefile",overwrite=T)

# site <- data.frame(siteNpc$Site_amerifulx,siteNpc$Long,siteNpc$Lat)
# colnames(site) <- c('id','xcrd','ycrd')
# xy   <- site[,c(2,3)]
# bb   <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=geog_crs)
# bb   <- spTransform(bb,shpNA@proj4string)
# 
# setwd('/projectnb/modislc/users/mkmoon/Planet/data_paper/shp/')
# writeOGR(bb,".","site_point_n_pc1",driver="ESRI Shapefile",overwrite=T)

plot(bb,add=T)


########################################
shpNA <- shapefile('/projectnb/planet/PLSP/Shape_file/map_na_re.shp')

tile_list <- c('h08v07','h09v07','h10v07',
               'h07v03','h08v03','h09v03',
               'h09v02','h10v02','h11v02','h12v02','h13v02','h14v02',
               'h08v04','h08v05','h08v06',
               'h09v04','h09v05','h09v06','h10v03','h10v04',
               'h10v05','h10v06','h11v03','h11v04','h11v05',
               'h12v03','h12v04','h12v05','h13v03','h13v04','h14v03','h14v04',
               'h03v06','h03v07','h11v07')
tile_list <- c('h07v05','h07v06',
               'h12v01','h13v01','h14v01','h15v01',
               'h15v02','h11v08')
tile_list <- c('h16v01','h15v02','h11v06')
tile_list <- c('h16v00','h17v00','h18v00','h19v00')
tile_list <- c('h16v00','h17v00','h18v00','h19v00')
tile_list <- c('h16v02','h17v01','h17v02')


setwd('/projectnb/modislc/users/mkmoon/Planet/data_paper/data/lct/')

path <- '/projectnb/modislc/users/mkmoon/mcd12q1/e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/2019.01.01/'
for(i in 1:length(tile_list)){
  file <- list.files(path,pattern=glob2rx(paste0('MCD12Q1.A2019001.',tile_list[i],'*.hdf')),full.names=T)
  
  sds <- get_subdatasets(file)
  lct <- raster(sds[1])
  
  pr3 <- projectExtent(lct,crs(shpNA))
  res(pr3) <- 500
  lct <- projectRaster(lct,pr3,method='ngb')
  
  writeRaster(lct,filename=paste0('lct_leae_',tile_list[i],'.tif'), format="GTiff", overwrite=TRUE)
  
  print(i)
}






########################################
# Rename directory for "SiteMaps"
siteFull <- read.csv('/projectnb/modislc/users/mkmoon/Planet/data_paper/data/site_full_w_pc_hls_PCrevised_ext.csv')

# dirAll <- list.dirs('/projectnb/modislc/users/mkmoon/Planet/data_paper/SiteMaps',full.names=T,recursive=F)
dirAll <- list.files('/projectnb/modislc/users/mkmoon/Planet/data_paper/SiteMaps',pattern=glob2rx('*.html'),full.names=T,recursive=F)

for(i in 2:dim(siteFull)[1]){
  dirOld <- dirAll[i]
  
  str1 <- unlist(strsplit(dirOld,'/'))[9]
  str1 <- unlist(strsplit(str1,'[.]'))[1]
  str1 <- unlist(strsplit(str1,' '))
  str1 <- paste(str1,collapse='_')
  if(i==92) str1 <- 'Univ._of_Mich._Biological_Station'
  if(i==39) str1<- 'Lyndon_B._Johnson_National_Grassland_NEON'
  try({
    str1 <- gsub('[(]','_',str1)
    str1 <- gsub('[)]','_',str1)  
    str1 <- gsub('&','and',str1)  
    str1 <- gsub("'", "_", str1)
    str1 <- gsub("#", "_", str1)
    str1 <- gsub(",", "_", str1)
  })
  str2 <- siteFull$Site_full[which(siteFull$Site==str1)]
  afcd <- siteFull$Site_amerifulx[which(siteFull$Site==str1)]
  str2 <- unlist(strsplit(str2,' '))
  str2 <- paste(str2,collapse='_')
  strF <- paste0(afcd,'__',str2)
  
  dirNew <- paste0('/projectnb/modislc/users/mkmoon/Planet/data_paper/SiteMaps/',strF,'.html')
  
  print(dirOld)
  print(dirNew)
  
  file.rename(dirOld,dirNew)
  
  # print(paste(str1,'     ',strF))
}

# NGEE
# Lost Creek
# Konza Prairie LTER (KNZ).geojson_files
# Uti
