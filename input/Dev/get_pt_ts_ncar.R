library(raster)
library(rgdal)
library(gdalUtils)

library(rjson)
library(geojsonR)


########################################

args <- commandArgs()
print(args)

ss <- as.numeric(args[3])
# ss <- 63
print(ss)

metSites <- read.csv('/projectnb/modislc/users/mkmoon/Planet/data/ncar/site_list_neon_ncar.csv')

########################################
params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/PLSP_Parameters.json')
source(params$setup$rFunctions)


########################################
## Get site name and image directory
geojsonDir <- params$setup$geojsonDir

strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[metSites$numSite[ss]]
print(strSite)

ckDir <- paste0(params$setup$outDir,strSite,'/chunk')
print(ckDir)



########################################
## 
imgBase <- raster(paste0(params$setup$outDir,strSite,'/base_image.tif'))
numPix <- length(imgBase)
imgNum <- setValues(imgBase,1:numPix)
numChunks <- params$setup$numChunks
chunk <- numPix%/%numChunks


# Get site name, image directory and coordinate
# Point Shape file
geog_crs = CRS("+proj=longlat +datum=WGS84")
site <- data.frame(1,metSites$field_longitude[ss],metSites$field_latitude[ss])
colnames(site) <- c('id','xcrd','ycrd')
xy   <- site[,c(2,3)]
bb   <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=geog_crs)

utm_crs = imgBase@crs
bb   <- spTransform(bb,utm_crs)

site <- data.frame(1,bb@bbox[1],(bb@bbox[2]+30))
colnames(site) <- c('id','xcrd','ycrd')
xy   <- site[,c(2,3)]
bb1  <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=utm_crs)

pixNums <- unlist(extract(imgNum,bb1,buffer=6))

print(length(pixNums))

##
pixNum  <- pixNums[1]
ckNum <- sprintf('%03d',(pixNum%/%chunk+1))
file <- list.files(path=ckDir,pattern=glob2rx(paste0('*',ckNum,'.rda')),full.names=T)

load(file)

valVI0 <- matrix(NA,length(dates),length(pixNums)) # ndvi
valVI1 <- matrix(NA,length(dates),length(pixNums)) # evi2

##
for(i in 1:length(pixNums)){
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
  
  valVI0[,i] <-     (nir - red) / (nir +     red    ) # ndvi
  valVI1[,i] <- 2.5*(nir - red) / (nir + 2.4*red + 1) # evi2
  
  print(i)
}

valVI01 <- data.frame(dates,valVI0,valVI1)
valvi   <-na.omit(valVI01)

colnames(valvi) <- c('date',
                     paste0('ndvi_pixel_',1:length(pixNums)),
                     paste0('evi2_pixel_',1:length(pixNums)))
                     

##
setwd('/projectnb/modislc/users/mkmoon/Planet/data/ncar/')
write.csv(valvi,file=paste0('planet_vi_',metSites$field_site_id[ss],'_',metSites$field_site_name[ss],'.csv'))



      