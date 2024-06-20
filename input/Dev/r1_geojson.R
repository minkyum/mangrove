library(sp)
library(raster)
library(terra)
library(sf)

library(geojsonsf)


# #################################################
# imgBases <- list.files('/projectnb/planet/PLSP/Img_cliped',pattern=glob2rx('base*.tif'),full.names=T,recursive=T)
# siteFull <- read.csv('/projectnb/modislc/users/mkmoon/Planet/data_paper/data/site_full_w_pc_hls_PCrevised_ext.csv')
# 
# ss <- 1
# 
# for(ss in 1:104){
#   
#   imgBase <- raster(imgBases[ss])
#   
#   e <- imgBase@extent
#   p <- as(e, 'SpatialPolygons')
#   crs(p) <- crs(imgBase)
#   p <- as(p, 'SpatialPolygonsDataFrame')
#   
#   str1 <- unlist(strsplit(imgBases[ss],'/'))[6]
#   
#   str2 <- siteFull$Site_full[which(siteFull$Site==str1)]
#   afcd <- siteFull$Site_amerifulx[which(siteFull$Site==str1)]
#   str2 <- unlist(strsplit(str2,' '))
#   str2 <- paste(str2,collapse='_')
#   strF <- paste0(afcd,'__',str2)
#   
#   tf <- paste0('/projectnb/modislc/users/mkmoon/Planet/rawImage/geojson (copy)/',strF,'.geojson')
#   writeOGR(p,tf,"GeoJSON", driver="GeoJSON")
#   
# }
# 
# #############################
# # Toolik
# p <- st_read('/projectnb/modislc/users/mkmoon/cdsa/plsp/shp/cdsa_toolik.shp')
# p$id <- 'cdsa_toolik'
# crs_geog <- '+proj=longlat +datum=WGS84 +no_defs'
# p <- st_transform(p,crs(crs_geog))
# tf <- '/projectnb/modislc/users/mkmoon/cdsa/plsp/geojson/cdsa_toolik.geojson'
# st_write(p,tf,layer='cdsa_toolik',delete_dsn=T)



#############################
## For Andrew; US-CdM & US-Mpj
crs_geog <- '+proj=longlat +datum=WGS84 +no_defs'

# US-CdM
x <- c(-109.7471,37.5241)
p <- st_point(x)
p <- st_sfc(p,crs=crs_geog)
p <- st_transform(p,32612) # UTM 12N
pts <- list(matrix(c(p[[1]][1]-3000,p[[1]][2]-3000,
                     p[[1]][1]+3000,p[[1]][2]-3000,
                     p[[1]][1]+3000,p[[1]][2]+3000,
                     p[[1]][1]-3000,p[[1]][2]+3000,
                     p[[1]][1]-3000,p[[1]][2]-3000),ncol=2, byrow=TRUE))
p  <- st_polygon(pts)
p  <- st_sfc(p,crs=32612)
p  <- st_transform(p,4326)
p  <- st_sf(p)
gj <- sf_geojson(p)              


# US-Mpj
x <- c(-106.2377,34.4385)
p <- st_point(x)
p <- st_sfc(p,crs=crs_geog)
p <- st_transform(p,32613)
pts <- list(matrix(c(p[[1]][1]-3000,p[[1]][2]-3000,
                     p[[1]][1]+3000,p[[1]][2]-3000,
                     p[[1]][1]+3000,p[[1]][2]+3000,
                     p[[1]][1]-3000,p[[1]][2]+3000,
                     p[[1]][1]-3000,p[[1]][2]-3000),ncol=2, byrow=TRUE))
p  <- st_polygon(pts)
p  <- st_sfc(p,crs=32613)
p  <- st_transform(p,4326)
p  <- st_sf(p)
gj <- sf_geojson(p)   


c(-90.144333144514576, 45.850891001799639,
  -90.015466855485414, 45.850891001799639,
  -90.015466855485414, 45.760908998200364,
  -90.144333144514576, 45.760908998200364, 
  -90.144333144514576, 45.850891001799639)


#############################
## For Gavin; UIC
crs_geog <- '+proj=longlat +datum=WGS84 +no_defs'

# UIC
x <- c(-87.658846,41.870069)
p <- st_point(x)
p <- st_sfc(p,crs=crs_geog)
p <- st_transform(p,32616) # UTM 12N
pts <- list(matrix(c(p[[1]][1]-1500,p[[1]][2]-1500,
                     p[[1]][1]+1500,p[[1]][2]-1500,
                     p[[1]][1]+1500,p[[1]][2]+1500,
                     p[[1]][1]-1500,p[[1]][2]+1500,
                     p[[1]][1]-1500,p[[1]][2]-1500),ncol=2, byrow=TRUE))
p  <- st_polygon(pts)
p  <- st_sfc(p,crs=32616)
p  <- st_transform(p,4326)
p  <- st_sf(p)
gj <- sf_geojson(p,simplify=F)              

