#Load required libraries
library(raster)
library(rgdal)
library(gdalUtils)
library(rgeos)
library(maptools)
library(rasterVis)
library(ncdf4)

ss <- 25


########################################
inBase <- '/projectnb/planet/PLSP/Product_010/'
outBase <- '/projectnb/planet/PLSP/Product_netCDF/'

sites <- list.dirs(inBase,recursive=F,full.names=F)
site  <- sites[ss]
print(site)

########################################
# Get product layers info
productTable <- read.csv('/usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/PLSP_Layers.csv',header=T,stringsAsFactors = F)
# Get a base image to pull raster info from     
baseImage <- raster(paste0('/projectnb/planet/PLSP/Img_cliped/',site,'/base_image.tif'))

########################################
# Get extent, and then define pixel centers in the x and y direction
ext = extent(baseImage)
res = res(baseImage)[1]
x = seq(ext[1]+res/2,ext[2]-res/2, res)
y = seq(ext[3]+res/2,ext[4]-res/2, res)
#Define dimensions for netCDF file
dimx = ncdim_def(name = 'x', longname = 'x coordinate', units='m', vals = as.double(x))
dimy = ncdim_def(name = 'y', longname = 'y coordinate', units='m', vals = rev(as.double(y)))


########################################
yy =1
  
  # Load files
  files <- list.files(paste0(inBase,site,'/',(2016+yy)),pattern=glob2rx(paste('*.tif',sep='')),full.names=T)
  files <- files[1]
  
  outFold <- paste0(outBase, site, '/')
  if(!dir.exists(outFold)) {dir.create(outFold)}
  outFile <- paste0(outFold,'test_PLSP_',(2016+yy),'.nc')
        
  # Loop through all the layers, and create a variable for each 
  results      <-vector("list", 1+1) 
  results[[1]] <- ncvar_def("crs","",list(),prec="char")
  
  for(i in 1){
      lyr <- productTable[i,]   # Pull the info for this layer from the productTable 
          
      if (lyr$data_type == 'Int16') {precision <- "short"}  #All are int16, so this isn't necessary
          
      # Create the variable, add to the list. Define the short_name, units, fill_value, long_name, and precision from the product table 
      results[[i+1]] <- ncvar_def(lyr$short_name, lyr$units, list(dimx,dimy), lyr$fill_value, lyr$long_name, prec=precision, compression=2)  
  }
        
  # Now create the netCDF file with the defined variables
  if (file.exists(outFile)) {file.remove(outFile)}
  ncout <- nc_create(outFile,results,force_v4=T)
        
  # Now loop through the layers again, this time actually 
  # writing the image data to the file
  for (i in 1) {
      lyr <- productTable[i,]
          
      # Open the data          
      mat <- matrix(values(raster(files[i],varname=lyr$short_name)),length(x),length(y))
      # Put the image into the file
      ncvar_put(ncout,results[[i+1]], mat) 
      
      # Fill in the attributes for the layer from the product table 
      ncatt_put(ncout,lyr$short_name,"scale",lyr$scale)
      ncatt_put(ncout,lyr$short_name,"offset",lyr$offset)
      ncatt_put(ncout,lyr$short_name,"data_type",lyr$data_type)
      ncatt_put(ncout,lyr$short_name,"valid_min",lyr$valid_min)
      ncatt_put(ncout,lyr$short_name,"valid_max",lyr$valid_max)
      
      print(paste((2016+yy),';',i))
  }
  
  ## Write the projection info for the transverse_mercator variable
  # Get projection in wkt format
  wkt <- showWKT(projection(baseImage), morphToESRI = FALSE)  
  # Need to pull the central meridian from the wkt
  spt <- unlist(strsplit(gsub(']','',wkt),','))
  central_meridian    <- as.numeric(spt[which(spt == "PARAMETER[\"central_meridian\"")+1])
  
  # Fill in the info.
  ncatt_put(ncout,"crs","grid_mapping_name","transverse_mercator")
  ncatt_put(ncout,"crs","scale_factor_at_central_meridian",0.9996)
  ncatt_put(ncout,"crs","longitude_of_central_meridian",central_meridian)
  ncatt_put(ncout,"crs","false_easting",5e+05)
  ncatt_put(ncout,"crs","false_northing",0)
  
  ## Close the file
  nc_close(ncout)



###########
n1 <- nc_open('/projectnb/planet/PLSP/Product_netCDF/Harvard_Forest_EMS_Tower__HFR1_/test_PLSP_2017.nc')
n2 <- raster('/projectnb/planet/PLSP/Product_netCDF/Harvard_Forest_EMS_Tower__HFR1_/test_PLSP_2017.nc',varname='NumCycles')
n3 <- raster('/projectnb/planet/PLSP/Product_netCDF/Harvard_Forest_EMS_Tower__HFR1_/PLSP_2017.nc',varname='NumCycles')
n4 <- nc_open('/projectnb/modislc/data/climate/daymet/v4/daymet_v4_daily_na_dayl_2016.nc')
n4 <- raster('/projectnb/modislc/data/climate/daymet/v4/daymet_v4_daily_na_dayl_2016.nc',varname='dayl')
compareRaster(n2,n3,orig=T)

nn <- nc_open('/projectnb/planet/PLSP/Product_netCDF_fillVal/Abby_Road_NEON/PLSP_2017.nc')
nn <- raster('/projectnb/planet/PLSP/Product_netCDF_fillVal/Abby_Road_NEON/PLSP_2017.nc',varname='50PCGI')
sum(is.na(values(nn)))


files <- list.files('/projectnb/planet/PLSP/Product_netCDF/',pattern = glob2rx('*.nc'),recursive=T,full.names = T)
for(i in 1:length(files)){
  n1 <- raster(files[i],varname='QA')
  n2 <- raster(files[i],varname='QA_2')
  if(sum(is.na(values(n1)))!=0|sum(is.na(values(n2)))!=0){
    paste('NA here:',i)
  }else{
    if(i%%10==0) print(i)
  }
}


