#Load required libraries
library(raster)
library(rgdal)
library(gdalUtils)
library(rgeos)
library(maptools)
library(rasterVis)
library(ncdf4)

args <- commandArgs()
print(args)

ss <- as.numeric(args[3])
# ss <- 1

inBase <- '/projectnb/planet/PLSP/Product_010/'
outBase <- '/projectnb/planet/PLSP/Product_netCDF/'

sites <- list.dirs(inBase,recursive=F,full.names=F)
site  <- sites[ss]
print(site)

lyrs <- read.csv('/usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/PLSP_Layers.csv',header=T,stringsAsFactors = F)


#Get a base imageto pull raster info from     
baseImage <- raster(paste0('/projectnb/planet/PLSP/Img_cliped/',site,'/base_image.tif'))

#Get extent, and then define pixel centers in the x and y direction
ext = extent(baseImage)
res = res(baseImage)[1]
x = seq(ext[1]+res/2,ext[2]-res/2, res)
y = seq(ext[3]+res/2,ext[4]-res/2, res)

#Define dimensions for netCDF file
dimx = ncdim_def(name = 'x', longname = 'x coordinate', units='m', vals = as.double(x))
dimy = ncdim_def(name = 'y', longname = 'y coordinate', units='m', vals = rev(as.double(y)))




for(yy in 1:5){
  
  dimT = ncdim_def(name = 't', longname = 'year', units='year', calendar="standard", vals=(2016+yy))
  
  files <- list.files(paste0(inBase,site,'/',(2016+yy)),pattern=glob2rx(paste('*.tif',sep='')),full.names=T)
  
  outFold <- paste0(outBase, site, '/')
  if(!dir.exists(outFold)) {dir.create(outFold)}
  
  outFile <- paste0(outFold,'PLSP_',yrs[yy],'.nc')
        
  #Loop through all the layers, and create a variable for each 
  results      <-vector("list", dim(lyrs)[1]+1) 
  results[[1]] <- ncvar_def("transverse_mercator","",list(),prec="char")
  
  for(i in 1:dim(lyrs)[1]){
      lyr <- lyrs[i,]   #Pull the info for this layer from the lyrs table
          
      if (lyr$data_type == 'Int16') {precision <- "short"}  #All are int16, so this isn't necessary
          
      #Create the variable, add to the list. Define the short_name, units, fill_value, long_name, and precision from the lyrs table
      results[[i+1]] <- ncvar_def(lyr$short_name, lyr$units, list(dimx,dimy), lyr$fill_value, lyr$long_name, prec=precision, compression=2)  
  }
        
  # Now create the netCDF file with the defined variables
  if (file.exists(outFile)) {file.remove(outFile)};
  ncout <- nc_create(outFile,results,force_v4=T)
      
        
  #Now loop through the layers again, this time actually 
  #writing the image data to the file
  for (i in 1:dim(lyrs)[1]) {
      lyr <- lyrs[i,]
          
      #Open the data. Current solution to "WRITE_BOTTOMUP" is 
      #to manually flip the image 
      mat <- matrix(values(raster(files[i],varname=lyr$short_name)),length(x),length(y))
      
      ncvar_put(ncout,results[[i+1]], mat)      #Now put the image into the file
      
      #Fill in the attributes for the layer from the lyrs table
      ncatt_put(ncout,lyr$short_name,"scale",lyr$scale)
      ncatt_put(ncout,lyr$short_name,"offset",lyr$offset)
      ncatt_put(ncout,lyr$short_name,"data_type",lyr$data_type)
      ncatt_put(ncout,lyr$short_name,"valid_min",lyr$valid_min)
      ncatt_put(ncout,lyr$short_name,"valid_max",lyr$valid_max)
      
      print(paste(yrs[yy],';',i))
  }
  
  #Now we need to write the projection info the the transverse_mercator variable
  
  wkt <- showWKT(projection(baseImage), morphToESRI = FALSE)  #Get projection in wkt format
  #Need to pull the central meridian from the wkt (Hopefull this will always work?)
  spt <- unlist(strsplit(gsub(']','',wkt),','))
  central_meridian    <- as.numeric(spt[which(spt == "PARAMETER[\"central_meridian\"")+1])
  
  #Fill in the info. Manaully wrote values that don't change across UTM zones (Is that correct?)
  ncatt_put(ncout,"transverse_mercator","long_name","CRS definition")
  ncatt_put(ncout,"transverse_mercator","grid_mapping_name","transverse_mercator")
  ncatt_put(ncout,"transverse_mercator","longitude_of_central_meridian",central_meridian)
  ncatt_put(ncout,"transverse_mercator","false_easting",5e+05)
  ncatt_put(ncout,"transverse_mercator","false_northing",0)
  ncatt_put(ncout,"transverse_mercator","latitude_of_projection_origin",0)
  ncatt_put(ncout,"transverse_mercator","scale_factor_at_central_meridian",0.9996)
  ncatt_put(ncout,"transverse_mercator","longitude_of_prime_meridian",0)
  ncatt_put(ncout,"transverse_mercator","semi_major_axis",6378137)
  ncatt_put(ncout,"transverse_mercator","inverse_flattening",298.257223563)
  ncatt_put(ncout,"transverse_mercator","GeoTransform",paste(ext[1],res,0,ext[4],0,-res))
  ncatt_put(ncout,"transverse_mercator","spatial_ref",gsub("\\", "", wkt, fixed=TRUE))
  
  #Define global attributes
  ncatt_put(ncout,0,"title","Land Surface Phenology from PlanetScope (PLSP)")
  ncatt_put(ncout,0,"product_version","v001")
  ncatt_put(ncout,0,"summary","A High Spatial Resolution Land Surface Phenology for AmeriFlux and NEON sites")
  ncatt_put(ncout,0,"software_repository","git@github.com:BU-LCSC/PLSP.git")

  ncatt_put(ncout,0,"creator_name","Land Cover & Surface Climate Group, Department of Earth & Environment, Boston University")
  ncatt_put(ncout,0,"creator_type","group")
  ncatt_put(ncout,0,"creator_email","mkmoon@bu.edu")
  ncatt_put(ncout,0,"creator_institution","Boston University")
  
  ncatt_put(ncout,0,"contributor_name", "Minkyu Moon, Mark A. Friedl")
  ncatt_put(ncout,0,"contributor_role", "Developer, Principal Investigator")
  ncatt_put(ncout,0,"acknowledgement","Developed with funding from NASA LCLUC Grant #80NSSC18K0334.")
  
  #Put additional attributes on coordinates
  ncatt_put(ncout,"x","axis","projection_x_coordinate")
  ncatt_put(ncout,"y","axis","projection_y_coordinate")
  ncatt_put(ncout,"t","axis","year")
  
  nc_close(ncout)
}

}

###########
n1 <- nc_open('/projectnb/planet/PLSP/Product_netCDF/Abby_Road_NEON/PLSP_2017.nc')

lyrs$short_name
sites <- list.dirs(inBase,recursive=F,full.names=F)
for(ss in 1:20){
  site  <- sites[ss]
  for(yy in 2017:2021){
    for(i in 1:24){
      ii <- sprintf('%02d',i)
      n2 <- raster(paste0('/projectnb/planet/PLSP/Product_netCDF/',site,'/PLSP_',yy,'.nc'),varname=lyrs$short_name[i])
      n3 <- raster(paste0('/projectnb/planet/PLSP/Product_010/',site,'/',yy,'/',ii,'_',yy,'_',lyrs$short_name[i],'.tif'))
      ifelse(compareRaster(n2,n3,res=T,orig=T),'',print(paste(site,yy)))
    }
  }
  print(site)
}




