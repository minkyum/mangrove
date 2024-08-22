library(sp)
library(raster)
library(terra)
library(sf)
library(rjson)



###############################
args <- commandArgs()
print(args)

numSite <- as.numeric(substr(args[3],1,3))
yToDo   <- as.numeric(substr(args[3],4,7))
# numSite <- 1; yToDo <- 4



########################################
params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/mangrove/input/PLCM_Parameters.json')
source(params$setup$rFunctions)

productTable <- read.csv(params$setup$productTable,header=T,stringsAsFactors = F)
phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr



########################################
## Get site name, image directory and coordinate
strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

ckPheDir <- paste0(params$setup$outDir,strSite,'/chunk_feat')
print(ckPheDir)

files <- list.files(path=ckPheDir,pattern=glob2rx('*.rda'),full.names=T)
print(length(files))



########################################
imgBase <- rast(paste0(params$setup$outDir,strSite,'/base_image.tif'))
numPix <- ncol(imgBase)*nrow(imgBase)
imgNum <- setValues(imgBase,1:numPix)
numChunks <- params$setup$numChunks
chunk <- numPix%/%numChunks



########################################
# Load features
Fmat <- matrix(NA,numPix,7) 
for(i in 1:numChunks){
  cc <- sprintf('%03d',i)
  cfile <- paste0(ckPheDir,'/chunk_feat_',cc,'.rda') 
  log <- try(load(cfile),silent=F)
  if (inherits(log, 'try-error')) next 
    
  if(i==numChunks){chunks <- c((chunk*(i-1)+1):numPix)
  }else{chunks <- c((chunk*(i-1)+1):(chunk*i))}
    
  chunkStart <- chunks[1];  chunkEnd <- chunks[length(chunks)]
    
  Fmat[chunkStart:chunkEnd,] <- f_mat[,(58*(yToDo-1)+c(42:48))]
  
 print(i)
}



########################################
pheDir  <- paste0(params$setup$outDir,strSite,'/output')
map_lct <- rast(paste0(pheDir,'/map_re_',2017+yToDo,'.tif'))

map_p1  <- setValues(imgBase,Fmat[,2])

# Plot
plot(map_p1)



#######################################
# Save
pheDir <- paste0(params$setup$outDir,strSite,'/output')
if (!dir.exists(pheDir)) {dir.create(pheDir)}

writeRaster(pred_map,filename=paste0(pheDir,'/map_re_1_',2017+yToDo,'.tif'),overwrite=TRUE)
save(rf,file=paste0(pheDir,'/rf_pc_',2017+yToDo,'.rda'))



