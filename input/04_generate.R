library(sp)
library(raster)
library(terra)
library(sf)
library(rjson)


###############################
args <- commandArgs()
print(args)

numSite <- as.numeric(args[3])
# numSite <- 1


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



#Get all images to process
########################################
imgBase <- raster(paste0(params$setup$outDir,strSite,'/base_image.tif'))
numPix <- length(imgBase)
numChunks <- params$setup$numChunks
chunk <- numPix%/%numChunks



# Save
pheDir <- paste0(params$setup$outDir,strSite,'/feat')
if (!dir.exists(pheDir)) {dir.create(pheDir)}


for(yToDo in 1:length(phenYrs)){
  
  l01 <- matrix(NA,numPix,1);l02 <- matrix(NA,numPix,1);l03 <- matrix(NA,numPix,1);l04 <- matrix(NA,numPix,1)
  l05 <- matrix(NA,numPix,1);l06 <- matrix(NA,numPix,1);l07 <- matrix(NA,numPix,1);l08 <- matrix(NA,numPix,1)
  l09 <- matrix(NA,numPix,1);l10 <- matrix(NA,numPix,1);l11 <- matrix(NA,numPix,1);l12 <- matrix(NA,numPix,1)
  
  for(i in 1:numChunks){
    cc <- sprintf('%03d',i)
    cfile <- paste0(ckPheDir,'/chunk_feat_',cc,'.rda') 
    log <- try(load(cfile),silent=F)
    if (inherits(log, 'try-error')) next 
    
    if(i==numChunks){chunks <- c((chunk*(i-1)+1):numPix)
    }else{chunks <- c((chunk*(i-1)+1):(chunk*i))}
    
    chunkStart <- chunks[1];  chunkEnd <- chunks[length(chunks)]
    
    l01[chunkStart:chunkEnd,] <- f_mat[, (58*(yToDo-1)+41)];l02[chunkStart:chunkEnd,] <- f_mat[, (58*(yToDo-1)+42)];l03[chunkStart:chunkEnd,] <- f_mat[, (58*(yToDo-1)+43)]
    l04[chunkStart:chunkEnd,] <- f_mat[, (58*(yToDo-1)+44)];l05[chunkStart:chunkEnd,] <- f_mat[, (58*(yToDo-1)+45)];l06[chunkStart:chunkEnd,] <- f_mat[, (58*(yToDo-1)+46)]
    l07[chunkStart:chunkEnd,] <- f_mat[, (58*(yToDo-1)+48)];l08[chunkStart:chunkEnd,] <- f_mat[, (58*(yToDo-1)+37)];l09[chunkStart:chunkEnd,] <- f_mat[, (58*(yToDo-1)+37)]
    l10[chunkStart:chunkEnd,] <- f_mat[, (58*(yToDo-1)+40)];l11[chunkStart:chunkEnd,] <- f_mat[, (58*(yToDo-1)+39)];l12[chunkStart:chunkEnd,] <- f_mat[, (58*(yToDo-1)+58)]
  
  }
  
  r01 <- setValues(imgBase,l01); r02 <- setValues(imgBase,l02); r03 <- setValues(imgBase,l03); r04 <- setValues(imgBase,l04)
  r05 <- setValues(imgBase,l05); r06 <- setValues(imgBase,l06); r07 <- setValues(imgBase,l07); r08 <- setValues(imgBase,l08)
  r09 <- setValues(imgBase,l09); r10 <- setValues(imgBase,l10); r11 <- setValues(imgBase,l11); r12 <- setValues(imgBase,l12)
  
  
  # Save
  pheDirYear <- paste0(pheDir,'/',phenYrs[yToDo])
  if (!dir.exists(pheDirYear)) {dir.create(pheDirYear)}
  
  writeRaster(r01,filename=paste0(pheDirYear,'/01_',phenYrs[yToDo],'_',productTable$short_name[ 1],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r02,filename=paste0(pheDirYear,'/02_',phenYrs[yToDo],'_',productTable$short_name[ 2],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r03,filename=paste0(pheDirYear,'/03_',phenYrs[yToDo],'_',productTable$short_name[ 3],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r04,filename=paste0(pheDirYear,'/04_',phenYrs[yToDo],'_',productTable$short_name[ 4],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r05,filename=paste0(pheDirYear,'/05_',phenYrs[yToDo],'_',productTable$short_name[ 5],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r06,filename=paste0(pheDirYear,'/06_',phenYrs[yToDo],'_',productTable$short_name[ 6],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r07,filename=paste0(pheDirYear,'/07_',phenYrs[yToDo],'_',productTable$short_name[ 7],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r08,filename=paste0(pheDirYear,'/08_',phenYrs[yToDo],'_',productTable$short_name[ 8],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r09,filename=paste0(pheDirYear,'/09_',phenYrs[yToDo],'_',productTable$short_name[ 9],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r10,filename=paste0(pheDirYear,'/10_',phenYrs[yToDo],'_',productTable$short_name[10],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r11,filename=paste0(pheDirYear,'/11_',phenYrs[yToDo],'_',productTable$short_name[11],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r12,filename=paste0(pheDirYear,'/12_',phenYrs[yToDo],'_',productTable$short_name[12],'.tif'), format="GTiff", overwrite=TRUE)
  
  
  print(yToDo)
}









