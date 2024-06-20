library(sp)
library(raster)
library(terra)
library(sf)

library(rjson)
library(geojsonR)

library(doMC)
library(doParallel)


###############################
args <- commandArgs()
print(args)

numSite <- as.numeric(substr(args[3],1,3))
cc      <- as.numeric(substr(args[3],4,6))
# numSite <- 1; cc <- 50



###############################
params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/PLSP_Parameters.json')
source(params$setup$rFunctions)



########################################
## Get site name, image directory and coordinate
strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

ckDir <- paste0(params$setup$outDir,strSite,'/chunk')
# ckDir <- paste0('/projectnb/modislc/users/mkmoon/Planet/rawImage/chunks/',strSite)
print(ckDir)

ckNum <- sprintf('%03d',cc)
file <- list.files(path=ckDir,pattern=glob2rx(paste0('*',ckNum,'.rda')),full.names=T)

load(file)


# ## Fill the first half of 2021 using the first half of 2020
# repRo <- which(dates >= '2021-1-29' & dates <= '2021-7-1')
# repDa <- dates[which(dates >= '2021-1-29' & dates <= '2021-7-1')]+365
# 
# dates <- c(dates,repDa)
# 
# band1 <- cbind(band1,band1[,repRo])
# band2 <- cbind(band2,band2[,repRo])
# band3 <- cbind(band3,band3[,repRo])
# band4 <- cbind(band4,band4[,repRo])


## Load water mask
waterRaster <- raster(paste0(params$setup$outDir,strSite,'/water_mask_30_1.tif'))

numCk <- params$setup$numChunks
chunk <- length(waterRaster)%/%numCk
if(cc==numCk){
  chunks <- c((chunk*(cc-1)+1):length(waterRaster))
}else{
  chunks <- c((chunk*(cc-1)+1):(chunk*cc))
}
waterMask <- values(waterRaster)[chunks]


##########################################
numPix <- dim(band1)[1]
phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr

pheno_mat <- matrix(NA,numPix,24*length(phenYrs))

for (i in 1:numPix){

  pheno_mat[i,] <- DoPhenologyPlanet(band1[i,],band2[i,],band3[i,],band4[i,],dates,phenYrs,params,waterMask[i])
  
  if(i%%10000==0) print(i)
}


# Save
ckPheDir <- paste0(params$setup$outDir,strSite,'/chunk_phe')
# ckPheDir <- paste0('/projectnb/modislc/users/mkmoon/cdsa/plsp/ImgCliped/',strSite,'/chunk_phe')
if (!dir.exists(ckPheDir)) {dir.create(ckPheDir)}

save(pheno_mat,file=paste0(ckPheDir,'/chunk_phe_',ckNum,'.rda'))




