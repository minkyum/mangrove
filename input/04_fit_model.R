library(sp)
library(raster)
library(terra)
library(sf)
library(rjson)

library(randomForest)



###############################
args <- commandArgs()
print(args)

numSite <- as.numeric(substr(args[3],1,3))
yToDo   <- as.numeric(substr(args[3],4,7))
# numSite <- 1; yToDo <- 1


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
# Load training data
t1 <- vect('/projectnb/modislc/users/mkmoon/mangrove/PLCM/shp/training_1/training.shp')   # crop
t2 <- vect('/projectnb/modislc/users/mkmoon/mangrove/PLCM/shp/training_2/training_2.shp') # buildings
t3 <- vect('/projectnb/modislc/users/mkmoon/mangrove/PLCM/shp/training_3/training_3.shp') # terrestrial vegetation
t4 <- vect('/projectnb/modislc/users/mkmoon/mangrove/PLCM/shp/training_4/training_4.shp') # shallow water
t5 <- vect('/projectnb/modislc/users/mkmoon/mangrove/PLCM/shp/training_5/training_5.shp') # deep water
t6 <- vect('/projectnb/modislc/users/mkmoon/mangrove/PLCM/shp/training_6/training_6.shp') # mangrove

t1 <- project(t1,crs(imgBase)); t2 <- project(t2,crs(imgBase)); t3 <- project(t3,crs(imgBase))
t4 <- project(t4,crs(imgBase)); t5 <- project(t5,crs(imgBase)); t6 <- project(t6,crs(imgBase))

pt1 <- extract(imgNum,t1); pt2 <- extract(imgNum,t2); pt3 <- extract(imgNum,t3)
pt4 <- extract(imgNum,t4); pt5 <- extract(imgNum,t5); pt6 <- extract(imgNum,t6)


# Divide the points into training, validation, and testing; 8:1:1
set.seed(456789)
pt11s <- sample(1:nrow(pt1),round(nrow(pt1)*0.8)); pt11  <- pt1[pt11s,]; ptt   <- pt1[-pt11s,] # T1
pt12s <- sample(1:nrow(ptt),round(nrow(ptt)*0.5)); pt12  <- ptt[pt12s,]; pt13  <- ptt[-pt12s,]
pt21s <- sample(1:nrow(pt2),round(nrow(pt2)*0.8)); pt21  <- pt2[pt21s,]; ptt   <- pt1[-pt21s,] # T2
pt22s <- sample(1:nrow(ptt),round(nrow(ptt)*0.5)); pt22  <- ptt[pt22s,]; pt23  <- ptt[-pt22s,]
pt31s <- sample(1:nrow(pt3),round(nrow(pt3)*0.8)); pt31  <- pt3[pt31s,]; ptt   <- pt1[-pt31s,] # T3
pt32s <- sample(1:nrow(ptt),round(nrow(ptt)*0.5)); pt32  <- ptt[pt32s,]; pt33  <- ptt[-pt32s,]
pt41s <- sample(1:nrow(pt4),round(nrow(pt4)*0.8)); pt41  <- pt4[pt41s,]; ptt   <- pt1[-pt41s,] # T4
pt42s <- sample(1:nrow(ptt),round(nrow(ptt)*0.5)); pt42  <- ptt[pt42s,]; pt43  <- ptt[-pt42s,]
pt51s <- sample(1:nrow(pt5),round(nrow(pt5)*0.8)); pt51  <- pt5[pt51s,]; ptt   <- pt1[-pt51s,] # T5
pt52s <- sample(1:nrow(ptt),round(nrow(ptt)*0.5)); pt52  <- ptt[pt52s,]; pt53  <- ptt[-pt52s,]
pt61s <- sample(1:nrow(pt6),round(nrow(pt6)*0.8)); pt61  <- pt6[pt61s,]; ptt   <- pt1[-pt61s,] # T6
pt62s <- sample(1:nrow(ptt),round(nrow(ptt)*0.5)); pt62  <- ptt[pt62s,]; pt63  <- ptt[-pt62s,]



########################################
# Load features
pheDir <- paste0(params$setup$outDir,strSite,'/feat')
if (!dir.exists(pheDir)) {dir.create(pheDir)}

Fmat <- matrix(NA,numPix,49) 

for(i in 1:numChunks){
  cc <- sprintf('%03d',i)
  cfile <- paste0(ckPheDir,'/chunk_feat_',cc,'.rda') 
  log <- try(load(cfile),silent=F)
  if (inherits(log, 'try-error')) next 
    
  if(i==numChunks){chunks <- c((chunk*(i-1)+1):numPix)
  }else{chunks <- c((chunk*(i-1)+1):(chunk*i))}
    
  chunkStart <- chunks[1];  chunkEnd <- chunks[length(chunks)]
    
  Fmat[chunkStart:chunkEnd,] <- f_mat[,(58*(yToDo-1)+c(1:48,58))]
  
 print(i)
}

Fmat11 <- Fmat[pt11[,2],]; Fmat12 <- Fmat[pt12[,2],]; Fmat13 <- Fmat[pt13[,2],]
Fmat21 <- Fmat[pt21[,2],]; Fmat22 <- Fmat[pt22[,2],]; Fmat23 <- Fmat[pt23[,2],]
Fmat31 <- Fmat[pt31[,2],]; Fmat32 <- Fmat[pt32[,2],]; Fmat33 <- Fmat[pt33[,2],]
Fmat41 <- Fmat[pt41[,2],]; Fmat42 <- Fmat[pt42[,2],]; Fmat43 <- Fmat[pt43[,2],]
Fmat51 <- Fmat[pt51[,2],]; Fmat52 <- Fmat[pt52[,2],]; Fmat53 <- Fmat[pt53[,2],]
Fmat61 <- Fmat[pt61[,2],]; Fmat62 <- Fmat[pt62[,2],]; Fmat63 <- Fmat[pt63[,2],]

lct <- c(rep(dim(Fmat11)[1],1),rep(dim(Fmat21)[1],2),rep(dim(Fmat31)[1],3),
         rep(dim(Fmat41)[1],4),rep(dim(Fmat51)[1],5),rep(dim(Fmat61)[1],6))

mat_pca <- as.data.frame(rbind(Fmat11,Fmat21,Fmat31,Fmat41,Fmat51,Fmat61))
pca <- princomp(~.,data=mat_pca,cor=T)

pca$scores[1:5, ]
plot(pca$scores[,1],pca$scores[,2],col=lct)

rand

