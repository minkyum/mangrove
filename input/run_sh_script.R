###############################
library(rjson)
params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/PLSP_Parameters.json')

geojsonList <- list.files(params$setup$geojsonDir)
sites <- 1:length(geojsonList)


###############################
### 01_Make image mosaic
setwd(paste0(params$setup$logDir,'01'))
for(numSite in 1){
  system(paste('qsub -V -pe omp ',params$setup$numCores,' -l h_rt=12:00:00 ',params$setup$rScripts,'run_script_01.sh ',numSite,sep=''))
}


### 01_1_Make water mask raster
setwd(paste0(params$setup$logDir,'01_1'))
for(numSite in 1){
  nn <- sprintf('%03d',numSite)
  for(ck in 1:24){
    system(paste('qsub -V -pe omp 2 -l h_rt=12:00:00 ',params$setup$rScripts,'run_script_01_1.sh ',nn,ck,sep=''))
  }
}


## 02_Make image chunks
setwd(paste0(params$setup$logDir,'02'))
for(numSite in 1){
  nn <- sprintf('%03d',numSite)
  for(cc in 1:params$setup$numChunks){
    system(paste('qsub -V -l h_rt=12:00:00 ',params$setup$rScripts,'run_script_02.sh ',nn,cc,sep=''))
  }
}
# 2-1
# setwd(paste0(params$setup$logDir,'02'))
# for(numSite in 53:54){
#   # system(paste0('qsub -V -pe omp ',params$setup$numCores,' -l h_rt=12:00:00 ',params$setup$rScripts,'run_script_02.sh ',numSite))
#   system(paste0('qsub -V -pe omp 28 -l h_rt=12:00:00 ',params$setup$rScripts,'run_script_02.sh ',numSite))
# }


### 03_Make pheno-metrics chunks 
setwd(paste0(params$setup$logDir,'03'))
for(numSite in 1){
  nn <- sprintf('%03d',numSite)
  for(cc in 1:params$setup$numChunks){
    system(paste('qsub -V -l h_rt=12:00:00 ',params$setup$rScripts,'run_script_03.sh ',nn,cc,sep=''))
  }
}


### 04_Generate TIFF outputs
setwd(paste0(params$setup$logDir,'04'))
for(numSite in 1){
  system(paste('qsub -V -pe omp 4 -l h_rt=12:00:00 ',params$setup$rScripts,'run_script_04.sh ',numSite,sep=''))  
}
# Remove temporary files
for(numSite in c(1:104)){
  strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
  print(strSite)
  # tempDir <- paste0(params$setup$outDir,strSite,'/chunk')
  # system(paste0('rm -r ',tempDir))
  tempFiles <- list.files(paste0(params$setup$workDir,'Product_010/',strSite),pattern=glob2rx('*xml'),recursive=T,full.names=T)
  for(i in 1:length(tempFiles)){
    system(paste0('rm -r ',tempFiles[i]))
  }
}


### 05_Greate netCDF 
setwd(paste0(params$setup$logDir,'05'))
for(numSite in c(1:104)){
  system(paste('qsub -V -pe omp 2 -l h_rt=12:00:00 ',params$setup$rScripts,'run_script_05.sh ',numSite,sep=''))  
}
