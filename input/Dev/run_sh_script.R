###############################
### Download images for 2021
setwd('/projectnb/modislc/users/mkmoon/Planet/rawImage/runLogs/')
for(ss in c(4)){
  system(paste0('qsub -V -N Rsite',sprintf('%02d',ss),'_1 -l h_rt=36:00:00 /usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/Dev/download.sh ',(ss-1)))  
}
####
#
# Suspicious: 86
#
# Incomplete: 13 (a chunk copied from 12 due to incomplete)
#
#
library(rjson)
json <- fromJSON(file='/projectnb/planet/PLSP/Img_raw/rawImage/Bartlett_Experimental_Forest_NEON/order_Bartlett_Experimental_Forest_NEON_chunk_0_2021_2023.json')

### Reorganize data
setwd('/projectnb/modislc/users/mkmoon/Planet/rawImage/runLogs/')
dirs <- list.files('/projectnb/modislc/users/mkmoon/Planet/rawImage/ord_test/')
for(ss in c(1)){
  system(paste0('qsub -V -N OG',sprintf('%02d',ss),' -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/Dev/org_script.sh ',dirs[ss]))  
}


### Move data
setwd('/projectnb/modislc/users/mkmoon/Planet/rawImage/runLogs/')
for(ss in c(24)){
  system(paste0('qsub -V -N MV',sprintf('%02d',ss),' -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/Dev/move.sh ',ss))  
}



###############################
### QC images
setwd('/projectnb/planet/PLSP/Output_images/runLogs/')
for(ss in c(1:104)){
  for(yy in 2017:2021){
    numSite <- sprintf('%03d',ss)
    system(paste('qsub -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/Dev/run_script_01.sh ',numSite,yy,sep=''))  
  }
}



###############################
### Comparison with HLS
setwd('/projectnb/modislc/users/mkmoon/Planet/data_paper/data/comp_hls/runLogs/')
for(ss in c(100:104)){
  for(yy in 2017:2019){
    for(vv in 1:6){
      numSite <- sprintf('%03d',ss)
      system(paste('qsub -V -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/Dev/run_script_02.sh ',numSite,yy,vv,sep=''))    
    }
  }
}


### Comparison with PhenoCam
setwd('/projectnb/modislc/users/mkmoon/Planet/data_paper/data/comp_phenocam/runLogs/')
for(ss in 1:104){
  system(paste('qsub -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/Dev/run_script_03.sh ',ss,sep=''))    
}



###############################
### time series for ncar
setwd('/projectnb/modislc/users/mkmoon/Planet/data/ncar/runLogs/')
for(ss in 1:10){
  system(paste('qsub -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/Dev/run_script_03.sh ',ss,sep=''))    
}


###############################
### file size
setwd('/projectnb/planet/PLSP/runLogs/')
for(ss in c('amflx','neon')){
  system(paste('qsub -V -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/Dev/get_file_size.sh ',ss,sep=''))    
}

###############################
### zip files
setwd('/projectnb/planet/PLSP/runLogs/')
dirs <- list.files('/projectnb/planet/PLSP/Product_010')
for(ss in c(2:99)){
  system(paste('qsub -V -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/Dev/zip_dirs.sh ',dirs[ss],sep=''))    
}

###############################
### PhenoCam sites - Supporting table
setwd('/projectnb/modislc/users/mkmoon/Planet/data_paper/data/comp_phenocam/runLogs/')
for(ss in 1:104){
  system(paste('qsub -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/Dev/run_script_04.sh ',ss,sep=''))    
}



###############################
### CDSA PLSP and MSLSP comparison
setwd('/projectnb/modislc/users/mkmoon/cdsa/plsp/runLogs/')
for(ss in 1:4){
  for(yy in 2019){
    for(vv in 1:10){
      numSite <- sprintf('%03d',ss)
      system(paste('qsub -V -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/Dev/run_script_11.sh ',numSite,yy,vv,sep=''))    
    }
  }
}


###############################
### Korea points
setwd('/projectnb/modislc/users/mkmoon/KoreaProject/runLogs/')
ctcls <- c('dbf_class1','dbf_class2','dbf_class3','dbf_class4','dbf_class5',
           'enf_class1','enf_class2','enf_class3','enf_class4','enf_class5')
for(i in 1:10){
  system(paste('qsub -V -pe omp 28 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/Dev/run_script_20.sh ',ctcls[i],sep=''))    
}






