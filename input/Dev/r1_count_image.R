metTable <- read.csv('/projectnb/modislc/users/mkmoon/Planet/data_paper/data/site_full_w_pc_hls_PCrevised_ext.csv')

numImage <- matrix(NA,104,1)
for(i in 1:104){
  a1 <- list.files(paste0('/projectnb/planet/neon/',metTable$Site[i]),pattern='*_SR*',recursive=T)
  a2 <- list.files(paste0('/projectnb/planet/amflx/',metTable$Site[i]),pattern='*_SR*',recursive=T)
  b1 <- list.files(paste0('/projectnb/modislc/users/mkmoon/Planet/rawImage/rawImage/',metTable$Site[i]),pattern='*_SR*',recursive=T)
  
  numImage[i] <- sum(length(a1),length(a2),length(b1))
  print(paste(i,';',numImage[i]))
}
summary(numImage)
