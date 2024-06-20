
par(mfrow=c(3,4),oma=c(2,2,0,0),mar=c(1,1,1,1))
for(ss in c(1:104)){
  dirs <- list.files('/projectnb/planet/PLSP/Img_cliped',full.names=T)
  files <- list.files(paste0(dirs[ss],'/mosaic'),pattern=glob2rx('*mosaic.tif'))
  yy <- as.Date(paste0(substr(files,1,4),'-',substr(files,5,6),'-',substr(files,7,8))) 
  plot(1:length(yy),as.Date(sort(yy),origin='1970-1-1'),ylim=c(16800,19030))
  abline(h=19015)
  legend('topleft',c(strsplit(dirs[ss],'/')[[1]][6],length(yy)),cex=1.2,bty='n')
  
  print(ss)  
}


