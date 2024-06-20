library(sp)
library(raster)
library(terra)
library(sf) 

library(RColorBrewer)
library(viridis)


########################################
args <- commandArgs()
print(args)

ss <- as.numeric(substr(args[3],1,3))
yy <- as.numeric(substr(args[3],4,7))
vv <- as.numeric(substr(args[3],8,9))

# ss <- 4; yy <- 2019; vv <- 6


vari <- c('OGI','50PCGI','OGMx','Peak','OGD','50PCGD','OGMn','EVImax','EVIamp','EVIarea')

options(warn=-1)


########################################
###
tileHLS <- c('06WVT','06VWR','06VUR','06WVB')
path <- paste0('/projectnb/modislc/users/mkmoon/MuSLI/V1_0/From_AWS/product_v011/',tileHLS[ss])
file <- list.files(path,pattern=glob2rx(paste0('*',yy,'*')),full.names=T)
rastH <- raster(file[1],varname=vari[vv])

###
SS <- c(14,19,27,87)
ss <- SS[ss]
imgs <- list.files('/projectnb/planet/PLSP/Product_010',pattern=glob2rx(paste0('*',yy,'_',vari[vv],'.tif')),full.names=T,recursive=T)
site <- unlist(strsplit(imgs[ss],'/'))[6]
print(site)

rastP <- raster(imgs[ss])

###
rastH1 <- crop(rastH,rastP)

rastP1 <- aggregate(rastP,3,fun=median)

r1 <- rastP1
r1[r1 < quantile(rastP1,0.05)] <- quantile(rastP1,0.05)
r1[r1 > quantile(rastP1,0.95)] <- quantile(rastP1,0.95)

r2 <- rastH1
r2[r2 < quantile(rastP,0.05)] <- quantile(rastP,0.05)
r2[r2 > quantile(rastP,0.95)] <- quantile(rastP,0.95)


########################################
setwd('/projectnb/modislc/users/mkmoon/cdsa/plsp/figures/')
png(filename=paste0('phe_',site,'_',vv,'_',vari[vv],'.png'),width=11,height=5,unit='in',res=300)

par(mfrow=c(1,2),oma=c(0,0,0,1),mar=c(1,1,1.5,1),mgp=c(2.5,1,0))

mycolRamp = colorRampPalette(c('#640065ff','#8200bcff','#1f00ffff','#00ffebff','#52ff00ff',
                               '#ffff00ff','#ff4200ff','#ff0000ff','#c90000ff','#8d0000ff','#610000ff'))
mycolRamp = colorRampPalette(rev(brewer.pal(11,'Spectral')))

plot(r1,col=mycolRamp(100),interpolate=T,legend=F,axes=F,box=F,main='Planet')
plot(r2,col=mycolRamp(100),interpolate=T,legend=T,axes=F,box=F,main='HLS')
  
dev.off()


