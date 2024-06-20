library(gdalUtils)
library(raster)
library(rgdal)

library(RColorBrewer)




###
mycolRamp = colorRampPalette(c('#640065ff','#8200bcff','#1f00ffff','#00ffebff','#52ff00ff',
                               '#ffff00ff','#ff4200ff','#ff0000ff','#c90000ff','#8d0000ff','#610000ff'))
setwd('/projectnb/modislc/users/mkmoon/Planet/data_paper/figure/')

png(filename='site_eg_legend_2.png',width=5,height=2,unit='in',res=600,bg=NA)

plot(raster(matrix(seq(0,1,length.out=100),10,10)),
     legend.only=T,
     col=mycolRamp(300),
     zlim=c(60,240),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=T,
     smallplot=c(0.1,0.9,0.5,0.55),
     axis.args=list(at=c(60,120,180,240),
                    labels=c(60,120,180,240),
                    cex.axis=1.5,font=2,col.axis='black',col='black'),
     # legend.args=list(text='Thermal forcing',side=2,font=1,line=0.5,cex=1.7)
     )

dev.off()



###
spec <- brewer.pal(9,'Greens')
mycolRamp = colorRampPalette(spec)

setwd('/projectnb/modislc/users/mkmoon/Planet/data_paper/figure/')
png(filename='site_eg_legend_1_evi.png',width=5,height=2,unit='in',res=600,bg=NA)

plot(raster(matrix(seq(0,1,length.out=100),10,10)),
     legend.only=T,
     col=mycolRamp(300),
     zlim=c(1,6),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=T,
     smallplot=c(0.1,0.9,0.5,0.55),
     axis.args=list(at=c(1,2,3,4,5,6),
                    labels=c(0.1,0.2,0.3,0.4,0.5,0.6),
                    cex.axis=1.5,font=2),
     # legend.args=list(text='Thermal forcing',side=2,font=1,line=0.5,cex=1.7)
)

dev.off()











