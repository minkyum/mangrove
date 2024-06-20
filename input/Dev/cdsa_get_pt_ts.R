library(sp)
library(raster)
library(terra)
library(sf)

library(rjson)
library(geojsonR)

library(doMC)
library(doParallel)

library(RColorBrewer)

###############################
args <- commandArgs()
print(args)

numSite <- as.numeric(args[3])
# numSite <- 3

NumSite <- c(14,19,28,87)
numSite <- NumSite[numSite]

########################################
params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/PLSP_Parameters.json')
source(params$setup$rFunctions)


########################################
## Get site name and image directory
geojsonDir <- params$setup$geojsonDir

strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

# ckDir <- paste0(params$setup$outDir,strSite,'/chunk')
ckDir <- paste0('/projectnb/modislc/users/mkmoon/Planet/rawImage/chunks/',strSite)
print(ckDir)



########################################
## 
imgBase <- raster(paste0(params$setup$outDir,strSite,'/base_image.tif'))
numPix <- length(imgBase)
imgNum <- setValues(imgBase,1:numPix)
numChunks <- params$setup$numChunks
chunk <- numPix%/%numChunks

nn <- sprintf('%03d',numSite)
ptShp <- shapefile(paste0('/projectnb/modislc/users/mkmoon/cdsa/plsp/shp/pt_',numSite-1,'.shp'))
ptShp <- spTransform(ptShp,crs(imgBase))

pixNums <- extract(imgNum,ptShp)


##
setwd('/projectnb/modislc/users/mkmoon/cdsa/plsp/figures/')
png(filename=paste0('ts_',strSite,'_eg_1.png'),width=12,height=6.5,unit='in',res=300)

par(mfrow=c(2,2),oma=c(0,1,0,0),mar=c(2,4,1,1),mgp=c(2.0,1,0))
for(ppp in 1:4){
  pixNum  <- pixNums[ppp]

  ckNum <- sprintf('%03d',(pixNum%/%chunk+1))
  file <- list.files(path=ckDir,pattern=glob2rx(paste0('*',ckNum,'.rda')),full.names=T)

  load(file)

  # repRo <- which(dates >= '2020-1-1' & dates <= '2020-7-1')
  # repDa <- dates[which(dates >= '2020-1-1' & dates <= '2020-7-1')]+366
  # 
  # dates <- c(dates,repDa)
  # 
  # band1 <- cbind(band1,band1[,repRo])
  # band2 <- cbind(band2,band2[,repRo])
  # band3 <- cbind(band3,band3[,repRo])
  # band4 <- cbind(band4,band4[,repRo])

  blue  <- band1[pixNum%%chunk,]
  green <- band2[pixNum%%chunk,]
  red   <- band3[pixNum%%chunk,]
  nir   <- band4[pixNum%%chunk,]
  phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr

  blue <- blue/10000; green <- green/10000; red <- red/10000; nir <- nir/10000
  
  pheno_pars <- params$phenology_parameters
  qa_pars    <- params$qa_parameters
  
  
  vi     <- 2.5*(nir - red) / (nir + 2.4*red + 1)
  
  # # Potential water
  # if( sum(vi<0,na.rm=T)/sum(!is.na(vi)) > pheno_pars$sumNegVIthresh & waterMask > pheno_pars$waterOccuThresh ){
  #   return(rep(c(NA,rep(NA,10),4,rep(NA,10),4,NA),length(phenYrs)))} 
  
  # plot(dates,vi,ylim=c(-0.1,0.6),axes=T,ylab='EVI2',
  #      cex.lab=1.2,cex.axis=1.2)
  
  
  # # Negative VIs
  # vi[vi<0]   <- NA

  # Spikes
  spikes <- CheckSpike_MultiBand(blue, red, vi, dates, pheno_pars)
  vi[spikes] <- NA

  #
  dormIms <- dates >= pheno_pars$dormStart & dates <= pheno_pars$dormEnd
  vi_dorm <- quantile(vi[dormIms & vi>0],probs=pheno_pars$dormantQuantile,na.rm=T)   #Calc vi dormant value

  vi[vi < vi_dorm] <- vi_dorm

  ######################################################
  # Gap fill
  splineStart <- as.Date(as.Date(paste0(phenYrs,'-01-01')) - pheno_pars$splineBuffer)
  numDaysFit  <-  365 + (pheno_pars$splineBuffer * 2)
  splineEnd <- splineStart+(numDaysFit-1)

  daysVec <- 1:numDaysFit
  inYear <- daysVec > pheno_pars$splineBuffer & daysVec <= (pheno_pars$splineBuffer+365)


  #Determine gaps that require filling
  gDates <- dates[!is.na(vi)]
  dDiff <- diff(gDates) > pheno_pars$gapLengthToFill     #Gaps greater than 20 days will be filled
  dStart <- gDates[c(dDiff,FALSE)]
  dEnd <- gDates[c(FALSE,dDiff)]


  #Locate gaps in date vector
  all_dates <- seq(min(splineStart), max(splineEnd), by="day")

  fill_locations <- matrix(FALSE,length(all_dates))
  for (d in 1:length(dStart)) {
    fill_locations[all_dates >= dStart[d] & all_dates < dEnd[d]] <- TRUE}

  fill_dates <- all_dates[fill_locations]

  yToDo <- 1:length(phenYrs)
  yrsWithGaps <- c()
  for (y in yToDo) {
    pred_dates <- seq(splineStart[y], splineEnd[y], by="day")
    if (sum(pred_dates %in% fill_dates) > 0) {yrsWithGaps <- c(yrsWithGaps,phenYrs[y])}
  }


  #If there are gaps to be filled, then we will spline all years.
  #If not, just spline product years
  yrs <- phenYrs

  numYrs <- length(yrs)
  daysVec <- 1:numDaysFit
  vecLength <- numDaysFit*numYrs


  #First, we will fit splines to each year invidually
  #To line up observations from each year, we will create a matrix for vi and each band (numDaysFit x numYears)

  smoothMat <- matrix(NA, numDaysFit, numYrs)
  maskMat <- matrix(0, numDaysFit, numYrs)
  fillMat <- smoothMat
  baseWeights <- maskMat

  for (y in 1:numYrs) {
    #Use try statement, because we don't want to stop processing if only an error in one year
    try({

      dateRange <- dates >= splineStart[y] & dates <= splineEnd[y] & !is.na(vi)

      dateSub <- dates[dateRange]
      viSub <- vi[dateRange]

      #Get weights
      weights <- matrix(1,length(dateSub))

      pred_dates <- seq(splineStart[y], splineEnd[y], by="day")

      #Assign weights and run cubic spline
      smoothed <- Smooth_VI(viSub, dateSub, pred_dates, weights, pheno_pars, vi_dorm)


      #Mask spline in gaps, and before/after first/last image
      maskMat[fill_locations[all_dates %in% pred_dates],y] <- 1    #Mask spline in gaps
      maskMat[pred_dates < dateSub[1],y] <- 1                      #Mask spline before first image and after last image
      maskMat[pred_dates > dateSub[length(dateSub)],y] <- 1

      #Mask spline in the buffer years (only interested in comparing splines in target year)
      maskMat[format(pred_dates,'%Y') != yrs[y],y]  <- 1

      fillDs <- pred_dates %in% dateSub

      smoothMat[,y] <- smoothed
      baseWeights[fillDs,y] <- weights
      fillMat[fillDs,y] <- viSub

    },silent=TRUE)
  }


  xs <- rep(daysVec,numYrs)
  ys <- matrix(fillMat,vecLength)
  ysGood <- !is.na(ys)
  baseW <- matrix(baseWeights,vecLength)   #Base Weights are 1=clear observation

  smoothMat_Masked <- smoothMat
  maskMat <- as.logical(maskMat)
  smoothMat_Masked[maskMat] <- NA


  #Loop through years, compare spline to other years, weight each year based on similarity, fit spline, calculate phenology

  weightArray <- calculateWeights(smoothMat_Masked, numDaysFit, numYrs, pheno_pars)

  prevYear <- daysVec <= pheno_pars$splineBuffer
  inYear <- daysVec > pheno_pars$splineBuffer & daysVec <= (pheno_pars$splineBuffer+365)
  nextYear <- daysVec > (pheno_pars$splineBuffer+365)




  y <- 4
  # for(y in 1:length(splineStart)){
    pred_dates <- seq(splineStart[y], splineEnd[y], by="day")


    if (yrs[y] %in% yrsWithGaps) {

      indPrev <- y-1; indPrev[indPrev<1] <- 1
      indNext <- y+1; indNext[indNext>numYrs] <- numYrs

      weights <- rbind(weightArray[prevYear,,indPrev],
                       weightArray[inYear,,y],
                       weightArray[nextYear,,indNext])

      #Where are the gaps?
      toFill <- fill_locations[all_dates %in% pred_dates]

      weights[!toFill,] <- 0     #Set weight to zero for observations that aren't in a gap
      weights[,y] <- 1           #Set weights in target year to 1


      #Now that we have weights, calculate phenology
      #######################
      weights <- matrix(weights,vecLength) * baseW   #Multiple weights by base weight
      theInds <- ysGood & weights > 0
      xs_sub <- xs[theInds]; w_sub <- weights[theInds]
      smoothed_vi <- Smooth_VI(ys[theInds], xs_sub, daysVec, w_sub, pheno_pars, vi_dorm)  #Fit spline

    } else {

      #Variables needed for next steps if the above gap filling was not done
      theInds <- matrix(FALSE,length(ysGood))
      theInds[((y-1)*numDaysFit+1):(y*numDaysFit)] <- TRUE
      xs_sub <- xs[theInds]; w_sub <- baseW[theInds]

      smoothed_vi <- smoothMat[,y]   #if no gaps to fill, just use existing spline
    }


    # Number of clear observation
    filled_vi <- fillMat[,y]
    filled_vi[baseWeights[,y] < 1] <- NA    #If weight is less than 1, implies it is a snow-fill, and we don't want to count snow-filled as a valid observation. So set to NA.
    numObs <- sum(!is.na(filled_vi) & inYear)   #Number of observations in year

    #
    viSub   <- filled_vi[!is.na(filled_vi)]
    dateSub <- pred_dates[!is.na(filled_vi)]

    # print(sum(smoothed_vi<0.03)/length(smoothed_vi))

    plot(18262:18627,seq(0,0.6,length.out=366),col='white',
         ylim=c(0,0.6),xlab='',ylab='EVI2',cex.lab=1.2,axe=F)
    axis(1,at=c(18262,18262+90,18262+180,18262+270,18262+360),
         c('Jan','Apr','Jul','Oct','Dec'),cex.axis=1.2)
    axis(2,c(0,0.2,0.4,0.6),cex.axis=1.2)
    points(dateSub,viSub)
    points(pred_dates,smoothed_vi,type='l',lwd=2)

    # points(dateSub,viSub,col='red')
    # points(pred_dates,smoothed_vi,col='red')
    # # text(18090,0.75,round(vi_dorm1,3),pos=4,col='red',cex=1.5)
    # text(18090,0.65,round(vi_dorm,3),pos=4,cex=1.5)
    # text(18090,0.55,round(max(smoothed_vi)-min(smoothed_vi),3),pos=4,cex=1.5)
  # }

    
}

dev.off()





###
mycolRamp = colorRampPalette(c('#640065ff','#8200bcff','#1f00ffff','#00ffebff','#52ff00ff',
                               '#ffff00ff','#ff4200ff','#ff0000ff','#c90000ff','#8d0000ff','#610000ff'))
mycolRamp = colorRampPalette(brewer.pal(11,'Spectral'))

setwd('/projectnb/modislc/users/mkmoon/cdsa/plsp/figures/')
png(filename='site_eg_legend_1.png',width=6.5,height=2,unit='in',res=600,bg=NA)

plot(raster(matrix(seq(0,1,length.out=100),10,10)),
     legend.only=T,
     col=mycolRamp(300),
     zlim=c(135,160),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=T,
     smallplot=c(0.1,0.9,0.5,0.55),
     axis.args=list(at=c(135,140,145,150,155,160),
                    labels=c(135,140,145,150,155,160),
                    cex.axis=1.5,font=1),
     # legend.args=list(text='Thermal forcing',side=2,font=1,line=0.5,cex=1.7)
)
dev.off()


mycolRamp = colorRampPalette(brewer.pal(9,'Greens'))

setwd('/projectnb/modislc/users/mkmoon/cdsa/plsp/figures/')
png(filename='site_eg_legend_2.png',width=6.5,height=2,unit='in',res=600,bg=NA)

plot(raster(matrix(seq(0,1,length.out=100),10,10)),
     legend.only=T,
     col=mycolRamp(300),
     zlim=c(2000,5000),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=T,
     smallplot=c(0.1,0.9,0.5,0.55),
     axis.args=list(at=c(2000,3000,4000,5000),
                    labels=c(0.2,0.3,0.4,0.5),
                    cex.axis=1.5,font=1),
     # legend.args=list(text='Thermal forcing',side=2,font=1,line=0.5,cex=1.7)
)
dev.off()
