library(sp)
library(raster)
library(terra)

library(rjson)
library(geojsonR)

library(doMC)
library(doParallel)

###############################
args <- commandArgs()
print(args)

numSite <- as.numeric(args[3])
# numSite <- 1


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
ptShp <- shapefile(paste0(params$setup$workDir,'Shape_file/pt_',nn,'.shp'))
ptShp <- spTransform(ptShp,crs(imgBase))

pixNums <- extract(imgNum,ptShp)


##
par(mfrow=c(2,2))
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
  
  plot(dates,vi)
  
  # # # Negative VIs
  # # vi[vi<0]   <- NA
  # 
  # # Spikes
  # spikes <- CheckSpike_MultiBand(blue, red, vi, dates, pheno_pars)
  # vi[spikes] <- NA
  # 
  # #
  # dormIms <- dates >= pheno_pars$dormStart & dates <= pheno_pars$dormEnd
  # vi_dorm <- quantile(vi[dormIms & vi>0],probs=pheno_pars$dormantQuantile,na.rm=T)   #Calc vi dormant value
  # 
  # vi[vi < vi_dorm] <- vi_dorm
  # 
  # ######################################################
  # # Gap fill
  # splineStart <- as.Date(as.Date(paste0(phenYrs,'-01-01')) - pheno_pars$splineBuffer) 
  # numDaysFit  <-  365 + (pheno_pars$splineBuffer * 2)    
  # splineEnd <- splineStart+(numDaysFit-1)
  # 
  # daysVec <- 1:numDaysFit
  # inYear <- daysVec > pheno_pars$splineBuffer & daysVec <= (pheno_pars$splineBuffer+365)
  # 
  # 
  # #Determine gaps that require filling
  # gDates <- dates[!is.na(vi)]  
  # dDiff <- diff(gDates) > pheno_pars$gapLengthToFill     #Gaps greater than 20 days will be filled
  # dStart <- gDates[c(dDiff,FALSE)]
  # dEnd <- gDates[c(FALSE,dDiff)]
  # 
  # 
  # #Locate gaps in date vector
  # all_dates <- seq(min(splineStart), max(splineEnd), by="day")
  # 
  # fill_locations <- matrix(FALSE,length(all_dates))
  # for (d in 1:length(dStart)) {
  #   fill_locations[all_dates >= dStart[d] & all_dates < dEnd[d]] <- TRUE}
  # 
  # fill_dates <- all_dates[fill_locations]
  # 
  # yToDo <- 1:length(phenYrs)
  # yrsWithGaps <- c()
  # for (y in yToDo) {
  #   pred_dates <- seq(splineStart[y], splineEnd[y], by="day")
  #   if (sum(pred_dates %in% fill_dates) > 0) {yrsWithGaps <- c(yrsWithGaps,phenYrs[y])}
  # }
  # 
  # 
  # #If there are gaps to be filled, then we will spline all years. 
  # #If not, just spline product years
  # yrs <- phenYrs
  # 
  # numYrs <- length(yrs)
  # daysVec <- 1:numDaysFit
  # vecLength <- numDaysFit*numYrs
  # 
  # 
  # #First, we will fit splines to each year invidually
  # #To line up observations from each year, we will create a matrix for vi and each band (numDaysFit x numYears)
  # 
  # smoothMat <- matrix(NA, numDaysFit, numYrs)
  # maskMat <- matrix(0, numDaysFit, numYrs)
  # fillMat <- smoothMat
  # baseWeights <- maskMat
  # 
  # for (y in 1:numYrs) {
  #   #Use try statement, because we don't want to stop processing if only an error in one year
  #   try({
  #     
  #     dateRange <- dates >= splineStart[y] & dates <= splineEnd[y] & !is.na(vi)   
  #     
  #     dateSub <- dates[dateRange]
  #     viSub <- vi[dateRange]
  #     
  #     #Get weights
  #     weights <- matrix(1,length(dateSub))
  #     
  #     pred_dates <- seq(splineStart[y], splineEnd[y], by="day")
  #     
  #     #Assign weights and run cubic spline
  #     smoothed <- Smooth_VI(viSub, dateSub, pred_dates, weights, pheno_pars, vi_dorm)
  #     
  #     
  #     #Mask spline in gaps, and before/after first/last image
  #     maskMat[fill_locations[all_dates %in% pred_dates],y] <- 1    #Mask spline in gaps
  #     maskMat[pred_dates < dateSub[1],y] <- 1                      #Mask spline before first image and after last image
  #     maskMat[pred_dates > dateSub[length(dateSub)],y] <- 1
  #     
  #     #Mask spline in the buffer years (only interested in comparing splines in target year)
  #     maskMat[format(pred_dates,'%Y') != yrs[y],y]  <- 1
  #     
  #     fillDs <- pred_dates %in% dateSub
  #     
  #     smoothMat[,y] <- smoothed
  #     baseWeights[fillDs,y] <- weights
  #     fillMat[fillDs,y] <- viSub
  #     
  #   },silent=TRUE)
  # }
  # 
  # 
  # xs <- rep(daysVec,numYrs)
  # ys <- matrix(fillMat,vecLength)
  # ysGood <- !is.na(ys)
  # baseW <- matrix(baseWeights,vecLength)   #Base Weights are 1=clear observation
  # 
  # smoothMat_Masked <- smoothMat
  # maskMat <- as.logical(maskMat)
  # smoothMat_Masked[maskMat] <- NA
  # 
  # 
  # #Loop through years, compare spline to other years, weight each year based on similarity, fit spline, calculate phenology
  # 
  # weightArray <- calculateWeights(smoothMat_Masked, numDaysFit, numYrs, pheno_pars) 
  # 
  # prevYear <- daysVec <= pheno_pars$splineBuffer
  # inYear <- daysVec > pheno_pars$splineBuffer & daysVec <= (pheno_pars$splineBuffer+365)
  # nextYear <- daysVec > (pheno_pars$splineBuffer+365)
  # 
  # 
  # 
  # 
  # y <- 1
  # 
  # 
  #   pred_dates <- seq(splineStart[y], splineEnd[y], by="day")
  #   
  #   
  #   if (yrs[y] %in% yrsWithGaps) {
  #     
  #     indPrev <- y-1; indPrev[indPrev<1] <- 1
  #     indNext <- y+1; indNext[indNext>numYrs] <- numYrs
  #     
  #     weights <- rbind(weightArray[prevYear,,indPrev],
  #                      weightArray[inYear,,y],
  #                      weightArray[nextYear,,indNext])
  #     
  #     #Where are the gaps?
  #     toFill <- fill_locations[all_dates %in% pred_dates]
  #     
  #     weights[!toFill,] <- 0     #Set weight to zero for observations that aren't in a gap
  #     weights[,y] <- 1           #Set weights in target year to 1
  #     
  #     
  #     #Now that we have weights, calculate phenology
  #     #######################
  #     weights <- matrix(weights,vecLength) * baseW   #Multiple weights by base weight
  #     theInds <- ysGood & weights > 0
  #     xs_sub <- xs[theInds]; w_sub <- weights[theInds]
  #     smoothed_vi <- Smooth_VI(ys[theInds], xs_sub, daysVec, w_sub, pheno_pars, vi_dorm)  #Fit spline
  #     
  #   } else {
  #     
  #     #Variables needed for next steps if the above gap filling was not done
  #     theInds <- matrix(FALSE,length(ysGood))
  #     theInds[((y-1)*numDaysFit+1):(y*numDaysFit)] <- TRUE
  #     xs_sub <- xs[theInds]; w_sub <- baseW[theInds]
  #     
  #     smoothed_vi <- smoothMat[,y]   #if no gaps to fill, just use existing spline
  #   }
  #   
  #   
  #   # Number of clear observation
  #   filled_vi <- fillMat[,y]
  #   filled_vi[baseWeights[,y] < 1] <- NA    #If weight is less than 1, implies it is a snow-fill, and we don't want to count snow-filled as a valid observation. So set to NA.
  #   numObs <- sum(!is.na(filled_vi) & inYear)   #Number of observations in year
  #   
  #   #
  #   viSub   <- filled_vi[!is.na(filled_vi)]
  #   dateSub <- pred_dates[!is.na(filled_vi)]
  # 
  # # print(sum(smoothed_vi<0.03)/length(smoothed_vi))
  # 
  # # plot(dateSub,viSub)
  # # points(pred_dates,smoothed_vi)
  # 
  # # points(dateSub,viSub,col='red')
  # # points(pred_dates,smoothed_vi,col='red')
  # # # text(18090,0.75,round(vi_dorm1,3),pos=4,col='red',cex=1.5)
  # # text(18090,0.65,round(vi_dorm,3),pos=4,cex=1.5)
  # # text(18090,0.55,round(max(smoothed_vi)-min(smoothed_vi),3),pos=4,cex=1.5)
}






##
ckNum <- 100
file <- list.files(path=ckDir,pattern=glob2rx(paste0('*',ckNum,'.rda')),full.names=T)

load(file)

Dates <- dates

qcOutDir <- paste0(params$setup$workDir,'Output_images/ts')
if (!dir.exists(qcOutDir)) {dir.create(qcOutDir)}

sam <- sample(1:dim(band1)[1],300)

for(ppp in 1:300){

  png(filename=paste0(qcOutDir,'/ts_',strSite,'_',ppp,'.png'),width=14,height=8,units='in',res=100)
  par(oma=c(2,2,1,1),mar=c(4,4,2,2))

  blue  <- band1[sam[ppp],]
  green <- band2[sam[ppp],]
  red   <- band3[sam[ppp],]
  nir   <- band4[sam[ppp],]
  phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr


  blue <- blue/10000; green <- green/10000; red <- red/10000; nir <- nir/10000

  vi    <- 2.5*(nir - red) / (nir + 2.4*red + 1)
  ndvi  <- (nir - red) / (nir + red)
  ndsi  <- (green-nir)/(green+nir)
  
  plot(Dates,ndvi,ylim=c(-0.2,0.9),pch=19,cex=0.5,col='blue')
  dormIms   <- Dates >= pheno_pars$dormStart & Dates <= pheno_pars$dormEnd
  ndvi_dorm <- quantile(ndvi[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T)   #Calc vi dormant value
  

    pheno_pars <- params$phenology_parameters
    qa_pars    <- params$qa_parameters

    #filters
    spikes <- CheckSpike_MultiBand(blue, red, vi, Dates, pheno_pars)
    vi[spikes] <- NA
    weights <- rep(1,length(Dates))


    # 3-day composite
    if(pheno_pars$do3dayComposites){
      dates3c <- c()
      vi3c <- c()
      weights3c <- c()

      allStart <- as.Date(as.Date(pheno_pars$dormStart) - pheno_pars$splineBuffer)
      allEnd <- as.Date(as.Date(pheno_pars$dormEnd) + pheno_pars$splineBuffer)
      all_dates <- seq(allStart, allEnd, by="day")

      dateNew <- 1
      for(dateSeq in seq(1,length(all_dates),3)){
        ind <- which(Dates==all_dates[dateSeq]|Dates==all_dates[dateSeq+1]|Dates==all_dates[dateSeq+2])
        if(length(ind)>0 & sum(!is.na(vi[ind]))>0){
          dates3c[dateNew] <- all_dates[dateSeq+1]
          vi3c[dateNew]   <- max(vi[ind],na.rm=T)

          weights[ind][ndsi[ind]>  0] <- 0.5
          weights3c[dateNew] <-  mean(weights[ind],na.rm=T)

          dateNew <- dateNew + 1
        }
      }
      dates    <- as.Date(dates3c,origin='1970-1-1')
      vi      <- vi3c
      weights <- weights3c
    }

    points(dates,vi,ylim=c(-0.2,0.8))


    dormIms <- dates >= pheno_pars$dormStart & dates <= pheno_pars$dormEnd
    vi_dorm1 <- quantile(vi[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T)   #Calc vi dormant value

    vi[vi<0]   <- NA
    vi_dorm <- quantile(vi[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T)   #Calc vi dormant value

    vi[vi < vi_dorm] <- vi_dorm

    #
    numDaysFit  <-  365 + (pheno_pars$splineBuffer * 2)
    daysVec <- 1:numDaysFit
    inYear <- daysVec > pheno_pars$splineBuffer & daysVec <= (pheno_pars$splineBuffer+365)



  yToDo <- 4


      splineStart <- as.Date(as.Date(paste0(phenYrs[yToDo],'-01-01')) - pheno_pars$splineBuffer)
      splineEnd <- as.Date(as.Date(paste0(phenYrs[yToDo],'-12-31')) + pheno_pars$splineBuffer)
      pred_dates <- seq(splineStart, splineEnd, by="day")

      dateRange <- dates >= splineStart & dates <= splineEnd & !is.na(vi)

      dateSub    <- dates[dateRange]
      viSub      <- vi[dateRange]
      weightsSub <- weights[dateRange]

      inyear <- as.numeric(format(dateSub,'%Y')) == phenYrs[yToDo]
      viObs_inyear  <- viSub[inyear]
      numObs <- sum(!is.na(viObs_inyear))


      smoothed_vi <- Smooth_VI(viSub, dateSub, pred_dates, weightsSub, pheno_pars, vi_dorm)




  points(pred_dates,smoothed_vi)
  abline(h=ndvi_dorm,col='blue',lwd=2)
  abline(h=vi_dorm1,col='red',lwd=2)
  abline(h=vi_dorm)
  text(18090,0.85,round(ndvi_dorm,3),pos=4,col='blue',cex=2)
  text(18090,0.75,round(vi_dorm1,3),pos=4,col='red',cex=2)
  text(18090,0.65,round(vi_dorm,3),pos=4,cex=1.5)
  text(18090,0.55,round(max(smoothed_vi)-min(smoothed_vi),3),pos=4,cex=1.5)

  dev.off()
  #
  # points(dateSub,viSub,col='red')
  # points(pred_dates,smoothed_vi,col='red')
}
      