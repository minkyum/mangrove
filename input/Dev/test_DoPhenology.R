library(raster)
library(rgdal)
library(gdalUtils)

library(rjson)

numSite <- 36; cc <- 50



########################################
params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/PlanetLSP/data_paper/PLSP_Parameters.json')
source(params$setup$rFunctions)

########################################
## Get site name, image directory and coordinate
strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

ckDir <- paste0(params$setup$outDir,strSite,'/chunk')
ckDir <- paste0('/projectnb/modislc/users/mkmoon/Planet/rawImage/chunks/',strSite)
print(ckDir)

ckNum <- sprintf('%03d',cc)
file <- list.files(path=ckDir,pattern=glob2rx(paste0('*',ckNum,'.rda')),full.names=T)

load(file)


pp <- 31

blue  <- band1[pp,]
green <- band2[pp,]
red   <- band3[pp,]
nir   <- band4[pp,]
phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr

## Load water mask
waterRater <- raster(paste0(params$setup$outDir,'Allequash_Creek_Site/water_mask_30.tif'))

cc <- 50
numCk <- params$setup$numChunks
chunk <- length(waterRater)%/%numCk
if(cc==numCk){
  chunks <- c((chunk*(cc-1)+1):length(waterRater))
}else{
  chunks <- c((chunk*(cc-1)+1):(chunk*cc))
}
waterMask <- values(waterRater)[chunks][pp]




########################################
DoPhenologyPlanet <- function(blue, green, red, nir, dates, phenYrs, params, waterMask){
  
  # Despike, calculate dormant value, fill negative VI values with dormant value
  log <- try({    
    
    pheno_pars <- params$phenology_parameters
    qa_pars    <- params$qa_parameters
    
    blue <- blue/10000; green <- green/10000; red <- red/10000; nir <- nir/10000
    vi   <- 2.5*(nir - red) / (nir + 2.4*red + 1)
    
    
    # Potential water
    if( sum(vi<0,na.rm=T)/sum(!is.na(vi)) > pheno_pars$sumNegVIthresh & waterMask > pheno_pars$waterOccuThreshLow ){
      return(rep(c(NA,rep(NA,10),4,rep(NA,10),4,NA),length(phenYrs)))} 
    
    
    # Spikes check, and remove
    spikes     <- CheckSpike_MultiBand(blue, red, vi, dates, pheno_pars)
    vi[spikes] <- NA
    
    # Replace negative VIs with dormant value
    dormIms <- dates >= pheno_pars$dormStart & dates <= pheno_pars$dormEnd
    vi_dorm <- quantile(vi[dormIms & vi>0],probs=pheno_pars$dormantQuantile,na.rm=T)   # Calc vi dormant value using non-negative VIs
    vi[vi < vi_dorm] <- vi_dorm
    
    #
    splineStart <- as.Date(as.Date(paste0(phenYrs,'-01-01')) - pheno_pars$splineBuffer) 
    numDaysFit  <- 365 + (pheno_pars$splineBuffer * 2)    
    splineEnd   <- splineStart+(numDaysFit-1)
    all_dates   <- seq(min(splineStart), max(splineEnd), by="day")
    
    daysVec <- 1:numDaysFit
    inYear  <- daysVec > pheno_pars$splineBuffer & daysVec <= (pheno_pars$splineBuffer+365)
    prevYear <- daysVec <= pheno_pars$splineBuffer
    nextYear <- daysVec > (pheno_pars$splineBuffer+365)
    
    
    ## 3-day composite
    if(pheno_pars$do3dayComposites){
      dates3c <- c(); vi3c <- c()
      
      dateNew <- 1
      for(dateSeq in seq(1,length(all_dates),3)){
        ind <- which(dates==all_dates[dateSeq]|dates==all_dates[dateSeq+1]|dates==all_dates[dateSeq+2])
        if(length(ind)>0 & sum(!is.na(vi[ind]))>0){
          dates3c[dateNew] <- all_dates[dateSeq+1]
          vi3c[dateNew]   <- max(vi[ind],na.rm=T)
          
          dateNew <- dateNew + 1
        }
      }
      dates    <- as.Date(dates3c,origin='1970-1-1')
      vi      <- vi3c
    }
    
    
    
    ## Gap filling
    #Determine gaps that require filling
    gDates <- dates[!is.na(vi)]  
    dDiff <- diff(gDates) > pheno_pars$gapLengthToFill     #Gaps greater than 20 days will be filled
    dStart <- gDates[c(dDiff,FALSE)]
    dEnd <- gDates[c(FALSE,dDiff)]
    
    #Locate gaps in date vector
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
    numYrs <- length(phenYrs)
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
        maskMat[format(pred_dates,'%Y') != phenYrs[y],y]  <- 1
        
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
    
    
  },silent=TRUE)
  #If there is an error despiking or other initial steps, return NAs
  if(inherits(log, "try-error")){return(rep(c(NA,rep(NA,10),4,rep(NA,10),4,NA),length(phenYrs)))} 
  
  
  outAll <- c()
  for(y in yToDo){
    log <- try({    
      
      pred_dates <- seq(splineStart[y], splineEnd[y], by="day")
      
      
      if (phenYrs[y] %in% yrsWithGaps) {
        
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
        
        # #Variables needed for next steps if the above gap filling was not done
        # theInds <- matrix(FALSE,length(ysGood))
        # theInds[((y-1)*numDaysFit+1):(y*numDaysFit)] <- TRUE
        # xs_sub <- xs[theInds]; w_sub <- baseW[theInds]
        
        smoothed_vi <- smoothMat[,y]   #if no gaps to fill, just use existing spline
      }
      
      
      # Number of clear observation
      filled_vi <- fillMat[,y]
      filled_vi[baseWeights[,y] < 1] <- NA    
      numObs <- sum(!is.na(filled_vi) & inYear)   #Number of observations in year
      
      #
      viSub   <- filled_vi[!is.na(filled_vi)]
      dateSub <- pred_dates[!is.na(filled_vi)]
      
      
      
      ################################################
      #Fit phenology
      peaks <- FindPeaks(smoothed_vi)
      if (all(is.na(peaks))) {outAll <- c(outAll,annualMetrics(viSub,dateSub,smoothed_vi,pred_dates,phenYrs[y],pheno_pars,vi_dorm,waterMask));next}
      
      #Find full segments
      full_segs <- GetSegs(peaks, smoothed_vi, pheno_pars)
      if (is.null(full_segs)) {outAll <- c(outAll,annualMetrics(viSub,dateSub,smoothed_vi,pred_dates,phenYrs[y],pheno_pars,vi_dorm,waterMask));next}
      
      #Only keep segments with peaks within year *****
      full_segs <- full_segs[inYear[sapply(full_segs, "[[", 2)] ]  #check if peaks are in the year
      if (length(full_segs)==0) {outAll <- c(outAll,annualMetrics(viSub,dateSub,smoothed_vi,pred_dates,phenYrs[y],pheno_pars,vi_dorm,waterMask));next}
      
      #Get PhenoDates
      pheno_dates <- GetPhenoDates(full_segs, smoothed_vi, pred_dates, pheno_pars)
      phen <- unlist(pheno_dates, use.names=F)
      phen <- phen - as.numeric(as.Date(paste0((as.numeric(phenYrs[y])-1),'-12-31')))
      if (all(is.na(phen))) {outAll <- c(outAll,annualMetrics(viSub,dateSub,smoothed_vi,pred_dates,phenYrs[y],pheno_pars,vi_dorm,waterMask));next}
      
      #EVI layers
      seg_metrics <- lapply(full_segs, GetSegMetrics,smoothed_vi,viSub,pred_dates,dateSub) #full segment metrics
      un <- unlist(seg_metrics, use.names=F)
      ln <- length(un)
      seg_amp <- un[seq(1, ln, by=9)] * 10000
      seg_max <- un[seq(2, ln, by=9)] * 10000
      seg_int <- un[seq(3, ln, by=9)] * 100  
      gup_rsq <- un[seq(4, ln, by=9)] * 10000
      gup_maxgap <- un[seq(6, ln, by=9)] 
      gdown_rsq <- un[seq(7, ln, by=9)] * 10000
      gdown_maxgap <- un[seq(9, ln, by=9)]
      
      
      ##
      theOrd <- order(seg_amp,decreasing=T)   
      
      # Filter for bad EVI layers
      if(seg_max[theOrd[1]] > 10000 | seg_max[theOrd[1]] < 0 | seg_amp[theOrd[1]] > 10000 | seg_amp[theOrd[1]] < 0){
        outAll <- c(outAll,c(NA,rep(NA,10),4,rep(NA,10),4,NA));next}
      
      # Filter for potential water
      if(vi_dorm < pheno_pars$VIdormThresh & seg_amp[theOrd[1]] < (pheno_pars$VIampThreshHigh*10000) & waterMask > pheno_pars$waterOccuThreshHigh){
        outAll <- c(outAll,c(NA,rep(NA,10),4,rep(NA,10),4,NA));next}
      
      if(vi_dorm < pheno_pars$VIdormThresh & seg_amp[theOrd[1]] < (pheno_pars$VIampThreshLow*10000) & waterMask > pheno_pars$waterOccuThreshLow){
        outAll <- c(outAll,c(NA,rep(NA,10),4,rep(NA,10),4,NA));next}
      
      
      numRecords <- length(seg_amp)  #how many cycles were recorded
      naCheck <- is.na(seg_amp)
      numCyc <- sum(naCheck == 0)  #how many cycles have good data (seg metrics has valid observations)
      
      
      # QA
      if(length(full_segs)==1){
        qual_1 <- GetQAs(gup_rsq, gdown_rsq, gup_maxgap, gdown_maxgap, theOrd, qa_pars)
      }else{
        qual_1 <- GetQAs(gup_rsq, gdown_rsq, gup_maxgap, gdown_maxgap, theOrd, qa_pars)[[1]][1]
        qual_2 <- GetQAs(gup_rsq, gdown_rsq, gup_maxgap, gdown_maxgap, theOrd, qa_pars)[[2]][1]
      } 
      
      
      ################################################
      if(numCyc == 0){outAll <- c(outAll,annualMetrics(viSub,dateSu,smoothed_vi,pred_dates,phenYrs[y],pheno_pars,vi_dorm,waterMask));next}
      
      if(numRecords == 1) {
        out <- c(1,phen,seg_max,seg_amp,seg_int,qual_1,c(rep(NA,10),4),numObs)
      }else{
        phen1 <- phen[seq(theOrd[1], length(phen), by = numRecords)]
        phen2 <- phen[seq(theOrd[2], length(phen), by = numRecords)]
        if(naCheck[theOrd[2]]){
          out <- c(numCyc,phen1,seg_max[theOrd[1]],seg_amp[theOrd[1]],seg_int[theOrd[1]],qual_1,
                   c(rep(NA,10),4),numObs)  
        }else{
          out <- c(numCyc,phen1,seg_max[theOrd[1]],seg_amp[theOrd[1]],seg_int[theOrd[1]],qual_1,
                   phen2,seg_max[theOrd[2]],seg_amp[theOrd[2]],seg_int[theOrd[2]],qual_2,numObs)  
        }
      }
      
    },silent=TRUE) #End of the try block
    if(inherits(log, "try-error")){outAll <- c(outAll,NA,c(rep(NA,10),4),c(rep(NA,10),4),NA)
    }else{outAll <- c(outAll,out);remove(out)}
    
  }
  
  return(outAll)
}



