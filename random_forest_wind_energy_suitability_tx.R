#
# WIND ENERGY DEVELOPMENT SUITABILITY MODEL
#
# Author: Kyle Taylor, Playa Lakes Joint Venture
# Date : 9/2016
#
# Latest iteration of the wind energy dev suitability model.  I have dropped the weighted-absence selection,
# because it didn't significantly improve (and in some cases, worsened) the accuracy of the model vs. a NULL
# model of random background pseudo-absences.  The new approach uses sensoring within a buffered distance
# consistent with the nearest-neighbor distance between turbines on the landscape (i.e., any location outside of a wind farm).
#

#
# include()
# wrapper function for require that will rudely attempt to install missing packages.
#
# Author: Kyle Taylor
#
.include <- function(x,from="cran",repo=NULL){
  if(from == "cran"){
    if(!do.call(require,as.list(x))) install.packages(x, repos=c("http://cran.revolutionanalytics.com","http://cran.us.r-project.org"));
      if(!do.call(require,as.list(x))) stop("auto installation of package ",x," failed.\n")
  } else if(from == "github"){
    if(!do.call(require,as.list(x))){
      if(!do.call(require,as.list('devtools'))) install.packages('devtools', repos=c("http://cran.revolutionanalytics.com","http://cran.us.r-project.org"));
      require('devtools');
      install_github(paste(repo,x,sep="/"));
    }
  } else{
    stop(paste("could find package:",x))
  }
}

.include('rgdal')
.include('raster')
.include('FedData')
.include('rgeos')
.include('randomForest')
.include('landscapeAnalysis')

topographic_variables <-
c(
  "elevation_tx.tif",
  "elevationStdDev3_tx.tif",
  "aspect_tx.tif",
  "rough27_tx.tif",
  "rough3_tx.tif",
  "slope_tx.tif"
)

# We lack average wind speed at 80 and 100 meters.
# and must use NREL's production classes for now
wind_production_variables <-
c(
  "lower_48_wind_production.tif",
  "us_wind_50m_no_exclusions.tif"
)
# proximity to transmission surfaces
transmission_capacity_variables <-
c(
  "transmission_capacity_step_up.tif",
  "transmission_capacity_100.tif",
  "transmission_capacity_100_161.tif",
  "transmission_capacity_230_287.tif",
  "transmission_capacity_345.tif",
  "transmission_capacity_500.tif"
)

recursiveFindFile <- function(name=NULL,root=Sys.getenv("HOME")){
  return(list.files(root,pattern=name,recursive=T,full.names=T))
}

fetchTopographicData <- function(x,dem=NULL){
  # parse list or individual raster object
  if(is.null(dem)){
    # calculate from a live DEM we fetch from NED
    if(!require(FedData)) stop("'fedData' package not available -- please install")
    # clean-up any lurking temp file space.  Sometimes get_net doesn't do this all the way.
    unlink("/tmp/1",recursive=T,force=T)
      unlink("/tmp/dem",recursive=T,force=T)
    x_ <- try(get_ned(template=x,res="1",label="dem",extraction.dir="/tmp",raw.dir="/tmp",force.redo=T))
      while(class(x_) == "try-error"){ x_ <- try(get_ned(template=x,res="1",label="dem",extraction.dir="/tmp",raw.dir="/tmp",force.redo=T)) }
        x <- x_; rm(x_)
  } else {
    x <- raster(dem)
  }
  topo_output <- list()
    topo_output[[1]] <- x # elevation
    topo_output[[2]] <- raster::focal(x,w=matrix(1,nrow=3,ncol=3),fun=sd,na.rm=T) # StdDevElev (3x3)
    topo_output[[3]] <- raster::terrain(x,opt='aspect',neighbors=8)
    topo_output[[4]] <- raster::focal(x, w=matrix(1,nrow=27,ncol=27), fun=function(x, ...) max(x) - min(x),pad=TRUE, padValue=NA, na.rm=TRUE)
    topo_output[[5]] <- raster::focal(x, w=matrix(1,nrow=3,ncol=3), fun=function(x, ...) max(x) - min(x),pad=TRUE, padValue=NA, na.rm=TRUE)
    topo_output[[6]] <- raster::terrain(x,opt='slope',neighbors=8)
  return(topo_output)
}

#
# MAIN
#

# fetch our turbine points and extent of the project area
s <- recursiveFindFile(root="/home/ktaylora/Incoming","tx.shp")
  extent <- landscapeAnalysis:::.readOGRfromPath(s[2])
    s <- landscapeAnalysis:::.readOGRfromPath(s[1])
      extent <- spTransform(extent,CRS(projection(s))) # snap to a consistent CRS

# generate some random absences and merge records into a training set for Random Forests
cat(" -- generating pseudo-absences and building a training data set for RF.\n")
sample_space <- rasterize(s,raster(resolution=(1/111319.9)*30,ext=extent(s),crs=CRS(projection(s))),field="layer") > 0
# Generate random, background pseudo-absences.  Censor all observations so that absences are at-least 343 meters away
# from a known wind turbine -- as identified by A. Daniels for K=1 NN distance
absences <- sampleRandom(sample_space,size=round(nrow(s)*1.2),sp=T,na.rm=F)
  absences <- absences[is.na(sp::over(absences,rgeos::gBuffer(s,byid=T,width=(1/111319.9)*343))[,1]),]
if(nrow(absences) > nrow(s)){
  absences <- absences[sample(1:nrow(absences), size=nrow(s)),] # subsample our absences to equal-prevalence with turbine locations so we have balanced classes
}

training_data <- rbind(s,absences)
  training_data@data <- data.frame(response=c(rep(1,nrow(s)),rep(0,nrow(absences))))

# calculate our topographic variables
if(length(recursiveFindFile(root="/home/ktaylora/Incoming",topographic_variables[1]))==0)){
  topographic_variables <- fetchTopographicData(sample_space)
}

# buffer Texas project area by 10 kilometers to account for boundary effect in calculating landscape metrics
extent <- rgeos::gBuffer(extent,byid=T,width=(1/111319.9)*10000)
