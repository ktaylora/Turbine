#
# NON-PARAMETRIC WIND ENERGY SUITABILITY MODEL
#
# Author: Kyle Taylor (kyle.taylor@pljv.org), with inspiration from the usual suspects (http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0089210).
#
# This is a (functional) work in progress.  Will work on increasing portability, speed of record generation, and
# dealing with class imbalances soon.
#
# Target platform: A multi-core 64-bit windows machine with ample ram to support large raster operations.
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
.include('spatstat')
.include('randomForest')

topographic_variables <-
c(
  "elevation.tif",
  "elevationStdDev3.tif",
  "aspect.tif",
  "rough27.tif",
  "rough3.tif",
  "slope.tif"
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
#
# splitExtent()
#
# When working with SDB queries for large areas, we consistently lose a lot of data... particularly for large
# counties.  This gets around that by splitting an extent object into adjacent quarters so that we can download and
# merge our raster segments later.  For small counties, this is inefficient.  But its better than just flatly attempting downloads
#
# Author: Kyle Taylor (kyle.taylor@pljv.org) [2016]
#
splitExtent <- function(e=NULL,multiple=2){
  .include('raster')
  # define our x/y vector ranges
  x <- rep(NA,multiple+1)
  y <- rep(NA,multiple+1)
  # define the x/y range for calculating the size of our extents
  xStep <- diff(c(e@xmin,e@xmax))/multiple
  yStep <- diff(c(e@ymin,e@ymax))/multiple
  # assign vertices to our product vectors
  for(i in 1:(multiple+1)){
    x[i] <- ifelse(i==1,
                   min(e@xmin),
                   x[i-1]+xStep)
    y[i] <- ifelse(i==1,
                   min(e@ymin),
                   y[i-1]+yStep)
  }
  # assign our vertices to extent objects
  extents <- as.list(rep(NA,multiple*multiple))
  # iterate over our extents, assigning as we go
  yStart <- i <- 1;
  while(i <= length(extents)){
    for(j in 1:multiple){ # stagger our y-values
      extents[i] <- extent(c(x[j],x[j+1],y[yStart],y[yStart+1]))
      i <- i+1;
    }
    yStart <- yStart+1;
  }
  return(extents)
}
#
# fetchTopographicData()
# wrapper function for FedData::get_ned()
#
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
# bulkCalculateTopographicVariables()
#
bulkCalculateTopographicVariables <- function(x){
    o <- list()
    o[[1]] <- x # elevation
    o[[2]] <- raster::focal(x,w=matrix(1,nrow=3,ncol=3),fun=sd,na.rm=T) # StdDevElev (3x3)
    o[[3]] <- raster::terrain(x,opt='aspect',neighbors=8)
    o[[4]] <- raster::focal(x, w=matrix(1,nrow=27,ncol=27), fun=function(x, ...) max(x) - min(x),pad=TRUE, padValue=NA, na.rm=TRUE)
    o[[5]] <- raster::focal(x, w=matrix(1,nrow=3,ncol=3), fun=function(x, ...) max(x) - min(x),pad=TRUE, padValue=NA, na.rm=TRUE)
    o[[6]] <- raster::terrain(x,opt='slope',neighbors=8)
    return(o)
}
#
# snapTo()
# Ensure absolute consistency between raster objects by cropping,projecting,snapping,and (if asked) resampling
# a raster object using a template
#
snapTo <- function(x,to=NULL,names=NULL,method='bilinear'){
  require(parallel)
  # set-up a cluster for parallelization
  cl <- makeCluster((parallel::detectCores()-4))
  # crop, reproject, and snap our raster to a resolution and projection consistent with the rest our explanatory data
  if(grepl(tolower(class(x)),pattern="character")){ lapply(x,FUN=raster) }
  e <- as(extent(to[[1]]),'SpatialPolygons')
    projection(e) <- CRS(projection(to[[1]]))
  if(class(x) == "list") {
    x <- parLapply(cl,x,fun=raster::crop,extent(spTransform(e,CRS(projection(x[[1]])))))
      x <- parLapply(cl,x,fun=raster::projectRaster,crs=CRS(projection(to[[1]])))
    extents <- lapply(x,alignExtent,to[[1]])
      for(i in 1:length(x)){ extent(x[[i]]) <- extents[[i]] }
    if(!is.null(method)){
      x <- parLapply(cl,x,fun=resample,y=to[[1]],method=method)
    }
  } else {
    x <- raster::crop(x,extent(spTransform(e,CRS(projection(x)))))
      x <- raster::projectRaster(x,crs=CRS(projection(to[[1]])))
    extent <- alignExtent(x,to[[1]])
      extent(x) <- extent
    if(!is.null(method)){
      x <- raster::resample(x,y=to[[1]],method=method)
    }
  }
  endCluster()
  return(x)
}

pathToRgdal <- function(x){
  .include('rgdal')
  path <- unlist(strsplit(x,split="\\\\"))
  layer <- unlist(strsplit(path[length(path)],split="[.]"))[1]
  dsn <- paste(path[1:(length(path)-1)],collapse="\\")
  return(readOGR(dsn,layer,verbose=F))
}

#
# MAIN
#

HOME <- Sys.getenv("HOME"); setwd(HOME);

if(file.exists("gplcc_region_extent.shp")){
  REGION_EXTENT <- readOGR(HOME,"gplcc_region_extent",verbose=F)
} else {
  REGION_EXTENT <- pathToRgdal(file.choose())
}

# parse our regional extent, if needed
if(!inherits(REGION_EXTENT,"Spatial")){
  stop("Error: Failed to parse a regional extent shapefile for our run.")
}

# pre-processing : process our raster source layers and check for spatial consistency(explanatory data)

# process topographic variables
explanatory_variables <- list.files(recursive=T,pattern=paste(topographic_variables,collapse="|"), full.names=T)

# are we missing topographic data?
if(length(explanatory_variables) < length(topographic_variables)){
  cat(" -- missing topographic variables. Will attempt to fetch/calculate missing from USGS.\n")
  e <- splitExtent(extent(REGION_EXTENT),multiple=5)
  fragments <- vector('list', length(e))
  # fetch raw DEM data from USGS
  for(i in 1:length(fragments)){
    cat(paste("[",i,"/",length(fragments),"]\n",sep=""))
    template <- raster(extent(e[[i]]),resolution=30,crs=CRS(projection(REGION_EXTENT)))
      fragments[[i]] <- fetchTopographicData(template) # fragments will be a list of lists with our six target topographic variables
  }
  # mosaic our fragments into each target topographic variables
  cat(" -- mosaicing topographic grid for:")
  for(j in 1:6){
    cat(topographic_variables[j],",")
    focal_variable <- lapply(fragments, function(x){ return(x[[j]]) })
      focal_variable <- do.call(mosaic,focal_variable)
    writeRaster(focal_variable, topographic_variables[j], overwrite=T)
  }
  explanatory_variables <- lapply(list.files(recursive=T,pattern=paste(topographic_variables,collapse="|"), full.names=T), FUN=raster)
} else {
  explanatory_variables <- lapply(explanatory_variables, FUN=raster)
}

# process transmission capacity data
if(length(list.files(recursive=T,full.names=T,pattern=paste(transmission_capacity_variables, collapse="|"))) == length(transmission_capacity_variables)){
  transmission_vars <- lapply(list.files(recursive=T,full.names=T,pattern=paste(transmission_capacity_variables, collapse="|")), FUN=raster)
  # ensure spatial consistency by snapping transmission vars to existing stack, then add to stack
  transmission_vars <- snapTo(transmission_vars,to=explanatory_variables[[1]])
    explanatory_variables <- append(explanatory_variables, transmission_vars)
} else {
  stop(paste("Error : Recursive find failed to find required transmission variables:",transmission_capacity_variables,sep=""))
}

# process wind production classes
if(length(list.files(recursive=T,full.names=T,pattern=wind_production_variables))>0){
  wind_production_vars <- lapply(list.files(recursive=T,full.names=T,pattern=wind_production_variables), FUN=raster)
    wind_production_vars <- snapTo(wind_production_vars,to=explanatory_variables[[1]])
      explanatory_variables <- append(explanatory_variables, wind_production_vars)
}
# process our turbine SpatialPoint data
if(file.exists("turbine_points.shp")){
  turbines <- readOGR(HOME,"turbine_points",verbose=F)
} else {
  cat(" -- please choose a shapefile containing our turbine point data")
  turbines <- pathToRgdal(file.choose())
}

# parse our regional extent, if needed
if(!inherits(turbines,"Spatial")){
  stop("Error: Failed to parse our turbine point shapefile for our run.")
}

# pre-process our points to ensure we only include turbines within the extent of our project area
turbines <- spTransform(turbines,CRS(projection(REGION_EXTENT)))
  turbines <- turbines[as.vector(!is.na(sp::over(turbines,REGION_EXTENT))),]

# generate point samples and forests
forests <- list();
while(length(forests)<100){ # build an arbitrary number of forests
  # generate a series of N kernel density surfaces with a downsampled p=0.85 density
  cat(" -- generating absence surface from KDE.\n")
  focal <- SpatialPoints(turbines[sample(1:nrow(turbines),size=floor(0.5*nrow(turbines)),replace=T),])
  w <- extent(focal)*1.5;
    w <- owin(xrange=w[c(1,2)],yrange=w[c(3,4)])
      focal <- unique(suppressWarnings(ppp(x=focal@coords[,1], y=focal@coords[,2], window=w)))
        presences <-  as(focal, 'SpatialPoints') # keep a copy of these as our "presences"
          projection(presences) <- projection(turbines)
  d <- raster(density.ppp(focal,edge=T,sigma=1.05))
    d <- round(d/cellStats(d,stat=max)*100)

  # contribute this kernel to the absence pool, weighting 'distant' locations as more likely to represent absence than 'close' locations
  SampleStratified <- function(d,count=500,split=70){
    s <- suppressWarnings(sampleStratified(d<split,size=c(count,1),sp=T))
      return(s[s$layer == 1,])
  }
  absences <- SampleStratified(d,count=900,split=50)
  for(i in 1:4){
    absences <- rbind(absences,SampleStratified(d,count=900,split=50+(i*10)))
  }

  # extract our presence/absence data
  if(length(presences)>nrow(absences)){
    presences <- presences[sample(1:length(presences),size=nrow(absences)),]
  } else {
    absences <- absences[sample(1:nrow(absences),size=length(presences)),]
  }
  cat(" -- extracting ...\n")
  t1<-Sys.time()
  presences <- data.frame(lapply(explanatory_variables, FUN=extract, y=presences))
    presences <- na.omit(cbind(resp=1,presences))
  absences  <- data.frame(lapply(explanatory_variables,FUN=extract,y=absences))
    absences <- na.omit(cbind(resp=0,absences))
  # balance our classes
  if(nrow(presences)>nrow(absences)){
    presences <- presences[sample(1:nrow(presences), size=nrow(absences)),]
  } else {
    absences <- absences[sample(1:nrow(absences), size=nrow(presences)),]
  }
  t <- rbind(presences,absences)
  # downscale to a consistent record density, sampling with replacement if necessary.
  if(nrow(t)<1500) {
    t <- t[sample(1:nrow(t), size=1500, replace=T),]
    warning("resampling extracted record table with replacement because n < 1500 -- this could lead to class imbalances")
  } else {
    t <- t[sample(1:nrow(t), size=1500, replace=F),]
  }
  # report our ETA
  t2<- Sys.time()
  cat(" -- eta: ~",as.numeric(t2-t1)*(100-(length(forests)+1)),"minutes remaining\n")
  # build our random forest
  forests[[length(forests)+1]] <- randomForest(as.factor(resp)~.,data=t,ntree=100,norm.votes=F,do.trace=T)
};cat(" -- done.\n");
# combine our forests into a single forest we can predict with
m_rf <- do.call(randomForest::combine, forests)
save.image("rfWindModelProjectWorkspace.rdata")
