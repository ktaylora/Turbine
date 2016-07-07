#
# NON-PARAMETRIC WIND ENERGY SUITABILITY MODEL
#
# Author: Kyle Taylor (kyle.taylor@pljv.org), with inspiration from the usual suspects (http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0089210).
#
# This is a (functional) work in progress.  Will work on increasing portability, speed of record generation, and
# dealing with class imbalances soon.
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

.include('randomForest')
.include('rgdal')
.include('raster')
.include('FedData')
.include('spatstat')

#
# LOCAL FUNCTIONS
#

topographic_variables <-
c(
  "elevation.tif",
  "elevationStdDev3.tif",
  "aspect.tif",
  "rough27.tif",
  "rough3.tif",
  "slope.tif"
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
fetchTopographicData <- function(x,useLocal=FALSE){
  # parse list or individual raster object
  if(useLocal){
    return(lapply(topographic_variables,FUN=raster))
  }
  # calculate from a live DEM we fetch from NED
  if(!require(FedData)) stop("'FedData' package not available -- please install")
  # clean-up any lurking temp file space.  Sometimes get_net doesn't do this all the way.
  unlink("tmp/1",recursive=T,force=T)
    unlink("tmp/dem",recursive=T,force=T)
  x_ <- try(get_ned(template=x,res="1",label="dem",extraction.dir="tmp",raw.dir="tmp",force.redo=T))
    while(class(x_) == "try-error"){ x_ <- try(get_ned(template=x,res="1",label="dem",extraction.dir="tmp",raw.dir="tmp",force.redo=T)) }
      x <- x_; rm(x_)
  topo_output <- bulkCalculateTopographicVariables(x)
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

#
# MAIN
#

HOME <- Sys.getenv("HOME");

if(file.exists("gplcc_region_extent.shp")){
  REGION_EXTENT <- readOGR(HOME,"gplcc_region_extent",verbose=F)
} else {
  REGION_EXTENT <- file.choose()
  path <- unlist(strsplit(REGION_EXTENT,split="\\\\"))
    layer <- path[length(path)]
      dsn <- paste(path[1:(length(path)-1)],collapse="\\")
  REGIONAL_EXTENT <- readOGR(layer,dsn,verbose=F)
}

# parse our regional extent, if needed
if(!inherits(REGION_EXTENT,"Spatial")){
  stop("Error: Failed to parse a regional extent shapefile for our run.")
}

# raster source layers (explanatory data)
explanatory_variables <- list.files(recursive=T,pattern=paste(topographic_variables,collapse="|"))

# are we missing topographic data?
if(length(explanatory_variables)<length(topographic_variables)){
  e <- splitExtent(extent(REGION_EXTENT),multiple=20)
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
    writeRaster(focal_variable, paste(topographic_variables[j],"tif",sep="."), overwrite=T)
  }

}

# delineate regional boundaries using the extent of BCR18/19
e <- extent(raster(paste(HOME,"/PLJV/DEM/elevationBcr1819.tif",sep="")))*0.85
  e<-as(e,'SpatialPolygons')
    projection(e) <- projection(raster(paste(HOME,"/PLJV/DEM/elevationBcr1819.tif",sep="")))
# read-in and crop our wellpad data to our extent
source_wellpad_pts <- paste(HOME,"/PLJV/wind_energy_likelihood_model/data/USGSWind_Turbine_03_2014.zip",sep="")
  unlink("/tmp/wp_pts",recursive=T, force=T); utils::unzip(source_wellpad_pts,exdir="/tmp/wp_pts")
     source_wellpad_pts <- spTransform(readOGR("/tmp/wp_pts/","Onshore_Industrial_Wind_Turbine_Locations_for_the_United_States_to_March_2014",verbose=F), CRS(projection(e)))
       source_wellpad_pts <- source_wellpad_pts[!is.na(as.vector(sp::over(source_wellpad_pts,e))),]
# generate point samples and forests
forests <- list();
while(length(forests)<50){ # build an arbitrary number of forests
  # generate a series of N kernel density surfaces with a downsampled p=0.85 density
  cat(" -- generating absence surface from KDE.\n")
  focal <- SpatialPoints(source_wellpad_pts[sample(1:nrow(source_wellpad_pts),size=floor(0.85*nrow(source_wellpad_pts)),replace=T),])
  w <- extent(focal)*1.5;
    w <- owin(xrange=w[c(1,2)],yrange=w[c(3,4)])
      focal <- unique(suppressWarnings(ppp(x=focal@coords[,1], y=focal@coords[,2], window=w)))
        presences <-  as(focal, 'SpatialPoints') # keep a copy of these as our "presences"
          projection(presences) <- projection(source_wellpad_pts)
  d <- raster(density.ppp(focal,edge=T,sigma=1.05))
    d <- abs((d/cellStats(d,stat=max))-1)
      d <- sampleRandom(d,size=floor(ncell(d)*0.75),sp=T)
        projection(d) <- projection(source_wellpad_pts)
  # contribute this kernel to the absence pool, weighting 'distant' locations as more likely to represent absence than 'close' locations
  absences <- NULL;
  focal <- d[d$layer<1 & d$layer>0.85,]
    absences <- focal[sample(1:nrow(focal),size=round(nrow(focal)*0.4)),]
  focal <- d[d$layer<0.85 & d$layer>0.7,]
    absences <- rbind(absences,focal[sample(1:nrow(focal),size=round(nrow(focal)*0.3)),])
  focal <- d[d$layer<0.7 & d$layer>0.5,]
    absences <- rbind(absences,focal[sample(1:nrow(focal),size=round(nrow(focal)*0.1)),])
  focal <- d[d$layer<0.5 & d$layer>0.3,]
    absences <- rbind(absences,focal[sample(1:nrow(focal),size=round(nrow(focal)*0.1)),])
  focal <- d[d$layer<0.3 & d$layer>0.01,]
    absences <- rbind(absences,focal[sample(1:nrow(focal),size=round(nrow(focal)*0.1)),])
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
