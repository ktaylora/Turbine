# 
# NON-PARAMETRIC WIND ENERGY SUITABILITY MODEL
#
# Author: Kyle Taylor (kyle.taylor@pljv.org), with inspiration from the usual suspects.
#
# This is a (functional) work in progress.  Will work on increasing portability, speed of record generation, and
# dealing with class imbalances soon.
#
require(randomForest)
require(rgdal)
require(raster)
require(spatstat)

HOME               <- Sys.getenv("HOME");
GPLCC_PILOT_REGION <- readOGR(paste(HOME,"/PLJV/boundaries/GPLCC Pilot/GPLCC_Pilot_Region/",sep=""),"GPLCC_pilot_region_boundary_aggregated",verbose=F)

# raster source layers (explanatory data)
explanatory_variables <- list.files("/Volumes/GREEN/oil_and_gas/input",pattern="tif$",full.names=T)
  explanatory_variables <- lapply(as.list(explanatory_variables),FUN=raster)
    names(explanatory_variables) <- unlist(lapply(explanatory_variables,FUN=names))

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
while(length(forests)<10){ # build an arbitrary number of forests
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
# combine our forests into a forest
m_rf <- do.call(randomForest::combine, forests)
save.image("rfWindModelProjectWorkspace.rdata")
