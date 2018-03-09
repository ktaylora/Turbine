#!/usr/bin/env R
#
# __author__ = "Joe Hack, Kyle Taylor"
# __copyright__ = "Copyright 2018, Playa Lakes Joint Venture"
# __credits__ = ["Joe Hack", "Kyle Taylor"]
# __license__ = "GPL"
# __version__ = "3"
# __maintainer__ = "Kyle Taylor"
# __email__ = "kyle.taylor@pljv.org"
# __status__ = "Testing"

#
# R session options
#

options(warn = -1, error=traceback)

#
# Runtime arguments
#

N_BAGGING <- 999 # how many times should we re-sample our input datasets?
     argv <- commandArgs(trailingOnly=T)

#
# Main
#

require(Turbine)

# 1. read wind turbine points
wind_pts <- OpenIMBCR:::readOGRfromPath(
    "/gis_data/Wind/faa_products/wind_turbines_feb_26_2018.shp"
  )

# mask all points outside of the PLJV region
boundary <- sp::spTransform(
      OpenIMBCR:::readOGRfromPath(
        "/gis_data/PLJV/PLJV_Boundary.shp"
      ),
      sp::CRS(raster::projection(wind_pts))
    )

overlapping <- is.na(sp::over(wind_pts, boundary)[,1])
wind_pts <- wind_pts[!overlapping,]

# 2. Generate Pseudo-absences
pseudo_abs_pool <- Turbine:::gen_pseudo_absences(
    wind_pts,
    boundary=boundary
  )

# trend extrapolation of regional wind build-out -- how many turbines have
# gone-up over the span of our dataset?

build_out <- data.frame(table(wind_pts$year))
  colnames(build_out) <- c("year", "n_built")

build_out$year <- as.numeric(as.vector(build_out$year))
# let's drop the current year from consideration -- sample size is off
# because it isn't over yet.
build_out <- build_out[1:(nrow(build_out)-1),]

plot(log(n_built) ~ year, data=build_out, type="l")
  grid(); grid();
  points(log(n_built) ~ year, data=build_out, pch=15)

# EDA 1. ) Does agricultural production in areas around wind farms change after
# wind goes up? How so?

# EDA 2.) Given energy production capacity (MWh), how much revenue has this
# generated?

# EDA 3.) How does energy production capacity for wind compare to traditional
# oil and gas for the region? How much more build-out would have to occur for
# wind to provide an adequate substitution for oil and gas?

# EDA 4.) Has the build-out of wind been spatially correllated with any bird species declines
# over the past 10-20 years?

#
# Let's do some staging for our k-folds cross validation by year
#
# we have 17 years of data, but we should drop the first year (1998) because
# build-out was way lower than observed in most years. That leaves us with four
# chunks to train/test with from 1999:2017 :
#
# (5[train] + 6[test]) + (3[train] + 3[test]) = 17

build_out <- build_out[2:nrow(build_out),] # drop 1998 -- it's crap

training_chunks <- testing_chunks <- list()

# From 1999:2004, let's bag (sample with replacement) our input dataset
cat(" -- rebagging timeseries (1999:2004):")

training_chunks[[1]] <-
  lapply(
    X=1:N_BAGGING,
    FUN=function(x){
      cat(".")
      return(
        merge_presences_absences_by_year(
          wind_pts,
          pseudo_abs_pool,
          years=1999:2004
          )
      )
    }
  )
cat("\n")

cat(" -- building evaluation dataset [1] (2005:2010)\n")

testing_chunks[[1]] <-
  merge_presences_absences_by_year(
      wind_pts,
      pseudo_abs_pool,
      years=2005:2010,
      bag=F
    )
# From 2010:2012, let's bag (sample with replacement) our input dataset

cat(" -- rebagging timeseries (2011:2014):")

training_chunks[[2]] <-
  lapply(
    X=1:N_BAGGING,
    FUN=function(x){
      cat(".")
      return(
        merge_presences_absences_by_year(
          wind_pts,
          pseudo_abs_pool,
          years=2011:2014
          )
      )
    }
  )
cat("\n")

cat(" -- building evaluation dataset [2] (2015:2017)\n")

testing_chunks[[2]] <-
  merge_presences_absences_by_year(
      wind_pts,
      pseudo_abs_pool,
      years=2015:2017,
      bag=F
    )

# 3.1 Extract topographic predictor variables:
topographic_variables <-
  c('aspect.tif',
    'slope.tif',
    'topographic_roughness_3x3.tif',
    'topographic_roughness_11x11.tif',
    'topographic_roughness_33x33.tif'
  )
stack <- stack(mapply(raster,topographic_variables))
# had to create two different stacks because topographic variables have different projections
# and extents than the other variables.
transmission_variables <-
  c("raster_69",
    "raster_115",
    "raster_138",
    "raster_230",
    "raster_345"
  )
ws <- raster("ws_30m_wgs")
raster_115 <- raster("raster_115")
ws <- snapTo(ws, to = raster_115)
big_stack <- stack(mapply(raster,transmission_variables))
big_stack <-stack(big_stack,ws)
names(big_stack)

#
# 3.2 Combine all bags and extract explanatory variables
#
chunks_train_test <- training_chunks
chunks_train_test[[1]] <- do.call(rbind, training_chunks[[1]])
f <- chunks_train_test[[1]]
t<- extract(stack, f, df = T, progress= 'text')

# 4. Run the model!:

# 4.1: variable selection:
#

t <- na.omit(t)
model_sel <- rf.modelSel(ydata=as.factor(t$response),
                         xdata=t[,!grepl(names(t),pattern="response|ID")],
                         parsimony=0.15,
                         imp.scale="se",
                         r=seq(0,1,0.05),
                         seed=25,
                         final.model=T)

varImpPlot(model_sel$rf.final)

snapTo <- function(x,to=NULL,names=NULL,method='bilinear'){
  require(parallel)
  # set-up a cluster for parallelization
  cl <- makeCluster((parallel::detectCores()-2))
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
# There was a slight difference in extent between ws_30m_wgs and the transmission rasters.
# snapTo interpolates a raster into an extent specified by another raster.
# The code below writes all relevant rasters to ".tif" files and creates a rasterstack of
# these tifs.
ws <- raster("ws_30m.tif")
ras_names <- c("raster_138", "raster_230", "raster_115", "raster_69")
stack_list <- list()
for (i in length(ras_names)){
  ras <- raster(ras_names[i])
  tif <- writeRaster(ras, paste(ras_names[i],".tif",sep = ""), progress = "text")
  stack_list[i] <- tif
}
big_stack <- stack(mapply(raster,ras_names))
big_stack <- stack(big_stack,ws)
names(big_stack) = c("raster_138","raster_230","raster_115","raster_69","wind_speed_100m")

# predict! finally
test = raster::predict(big_stack, model_sel$rf.final, index = 2, type = "prob", progress= "text")
writeRaster(test, "suitability_prob.tif", progress = "text")
