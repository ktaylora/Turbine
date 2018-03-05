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
options(warn = -1, error=traceback)
argv <- commandArgs(trailingOnly=T)

bs_sample <- function(pts=NULL){
    return(NULL)
  }
  
#' generate a large pool of random pseudo-absences within a geographic boundary 
#' using a set of user-specified wind turbine point locations. The pool will be 
#' arbitrarily large (nrow(pts)*iter) with the intention that the dataset will
#' be downsampled through a bagging procedure down-the-line.
gen_pseudo_absences <- function(pts=NULL, boundary=NULL, iter=100, buffer_width=1000){
  # let's reproject our data to a planar geography for rgeos
  original_crs <- sp::CRS(raster::projection(pts))
  pts <- sp::spTransform(
      pts, 
      sp::CRS(raster::projection("+init=epsg:2163"))
    )
  boundary <- sp::spTransform(
      boundary, 
      sp::CRS(raster::projection("+init=epsg:2163"))
    )
  # buffer our turbine occurrences by some user-specified distance
  buffered_turbine_occurrences <- rgeos::gBuffer(
      pts, byid=T, width=buffer_width
    )
  cat(" -- generating a large pool of pseudo-absences\n")
  wind_pts_raster <- raster::raster(
      resolution=30,
      ext=raster::extent(boundary),
      crs=sp::CRS(raster::projection(boundary))
    )
  # set our pool to zero  
  raster::values(wind_pts_raster) <- 0  
  # randomly generate a pool of pseudo-absences across n=iter steps
  cat(" -- iteritively re-sampling absence space:")
  pseudo_absences <- do.call(sp::rbind.SpatialPointsDataFrame, lapply(
    X=1:iter,
    FUN=function(r){
      cat(".")
      absences <- raster::sampleRandom(
          wind_pts_raster,
          size=round(nrow(wind_pts)*1.2),
          sp=T,
          na.rm=F
        )
      absences <- absences[
        is.na(sp::over(absences, buffered_turbine_occurrences)[,1]),]
      return(absences)
    }
  ))
  cat("\n")
  pseudo_absences@data <- data.frame(response=rep(0, nrow(pseudo_absences)))
  # reproject back to our native CRS and return to user
  return(sp::spTransform(pseudo_absences, original_crs))
}
select_records_by_years <- function(pts=NULL, years=NULL){

  }

#
# Main
#
  
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
pseudo_abs_pool <- gen_pseudo_absences(
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

# Does agricultural production in areas around wind farms change after wind goes up? 
# How so?  

# Given energy production capacity (MWh), how much revenue has this generated?

# How does energy production capacity for wind compare to traditional oil and gas
# for the region? How much more build-out would have to occur for wind to provide
# an adequate substitution for oil and gas?


# 3.1 Extract topographic predictor variables:
topographic_variables <-
  c('aspect.tif',
    'slope.tif',
    'topographic_roughness_3x3.tif',
    'topographic_roughness_11x11.tif',
    'topographic_roughness_33x33.tif'
  )
# line below creates a raster stack of all topographic rasters
stack <- stack(mapply(raster,topographic_variables))

t <- extract(stack,training_data,df=T)
# extract will extract datapoints from each raster in the stack for all ~ 32,000
# training points. How to validate that points were drawn from the same spatial location?
t$response <- training_data@data[,1]
# for some reason response was lost during the extraction process and had to be added back in.

# 3.2 Extract wind speed variables:
# Right now I only have data for windspeed at 100m. It looks like the average turbine height
# is around 300 ft, right between 100 and 80m, which NREL also has a map for. We might want to
# interpolate an 80m raster too, although it does look like the two maps are highly correllated.
setwd("/home/jhack/global_workspace/Joe_Work/turbine_workspace")
wind_speed_100m <- raster("ws_30m")
# Extract wind_speed_100m points. I will rewrite extraction steps later so that all predictor 
# variables in the topographic, windspeed and transmission categories are loaded at once.
f <- extract(wind_speed_100m, training_data, df=T)
# Merge f dataframe with t
t$wind_speed_100m <- f[,2]

# 3.3 Extract transmission data variables:
transmission_variables <-
  c("raster_69",
    "raster_115",
    "raster_138",
    "raster_230",
    "raster_345"
    )
q <- extract(stack(mapply(raster,transmission_variables)),training_data,df=T)
t$raster_69 <- q[,2]
t$raster_115 <- q[,3]
t$raster_138 <- q[,4]
t$raster_230 <- q[,5]
t$raster_345 <- q[,6]

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

# the command above will run the final model, but I believe Kyle wants to build many different forests.
# to resolve: cross-validation.
