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

# EDA 1. ) Does agricultural production in areas around wind farms change after 
# wind goes up? How so?  

# EDA 2.) Given energy production capacity (MWh), how much revenue has this 
# generated?

# EDA 3.) How does energy production capacity for wind compare to traditional 
# oil and gas for the region? How much more build-out would have to occur for 
# wind to provide an adequate substitution for oil and gas?

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

# From 1999:2003, let's bag (sample with replacement) our input dataset 
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
