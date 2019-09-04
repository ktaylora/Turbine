#!/usr/bin/Rscript

#
# Tubrine Tabular Asset Ingestion for PostgreSQL/PostGIS
#
# Author : Kyle Taylor (2019)
# License : GPL v.3 (2019)
#

#
# LOCAL FUNCTIONS
#

debug <- function(x) if(DEBUGGING) cat("DEBUG:", paste(x, collapse=" "), "\n")

est_runtime <- function(expression){
    elapsed = 3;
    return(as.vector(summary(system.time(
        expression
    )))[elapsed])
}

# slowly (but efficiently) estimate the number of turbines that intersect
# each US National Grid unit
counts <- sapply(
  1:nrow(units),
  FUN=function(i){
    return( sum(rowSums(!is.na(sp::over(units[i,], turbines, returnList=T)[[1]])) > 0) )
  }
)
# stage for attributing our topographic variables
topographic_vars <- raster::stack(
    "../raster/pljv_topographic_variables.tif"
  ) 
names(topographic_vars) <- readLines(
    "../raster/pljv_topographic_variables.txt"
  )
est_runtime(topographic_vars <- raster::extract(
    topographic_vars, 
    units, 
    fun=mean, 
    na.rm=T, 
    df=T, 
    progress='text'
  ))
