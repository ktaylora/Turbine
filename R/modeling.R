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
merge_presences_absences_by_year <- function(presences=NULL, absences=NULL, years=NULL, bag=T){
    presences  <- presences[presences$year %in% years, ]
    if(bag){
      presences <- presences[
          sample(1:nrow(presences), replace=T, size=nrow(presences)),
        ]
      absences <- absences[
        sample(1:nrow(absences), size=nrow(presences), replace=T),
      ]  
    } else {
      absences <- absences[
        sample(1:nrow(absences), size=nrow(presences), replace=F),
      ]
    }
    # merge and return to user
    presences@data <- data.frame(response=rep(1, nrow(presences)))
    absences@data <- data.frame(response=rep(0, nrow(absences)))
    return(sp::rbind.SpatialPointsDataFrame(presences,absences))
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

