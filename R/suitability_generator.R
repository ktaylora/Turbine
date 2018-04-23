# Title     : 'Turbine' Suitability Generator
# Objective : Use a user-specified fitted model to generate a continuous raster surface of wind suitability
# Created by: Kyle Taylor (kyle.taylor@pljv.org)
# Created on: 4/19/18

#' hidden function that will perform a min-max normalization on an input raster to ensure it is projected as 0-to-1
min_max_normalize <- function(d) (exp(d) - min(exp(d)))/(max(exp(d)) - min(exp(d)))
#' generate a spatially-consistent stack of explanatory variables and cache the results to disk for
#' redistribution
explanatory_vars_to_rdata_file <- function(rasters=NULL, filename=NULL, n=NULL){
  rasters <- try(raster::stack(rasters))
  # if we can't make a raster stack from our input data, raise an error
  if(class(rasters) == "try-error"){
    stop("couldn't make a raster stack out of our input rasters -- are they spatially consistent?")
  }
  # did the user specify names for our raster stack?
  if(!is.null(n)) {
    names(rasters) <- n
  } else {
    n <- names(rasters)
  }
  # set a default filename for our r data file
  if(is.null(filename)){
    filename="explanatory_vars.rdata"
  }
  # cache a copy of our stack to disk as a multi-band GeoTIFF
  raster::writeRaster(
      rasters,
      filename="explanatory_vars.tif",
      options="INTERLEAVE=BAND",
      overwrite=TRUE,
      progress='text'
  )
  # re-read our input data from the CWD and generate an rdata file for export
  final_stack <- raster::stack("./explanatory_vars.tif")
  if(file.exists(filename)){
    r_data_envir <- new.env()
    load(filename, envir=r_data_envir)
    r_data_envir$final_stack <- final_stack
    r_data_envir$n <- n
    save(list=c("final_stack", "n"), file=filename, compress=T, compression_level=9, envir=r_data_envir)
  } else {
    save(list=c("final_stack", "n"), file=filename, compress=T, compression_level=9)
  }
}
#' function that will merge presences and absences SpatialPoints data.frame using a 'year' attribute
#' @export
gen_suitability_raster <- function(m=NULL, explanatory_vars=NULL, write=NULL, quietly=F, normalize=T){
  # produce an integer-scaled raster surface
  predicted <- round( raster::predict(
    object=explanatory_vars,
    model=m,
    type="response",
    na.rm=T,
    inf.rm=T,
    datatype='FLT4S',
    progress=ifelse(quietly==F, 'text', NULL)
  ) * 100 )
  # optionally re-scale a raster from 0-to-1
  if(normalize) predicted <- round( min_max_normalize(predicted) * 100 )
  # write to disk?
  if(!is.null(write)){
    raster::writeRaster(
      x=predicted,
      filename=write,
      format="GTiff",
      datatype="INT1S",
      overwrite=T
    )
  }
  # return to user
  return(predicted)
}
