# Title     : 'Turbine' Suitability Generator
# Objective : Use a user-specified fitted model to generate a continuous raster surface of wind suitability
# Created by: Kyle Taylor (kyle.taylor@pljv.org)
# Created on: 4/19/18

#' hidden function that will perform a min-max normalization on an input raster to ensure it is projected as 0-to-1
min_max_normalize <- function(d) (exp(d) - min(exp(d)))/(max(exp(d)) - min(exp(d)))
#' generate a spatially-consistent stack of explanatory variables and cache the results to disk for
#' redistribution
explanatory_vars_to_rdata_file <- function(rasters=NULL, filename=NULL, names=NULL){
  rasters <- try(raster::stack(rasters))
  if(class(rasters) == "try-error"){
    stop("couldn't make a raster stack out of our input rasters -- are they spatially consistent?")
  }
  if(is.null(names)) {
    names <- names(rasters)
  } else {
    names(rasters) <- names
  }
  if(is.null(filename)){
    filename="explanatory_vars.rdata"
  }
  raster::writeRaster(
      rasters,
      filename="explanatory_vars.tif",
      options="INTERLEAVE=BAND",
      overwrite=TRUE,
      progress='text'
  )
  # re-read our input data and generate an rdata file
  final_stack <- raster::stack("explanatory_vars.tif")
  save(list="final_stack", file=filename, compress=T, compression_level=9)
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
