# Title     : 'Turbine' Suitability Generator
# Objective : Use a user-specified fitted model to generate a continuous raster surface of wind suitability
# Created by: Kyle Taylor (kyle.taylor@pljv.org)
# Created on: 4/19/18

#' hidden function that will perform a min-max normalization on an input raster to ensure it is projected as 0-to-1
min_max_normalize <- function(d) {
  d <- exp(d)
  min <- raster::cellStats(d, stat=min)
  max <- raster::cellStats(d, stat=max)
  return( ( d - min ) / ( max - min) )
}
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
    r_data_envir$n <- n
    save(list=c("n"), file=filename, compress=T, compression_level=9, envir=r_data_envir)
  } else {
    save(list=c("n"), file=filename, compress=T, compression_level=9)
  }
}
#' function that will merge presences and absences SpatialPoints data.frame using a 'year' attribute
#' @export
gen_suitability_raster <- function(m=NULL, explanatory_vars=NULL, write=NULL, quietly=F, normalize=T, n=NULL){
  # split-up a large raster into chunks that we can process
  # in parallel
  stopifnot(require(SpaDES))
  if(is.null(n)){
    n <- parallel::detectCores()-1
  }

  cl <- parallel::makeCluster(n)
  # hackish way of loading SpaDES package on our cluster
  stopifnot(sum(unlist(parallel::clusterApply(cl, x=rep(1,n), fun=function(x){require(SpaDES)}))) == n)
  # export our explanatory_vars
  parallel::clusterExport(cl, varlist="explanatory_vars")
  # produce an integer-scaled raster surface
  chunks <- parallel::parLapply(
    cl=cl,
    X=1:raster::nlayers(explanatory_vars),
    fun=function(band){
      band <- raster::subset(explanatory_vars, band)
      return(splitRaster(band))
  })
  # clean-up our cluster
  parallel::stopCluster(cl)
  rm(cl)

  band <- parallel::parLapply(
    cl=cl,
    X=chunks,
    fun=function(chunk){

      return(round( min_max_normalize( raster::predict(
        object=chunk,
        model=m,
        type="response",
        na.rm=T,
        inf.rm=T,
        progress=ifelse(quietly==F, 'text', NULL)
      ) ) * 100 ))
    }
  )
  predicted <- round( min_max_normalize( raster::predict(
    object=explanatory_vars,
    model=m,
    type="response",
    na.rm=T,
    inf.rm=T,
    progress=ifelse(quietly==F, 'text', NULL)
  ) ) * 100 )
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
