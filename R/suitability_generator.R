# Title     : 'Turbine' Suitability Generator
# Objective : Use a user-specified fitted model to generate a continuous raster
# surface of wind suitability
# Created by: Kyle Taylor (kyle.taylor@pljv.org)
# Created on: 4/19/18

#' hidden function that will perform a min-max normalization on an input raster to ensure it is projected as 0-to-1
min_max_normalize <- function(d) {
  d <- exp(d)
  min <- raster::cellStats(d, stat = min)
  max <- raster::cellStats(d, stat = max)
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
  if(file.exists(filename)){
    r_data_envir <- new.env()
    load(filename, envir=r_data_envir)
    r_data_envir$n <- n
    save(list=c("n"), file=filename, compress=T, compression_level=9, envir=r_data_envir)
  } else {
    save(list=c("n"), file=filename, compress=T, compression_level=9)
  }
}
#' function that will generate a suitability surface using a random forest model
#' @export
gen_rf_suitability_raster <- function(m=NULL, explanatory_vars=NULL, write=NULL){
  require(raster)
  require(randomForest)
  # raster::predict will return the probability of the "absence" (0) class
  # we are inverting the logic here so that "presence" is 1-absence
  p <- raster::predict(
        object=explanatory_vars,
        index=2,
        model=m,
        type='prob',
        datatype='FLT4S'
  )
  p <- round( p * 100 )
  # write to disk?
  if(!is.null(write)){
    raster::writeRaster(
      x=p,
      filename=write,
      format='GTiff',
      datatype='INT1U',
      overwrite=T
    )
    return(p)
  } else {
    return(p)
  }
}
#' function that will merge presences and absences SpatialPoints data.frame using a 'year' attribute
#' @export
gen_gam_suitability_raster <- function(
  m=NULL,
  explanatory_vars=NULL,
  write=NULL,
  quietly=T,
  normalize=T,
  n=NULL,
  TMP_PATH="/tmp/r_raster_tmp",
  MAX_CPUS=6,
  MAX_THREAD_RAM=400 # units here are number of cells to hold in an array (mem)
){
  # split-up a large raster into chunks that we can process
  # in parallel
  if(is.null(n)){
    n <- parallel::detectCores() - 1
  }
  cl <- parallel::makeCluster(n)
  # hackish way of loading SpaDES package on our cluster and apportioning tmp files
  parallel::clusterExport(
    cl,
    varlist=c("TMP_PATH"),
    envir=environment()
  )
  ret <- unlist(
      parallel::clusterApply(
          cl,
          x=rep(1,n),
          fun=function(x){
            stopifnot(require(raster));
            stopifnot(require(SpaDES));
            raster::rasterOptions(tmpdir=TMP_PATH);
            setPaths(
              paste(TMP_PATH, "/cache", sep=""),
              paste(TMP_PATH, "/input", sep=""),
              paste(TMP_PATH, "/modules", sep=""),
              paste(TMP_PATH, "/output", sep="")
            )
          })
  ); rm(ret);
  # how many bands are we working with?
  bands <- raster::nlayers(explanatory_vars)
  # grab the names of our raster variables
  names <- names(explanatory_vars)
  # export our explanatory_vars
  parallel::clusterExport(
    cl,
    varlist=c("explanatory_vars", "n", "bands"),
    envir=environment()
  )
  # produce a chunked raster surface (this takes ~ 1 hour)
  if(!quietly) cat(" -- chunking our input dataset into tiles (takes ~ 1 hr)\n")
  chunks <- parallel::parLapply(
    cl=cl,
    X=1:bands,
    fun=function(band){
      explanatory_vars <- raster::subset(explanatory_vars, band)
      return(splitRaster(
        explanatory_vars,
        nx=ceiling(sqrt(n)),
        ny=ceiling(sqrt(n)))
      )
  })
  # clean-up our cluster
  parallel::stopCluster(cl)
  rm(cl)
  # start from scratch -- be very conscious of RAM with our clustering here
  cl <- parallel::makeCluster(ifelse(n > MAX_CPUS, MAX_CPUS, n))
  # export our chunks
  parallel::clusterExport(
    cl,
    varlist=c("chunks", "m", "TMP_PATH", "MAX_THREAD_RAM"),
    envir=environment()
  )
  # clean-up our cluster in prep for our chunking operation
  ret <- parallel::clusterApply(
    cl,
    x=rep(1,n),
    fun=function(x){
      require(raster)
      raster::rasterOptions(tmpdir=TMP_PATH, maxmemory=MAX_THREAD_RAM);
    }
  ); rm(ret);

  # parallelize our raster prediction across our tiles (this takes ~17 hours)
  predicted_suitability <- parallel::parLapply(
    cl=cl,
    X=1:length(chunks[[1]]), # number of tiles per-chunk
    fun=function(tile){
      # get the focal tile across our bands (chunks)
      chunks <- raster::stack(
        lapply(chunks, FUN=function(x){ return(x[[tile]]) })
      )
      # return the normalized output for the focal tile
      return(round(raster::predict(
        object=chunks,
        model=m,
        type='response',
        na.rm=T,
        inf.rm=T,
        datatype='INT1U'
      )  * 100 ))
    }
  )
  # mosaic our tiles together
  predicted_suitability <- do.call(
    raster::mosaic,
    predicted_suitability
  )
  # clean-up our cluster
  parallel::stopCluster(cl)
  rm(cl)
  # write to disk?
  if(!is.null(write)){
    raster::writeRaster(
      x=predicted_suitability,
      filename=write,
      format='GTiff',
      datatype='INT1U',
      overwrite=T
    )
  } else {
    # return to user
    return(predicted_suitability)
  }
}
#' short-hand guassian smoother function to smooth out slivers --
#' useful for feature extraction from a noisy raster surface (e.g.,
#' from random forests)
#' @export
smoother <- function(
  r=NULL,
  smoothing_fun=mean,
  n_pixels=33,
  progress=null,
  ...
){
    return(raster::focal(
        r,
        w=raster::focalWeight(r, d=raster::res(r)[1]*n_pixels, type='Gauss'),
        progress=progress,
        fun=smoothing_fun
    ))
}
