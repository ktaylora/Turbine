#!/usr/bin/Rscript

#
# Takes a while to run, but this should run like a tank.
#
# Author : Kyle Taylor (2019)
# License : GPL v.3 (2019)
#

#
# RUNTIME ARGUMENTS
#

TMP_DIR = "/tmp/topographic_variables"
PRODUCTS_DIR = "/gis_data/Elevation/Deliverables"
DEBUGGING = T

#
# LOCAL FUNCTIONS
#

debug <- function(x) if(DEBUGGING) cat("DEBUG:", paste(x, collapse=" "), "\n")

#
# MAIN
#

unlink(
    TMP_DIR, 
    recursive=T, 
    force=T
  )

dir.create(TMP_DIR)

stopifnot(require(raster))

cat(" -- setting up project workspace\n")
rasterOptions(tmpdir=TMP_DIR)
if(!file.exists("products")){
  dir.create("products")
}

cat(" -- reading existing raster data\n")
existing <- list.files(PRODUCTS_DIR, pattern="tif$")
elev <- raster::raster("30m_elevation_data.tif");

cat(" -- processing moving windows (fun=sd) for topographic roughness:\n")
windows <- c(3,11,33,107)
for(x in windows){
  cat(" -- window:",x)
  outfile <- paste(c("topographic_roughness_",x,"x",x,".tif"),sep="",collapse="")
  if(!outfile %in% existing){
    focal <- raster::focal(elev, w=matrix(1,nrow=x, ncol=x), fun=sd, na.rm=T)
    raster::writeRaster(focal,paste(PRODUCTS_DIR,"/topographic_roughness_",x,"x",x,".tif",sep="",collapse=""), overwrite=T)
    cat(" [done]\n")
  } else {
    cat(" [previous]\n")
  }
}

suppressWarnings(rm(focal,outfile,windows));
gc();

cat(" -- calculating slope, aspect, flow direction, and topographic position index:\n")
result <- terrain(elev, opt=c("slope","aspect","flowdir","TPI"), progress='text')

writeRaster(
    result[[1]],
    paste(PRODUCTS_DIR,"/slope.tif",sep=""), 
    progress='text', 
    overwrite=T
  )

writeRaster(
    result[[2]],
    paste(PRODUCTS_DIR,"/aspect.tif",sep=""),
    progress='text',
    overwrite=T
  )

writeRaster(
    result[[3]],
    paste(PRODUCTS_DIR,"/flow_direction.tif",sep=""),
    progress='text',
    overwrite=T
  )

writeRaster(
    result[[4]],
    paste(PRODUCTS_DIR,"/topographic_position_index.tif",sep=""),
    progress='text', 
    overwrite=T
  )

unlink(
    TMP_DIR,
    recursive=T,
    force=T
  )
  
quit("no", status=0)
