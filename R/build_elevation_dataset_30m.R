# Author: Kyle Taylor (kyle.taylor@pljv.org)

#
# Please don't edit this script or move the resulting TIF files unless you know where they
# go -- other programs we've written depend on very specific paths for these elevation 
# products. KT (5-10-2018)
#
# see this thread : https://timwhitlock.info/blog/2010/04/google-maps-zoom-scales/
# for a discussion on going from AWS zoom-levels to meters
#

TMP_DIR = "/tmp/ned_data"

stopifnot(require(FedData))
#require(elevatr)
stopifnot(require(rgdal))
stopifnot(require(raster))
stopifnot(require(rgeos))
 
pljv_region_extent <- rgdal::readOGR(
    "/gis_data/PLJV/",
    "pljv_extent_bounding_box",
    verbose=F
  )

pljv_region_extent <- rgeos::gBuffer(pljv_region_extent, byid=F, width=2000)

#elevation_raster <- elevatr::get_elev_raster(
#    pljv_region_extent, 
#    src="aws", 
#    z=12
#  )

unlink(TMP_DIR, recursive=T, force=T)
dir.create(TMP_DIR)

elevation_raster <- get_ned(
    template=pljv_region_extent, 
    label="elev_workspace",
    extraction.dir=paste(TMP_DIR,"/extracted",sep=""),
    raw.dir=paste(TMP_DIR,"/raw",sep=""),
    force.redo=T,
    res=1 # 1/3 arc-second
  )

writeRaster(elevation_raster, "30m_elevation_data.tif", progress='text')
