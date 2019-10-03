"""
Performs statistical downscaling (using a 'robust linear estimator') of wind data from NREL's ~2 km grid to a 30 meter grid that we
can use to attribute wind data to the 1 KM US National Grid
"""

import logging

logger = logging.getLogger(__name__)

import elevation
from rasterstats import zonal_stats

from beatbox import Vector

def get_project_region_extent(config_file='config.json', sql_options={'table_name':'mboggie.project_region_buffered_1km'}):
    """
    Query our SQL database and determine project region extent
    """
    logger.debug("Querying database for project region extent")
    project_region = Vector(
        input=config_file, 
        options=sql_options)
    project_region.crs = '+init=epsg:2163'
    return project_region.to_geodataframe().to_crs("+init=epsg:4326").geometry.total_bounds

def get_ned_elevation_rasters(extent, res=30, **kwargs):
    """
    Wrapper for elevation.clip that will generate a GeoTIFF file for a given extent.
    :param output: Path to output file. Existing files will be overwritten.
    :param cache_dir: Root of the DEM cache folder.
    :param product: DEM product choice.
    """
    logger.debug("Fetching NED raster data for project region extent")
    elevation.clip(bounds=extent, **kwargs)
    elevation.clean()


if __name__ == '__main__':
    
    get_ned_elevation_rasters(get_project_region_extent(), output='elevation.tif', max_download_tiles=500)

# generate the usual suspects for topographic variables

# use a bilinear interpolation to downscale our attributed NREL grid to our 30 meter NED grid
# use a robust regression to estimate NREL station-level wind variable ~ f(bilinear interpolation + topographic variables)
# write 30 meter downscaled raster surfaces to disc for review
# aggregate mean + variance statistics of the downscaled grid to our US National grid units
# write US National Grid units to PostGIS
