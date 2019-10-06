"""
Performs statistical downscaling (using a 'robust linear estimator') of wind data from NREL's ~2 km grid to a 250 meter grid that we
can use to attribute wind data to the 1 KM US National Grid
"""

import logging
logger = logging.getLogger(__name__)

import warnings
warnings.filterwarnings("ignore")

import sys
import os

import elevation
import numpy as np

from rasterstats import zonal_stats

from beatbox import Vector, Raster

from beatbox.moving_windows import ndimage_filter
from beatbox.raster import slope, aspect

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

def get_ned_elevation_raster(extent, res=250, **kwargs):
    """
    Wrapper for elevation.clip that will generate a GeoTIFF file for a given project extent.
    :param output: Path to output file. Existing files will be overwritten.
    :param cache_dir: Root of the DEM cache folder.
    :param product: DEM product choice.
    """
    PRODUCT = {'30':elevation.PRODUCTS[0],'250':elevation.PRODUCTS[1]}
    logger.debug("Fetching NED raster data for project region extent")
    if kwargs.get('outfile') is None:
        raise AttributeError('outfile= argument is required by elevation.clip() and cannot be None')
    if os.path.isfile(kwargs.get('outfile')):
        logger.warning("Warning : It looks like outfile="+kwargs.get('outfile')+" already exists... skipping.")
        return kwargs.get('outfile')
    else:
        elevation.clip(bounds=extent, product=PRODUCT[str(res)], **kwargs)
        elevation.clean()
        return kwargs.get('outfile')


if __name__ == '__main__':

  logger.debug('Checking for elevation (NED) products')

  elev = 'raster/elevation.tif'
  if not os.path.isfile(elev):
      elev = get_ned_elevation_raster(get_project_region_extent(), output='elevation.tif', res=250, max_download_tiles=500)
  elev = Raster(input=elev, use_disc_caching=True)

  logger.debug('Generating slope / aspect / topographic roughness products')

  slp = Raster(use_disc_caching=True)
  asp = Raster(use_disc_caching=True)
  slp.array = slope(elev.array)
  asp.array = aspect(elev.array)

  topographic_roughness_3x3 = ndimage_filter(image=elevation)

# use a bilinear interpolation to downscale our attributed NREL grid to our 30 meter NED grid
# use a robust regression to estimate NREL station-level wind variable ~ f(bilinear interpolation + topographic variables)
# write 30 meter downscaled raster surfaces to disc for review
# aggregate mean + variance statistics of the downscaled grid to our US National grid units
# write US National Grid units to PostGIS
