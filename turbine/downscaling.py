"""
Performs statistical downscaling (using a 'robust linear estimator') of wind data from NREL's ~2 km grid to a 250 meter grid that we
can use to attribute wind data to the 1 KM US National Grid
"""
import warnings
warnings.filterwarnings("ignore")

import logging
logger = logging.getLogger(__name__)

import sys
import os

import elevation
import numpy as np

import rasterio as rio

from rasterstats import zonal_stats

from beatbox import Vector, Raster

from beatbox.moving_windows import ndimage_filter
from beatbox.raster import slope, aspect, NdArrayDiscCache, _NUMPY_TYPES

def write_raster(array=None, filename=None, **kwargs):
    
    kwargs['driver'] = kwargs.get('driver', 'GTiff')
    kwargs['dtype'] = kwargs.get('dtype', str(array.dtype))
    kwargs['width'] = kwargs.get('width', array.shape[0])
    kwargs['height'] = kwargs.get('height', array.shape[1])
    kwargs['count'] = kwargs.get('count', 1)

    try:
        with rio.open(os.path.abspath(filename), 'w', **kwargs) as dst:
            dst.write(array, 1)
    except FileNotFoundError:
        raise FileNotFoundError('in filename= argument of write_raster()')

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
    slp = 'raster/slope.tif'
    asp = 'raster/aspect.tif'
    tri_3 = 'raster/tri_3x3.tif'
    tri_11 = 'raster/tri_11x11.tif'

    if not os.path.isfile(elev):
        elev = get_ned_elevation_raster(get_project_region_extent(), output='elevation.tif', res=250, max_download_tiles=500)
        os.rename('elevation.tif', "raster/elevation.tif")

    elev = Raster(elev, use_disc_caching=True)

    logger.debug('Generating slope / aspect / topographic roughness products')

    if not os.path.isfile(slp):
        slp = slope(elev.array, use_disc_caching=True)
        write_raster(slp.array, 'raster/slope.tif')
    else:
        slp = Raster(input=slp)

    if not os.path.isfile(asp):
        asp = aspect(elev.array, use_disc_caching=True)
        write_raster(asp.array, 'raster/aspect.tif')
    else:
        asp = Raster(input=asp)

    if not os.path.isfile(tri_3):
        tri_3 = ndimage_filter(
            elev.array, 
            use_disc_caching=True, 
            function=np.std, 
            size=3
        )
        write_raster(tri_3.array, 'raster/tri_3x3.tif')
    else:
        tri_3 = Raster(input=tri_3)

    if not os.path.isfile(tri_11):
        tri_11 = ndimage_filter(
            elev.array, 
            use_disc_caching=True, 
            function=np.std, 
            size=3
        )
        write_raster(tri_11.array, 'raster/tri_11x11.tif')
    else:
        tri_11 = Raster(input=tri_11)

    # Fetch our NREL data from HSDS
    # use a bilinear interpolation to downscale our attributed NREL grid to our 30 meter NED grid
    # use a robust regression to estimate NREL station-level wind variable ~ f(bilinear interpolation + topographic variables)
    # write 30 meter downscaled raster surfaces to disc for review
    # aggregate mean + variance statistics of the downscaled grid to our US National grid units
    # write US National Grid units to PostGIS
