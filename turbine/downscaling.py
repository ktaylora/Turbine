"""
Performs statistical downscaling (using a 'robust linear estimator') of wind 
data from NREL's ~2 km grid to a 250 meter grid that we can use to attribute 
wind data to the 1 KM US National Grid
"""
import warnings
warnings.filterwarnings("ignore")

import logging
logger = logging.getLogger(__name__)

import sys
import os

import elevation
import numpy as np

from osgeo.gdal_array import NumericTypeCodeToGDALTypeCode as dtype_to_gdal

import rasterio as rio
from rasterstats import zonal_stats

from beatbox import Vector, Raster

from beatbox.moving_windows import ndimage_filter
from beatbox.raster import slope, aspect, NdArrayDiscCache

#from .wind_toolkit import generate_h5_grid_geodataframe
#from .wind_toolkit import attribute_and_bootstrap_timeseries

def _geot_to_affine(geot):
    c, a, b, f, d, e  = list(geot)
    return( rio.Affine(a,b,c,d,e,f) )
    
def write_raster(array=None, filename=None, template=None, **kwargs):
    """
    Wrapper for rasterio that can write NumPy arrays to disc using an optional
    Raster template object
    """
    
    kwargs['driver'] = kwargs.get('driver', 'GTiff')
    kwargs['dtype'] = kwargs.get('dtype', str(array.dtype))
    kwargs['width'] = kwargs.get('width', array.shape[0])
    kwargs['height'] = kwargs.get('height', array.shape[1])
    kwargs['count'] = kwargs.get('count', 1)
    kwargs['crs'] = kwargs.get('crs', None)
    kwargs['transform'] = kwargs.get('transform', None)

    if template is not None:
        kwargs['transform'] = _geot_to_affine(template.geot)
        kwargs['crs'] = template.projection.ExportToProj4()
    if kwargs['crs'] is None:
        logger.debug('crs= was not specified and cannot be determined from a '+
            'numpy array; Resulting GeoTIFF will have no projection.')
    if kwargs['transform'] is None:
        logger.debug('transform= was not specified; Resulting GeoTIFF will '+ 
            'have an undefined affine transformation.')
  
    try:
        with rio.open(filename, 'w', **kwargs) as dst:
            dst.write(array, 1)
            
        return(True)
        
    except FileNotFoundError:
        logger.exception('FileNotFoundError in filename= argument of write_raster():'+
            'This should not happen -- are you writing to a weird dir?')
        return(False)
    
    return(False)

def _gdal_rasterize(shapefile_path, template=None, outfile=None, **kwargs):
    """
    Rasterize a vector dataset stored as a shapefile on disk 
    using gdal.RasterizeLayer. This is slow and requires a lot of disk I/O.
    """
    kwargs['dst_filename'] = kwargs.get('dst_filename', outfile)
    
    kwargs['eType'] = kwargs.get('eType', dtype_to_gdal(template.array.dtype))
    kwargs['xsize'] = kwargs.get('xsize', abs(template.geot[1]))
    kwargs['ysize'] = kwargs.get('ysize', abs(template.geot[5]))
    kwargs['bands'] = kwargs.get('bands', 1)

    target_ds = gdal.GetDriverByName('GTiff').Create(**kwargs)
    layer = ogr.Open(shapefile_path).get_layer()

    gdal.RasterizeLayer(dataset=target_ds, bands=kwargs['bands'], 
        layer=layer, burn_values=1)
    
    target_ds.FlushCache()

def rasterize(obj, template=None, **kwargs):
    """
    Wrapper for rasterio.features.rasterize that will accept a geopandas
    dataframe of features and burn the vector geometries using a Raster 
    template object
    """
    
    kwargs['field'] = kwargs.get('field', 0)
    kwargs['dtype'] = kwargs.get('dtype', template.array.dtype)

    shapes = ((geom,value) for geom, value in zip(obj.geometry, obj[field]))

    test = rio.features.rasterize(
        shapes = shapes, 
        fill = 0, 
        out = np.zeros(shape=template.array.shape, dtype = kwargs['dtype']), 
        transform = _geot_to_affine(template.geot)) 
  
def get_project_region_extent(config_file='config.json', 
    sql_options={'table_name':'mboggie.project_region_buffered_1km'}):
    """
    Query our SQL database and determine project region extent
    """
    logger.debug("Querying database for project region extent")
    project_region = Vector(
        input=config_file,
        options=sql_options)
    project_region.crs = '+init=epsg:2163'
    return project_region.to_geodataframe().\
        to_crs("+init=epsg:4326").geometry.total_bounds

def get_ned_elevation_raster(extent, res=250, **kwargs):
    """
    Wrapper for elevation.clip that will generate a GeoTIFF file for a 
    given project extent.
    :param output: Path to output file. Existing files will be overwritten.
    :param cache_dir: Root of the DEM cache folder.
    :param product: DEM product choice.
    """

    PRODUCT = {
      '30':elevation.PRODUCTS[0],
      '250':elevation.PRODUCTS[1]
    }

    logger.debug("Fetching NED raster data for project region extent")

    if kwargs.get('output') is None:
        raise AttributeError('output= argument is required by '+
            'elevation.clip() and cannot be None')
    if os.path.isfile(kwargs.get('output')):
        logger.warning("Warning : It looks like output="+
            kwargs.get('output')+" already exists... skipping.")

    else:
        elevation.clip(bounds=extent, product=PRODUCT[str(res)], **kwargs)
        elevation.clean()

if __name__ == "__main__":

    logger.debug('Checking for elevation (NED) products')

    elev = 'raster/elevation.tif'
    slp = 'raster/slope.tif'
    asp = 'raster/aspect.tif'
    tri_3 = 'raster/tri_3x3.tif'
    tri_11 = 'raster/tri_11x11.tif'

    if not os.path.isfile(elev):
        get_ned_elevation_raster(get_project_region_extent(), 
            output='elevation.tif', res=250, max_download_tiles=500)
        os.rename(os.path.expanduser("~") +
            '.cache/elevation/SRTM3/elevation.tif', "raster/elevation.tif")

    elev = Raster(elev, use_disc_caching=False)

    logger.debug('Generating slope / aspect / topographic roughness products')

    if not os.path.isfile(slp):
        slp = slope(np.ndarray.astype(elev.array,np.float), 
            use_disc_caching=True)
        write_raster(slp.array, 'raster/slope.tif', template=elev)
    else:
        slp = Raster(input=slp)

    if not os.path.isfile(asp):
        asp = aspect(np.ndarray.astype(elev.array,np.float), 
            use_disc_caching=True)
        write_raster(asp.array, 'raster/aspect.tif', template=elev)
    else:
        asp = Raster(input=asp)

    if not os.path.isfile(tri_3):
        tri_3 = ndimage_filter(
            elev.array,
            use_disc_caching=False,
            function=np.std,
            size=3
        )
        write_raster(tri_3, 'raster/tri_3x3.tif', template=elev)
    else:
        tri_3 = Raster(input=tri_3)

    if not os.path.isfile(tri_11):
        tri_11 = ndimage_filter(
            elev.array,
            use_disc_caching=False,
            function=np.std,
            size=11
        )
        write_raster(tri_11, 'raster/tri_11x11.tif', template=elev)
    else:
        tri_11 = Raster(input=tri_11)

    logger.debug('Checking for NREL Wind Toolkit HSDS grid')
    
    nrel_grid = Vector(
      'config.json', 
      sql_options = {
          'table_name':'mboggie.project_region_buffered_1km'
    }).to_geodataframe()
      
    nrel_grid = generate_h5_grid_geodataframe(
      filter_by_intersection = nrel_grid)
    
    wind_products = [
        'windspeed_200m', 
        'windspeed_100m', 
        'windspeed_40m', 
        'winddirection_200m',
        'winddirection_100m',
        'winddirection_40m',
        'GHI']
        
    nrel_grid = attribute_and_bootstrap_timeseries(
      gdf = nrel_grid,
      datasets = wind_products,
      n_bootstrap_replicates = 30)
    
    # use a bilinear interpolation to downscale our attributed NREL grid to our 250 meter NED grid
    wind_rasters = [ rasterize(nrel_grid, template=elev, field=f, dtype=np.dtype('float32')) for f in wind_products ]
    wind_rasters = [ write_raster(array=wind_rasters[i], filename='raster/'+f+'.tif', template=elev) 
        for i, f in enumerate(wind_products) ]
    # use a robust regression to estimate NREL station-level wind variable ~ f(bilinear interpolation + topographic variables)
    # write 30 meter downscaled raster surfaces to disc for review
    # aggregate mean + variance statistics of the downscaled grid to our US National grid units
    # write US National Grid units to PostGIS
