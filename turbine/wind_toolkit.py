"""
Manage downloading and ingesting tasks for NREL's HSDS Service
"""

import warnings
warnings.filterwarnings('ignore')

import logging
logger = logging.getLogger(__name__)

from shapely.geometry import Point
from shapely.ops import cascaded_union

from geopandas import GeoSeries, GeoDataFrame
from pandas import DataFrame

import gdal, ogr
import h5pyd as h5

import math

_WIND_TOOLKIT_DEFAULT_EPSG = "+init=epsg:4326"
_HOURS_PER_MONTH = 730
_HSDS_CACHE_FILE_PATH = '/vector/h5_grid.shp'

_wind_toolkit_datasets = [
    "DIF",
    "DNI",
    "GHI",
    "inversemoninobukhovlength_2m",
    "precipitationrate_0m",
    "pressure_0m",
    "pressure_100m",
    "pressure_200m",
    "relativehumidity_2m",
    "temperature_100m",
    "temperature_10m",
    "temperature_120m",
    "temperature_140m",
    "temperature_160m",
    "temperature_200m",
    "temperature_2m",
    "temperature_40m",
    "temperature_60",
    "temperature_80m",
    "winddirection_100m",
    "winddirection_10m",
    "winddirection_120m",
    "winddirection_140m",
    "winddirection_160m",
    "winddirection_200m",
    "winddirection_40m",
    "winddirection_60m",
    "winddirection_80m",
    "windspeed_100m",
    "windspeed_10m",
    "windspeed_120m",
    "windspeed_140m",
    "windspeed_160m",
    "windspeed_200m",
    "windspeed_40m",
    "windspeed_60m",
    "windspeed_80m"
]

def _bootstrap_gaussian_sample(n_samples=10, mean=0, sd=2):
    pass

def rasterize(shapefile_path=None, template=None, outfile=None, **kwargs):
    """
    Quick and hackish approach to rasterizing a vector dataset stored as a shapefile on disk 
    using gdal.RasterizeLayer. I haven't found a faster way to do this.
    """
    kwargs['dst_filename'] = kwargs.get('dst_filename', outfile)
    kwargs['eType'] = kwargs.get('eType', gdal.GDT_Int32)
    kwargs['xsize'] = kwargs.get('xsize', abs(template.geot[1]))
    kwargs['ysize'] = kwargs.get('ysize', abs(template.geot[5]))
    kwargs['bands'] = kwargs.get('bands', 1)

    target_ds = gdal.GetDriverByName('GTiff').Create(**kwargs)
    layer = ogr.Open(shapefile_path).get_layer()

    gdal.RasterizeLayer(dataset=target_ds, bands=kwargs['bands'], layer=layer, burn_values=1)
    
    target_ds.FlushCache()


def h5_grid_to_geodataframe(datasets=None, filter_by_intersection=None, cache_file=_HSDS_CACHE_FILE_PATH):
    """
    Will parse NREL's Wind Toolkit as efficiently as possible and 
    """
    f = h5.File("/nrel/wtk-us.h5", 'r')
    
    if os.path.exists(cache_file):
        logger.debug("Using cached file for project region coordinates for the wind toolkit")
        gdf = GeoDataFrame().from_file(cache_file)

    else:

        logger.debug("Fetching coordinates from wind toolkit HSDS interface")
        n_rows, n_cols = f['coordinates'].shape
        coords = f['coordinates'][:].flatten()

        target_rows = i = list(range(len(coords)))
    
        target_rows_id = [ math.ceil(x/n_cols) for x in i ]
        target_cols_id = [ math.ceil( n_cols * ( float(x/n_cols) - math.floor(x/n_cols) ) ) for x in i ]

        gdf = GeoDataFrame({
            'geometry' : GeoSeries([Point(reversed(i)) for i in coords]),
            'id' : i,
            'x' :  target_cols_id,
            'y' : target_rows_id
        })

    	gdf.crs = _WIND_TOOLKIT_DEFAULT_EPSG
    	gdf.to_file(_HSDS_CACHE_FILE_PATH)

    if filter_by_intersection is not None:
        target_rows = list(gdf.loc[gdf.within(cascaded_union(filter_by_intersection.geometry))]['id'])
        if len(target_rows) is 0:
            raise AttributeError('filter_by_intersection= resulted in no intersecting geometries')
    
    gdf = gdf.iloc[target_rows,:]

    logger.debug("Merging in our selected attributes")

    for dataset in datasets:
        gdf = gdf.join(DataFrame({
            dataset : f[dataset][::_HOURS_PER_MONTH, gdf['y'], gdf['x']].flatten()
        }, index=target_rows))

    del i, target_rows, f, coords
    return(gdf)

