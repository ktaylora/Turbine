"""
Manage downloading and ingesting tasks for NREL's HSDS Service
"""

import warnings
warnings.filterwarnings('ignore')

import logging
logger = logging.getLogger(__name__)

from shapely.geometry import Point
from geopandas import GeoSeries, GeoDataFrame, overlay
from pandas import DataFrame

import gdal, ogr
import h5pyd as h5

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

def extent_to_bounding_box(extent=None):
    ne = (1,1)
    sw = (2,2)
    return (ne, sw)

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


def h5_to_geodataframe(dataset_name=None, filter_by_intersection=None):
    """
    Will parse NREL's Wind Toolkit 
    """
    f = h5.File("/nrel/wtk-us.h5", 'r')
    
    logger.debug("Fetching coordinates from wind toolkit...")
    coords = f['coordinates'][:].flatten()

    target_rows = list(range(len(coords)))

    gdf = GeoDataFrame({
        'geometry' : GeoSeries([Point(i) for i in coords]),
        'id' : target_rows
    })

    if filter_by_intersection is not None:
        target_rows = list(overlay(gdf, filter_by_intersection, how='intersection')['id'])
        if len(target_rows) is 0:
            raise AttributeError('filter_by_intersection= resulted in no intersecting geometries')
    
    gdf = gdf.iloc[target_rows,:]

    logger.debug("Merging in our selected attributes")

    for dataset in dataset_name:
        gdf = gdf.join(DataFrame({
            dataset : f[dataset][:].flatten()[target_rows]
        }, index=target_rows))

    del f, coords
    return(gdf)

