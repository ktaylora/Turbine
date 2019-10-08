"""
Manage downloading and ingesting tasks for NREL's HSDS Service
"""

import warnings
warnings.filterwarnings('ignore')

from shapely.geometry import Point
from geopandas import GeoSeries, GeoDataFrame

import gdal, ogr
import h5pyd as h5

def extent_to_hsds_bounding_box(extent=None):
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
    del target_ds

def h5_to_geodataframe(dataset_name=None):
    f = h5.File("/nrel/wtk-us.h5", 'r')
    ds = f[dataset_name]
    #x,y = zip(*f['coordinates'].flatten())
    coords = f['coordinates'][:].flatten()
    _gdf = GeoDataFrame({
        'geometry' : GeoSeries([Point(i) for i in coords])
    })
    # merge in our attributes
    _gdf = _gdf.join(ds)

