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

from scipy.stats import norm
from numpy import poly1d as polynomial_regression
from numpy import unique, polyfit, mean

import gdal, ogr
import h5pyd as h5

import math

_WIND_TOOLKIT_DEFAULT_EPSG = "+init=epsg:4326"
_HOURS_PER_MONTH = 730
_HSDS_CACHE_FILE_PATH = 'vector/h5_grid.shp'

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

def _bootstrap_normal_dist(n_samples=10, mean=0, variance=2, fun=None):
    samples = list(norm.rvs(loc=mean, scale=variance, size=n_samples))
    if fun is not None:
        return [fun(x) for x in samples]
    return samples

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


def generate_h5_grid_geodataframe(datasets=None, filter_by_intersection=None, cache_file=_HSDS_CACHE_FILE_PATH, bootstrap_timeseries=0):
    """
    Will parse NREL's Wind Toolkit as efficiently as possible and 
    """
    f = h5.File("/nrel/wtk-us.h5", 'r')
    
    if os.path.exists(cache_file):
        logger.debug("Using cached file for project region coordinates for the wind toolkit")
        
        gdf = GeoDataFrame().from_file(cache_file)
        target_rows = gdf['id']
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
    	gdf = gdf.iloc[target_rows,:] 
    	
    	gdf.to_file(filename=_HSDS_CACHE_FILE_PATH)
        
    if filter_by_intersection is not None:
        target_rows = list(gdf.loc[gdf.within(cascaded_union(filter_by_intersection.geometry))]['id'])
        if len(target_rows) is 0:
            raise AttributeError('filter_by_intersection= resulted in no intersecting geometries')
        gdf = gdf.iloc[target_rows,:] 
    
    f.close()
    del f # try and cleanly flush our toolkit session 
    
    return(gdf)

def fetch_and_attribute_timeseries(gdf=None, timeseries=_HOURS_PER_MONTH, datasets=None, boostrap_timeseries=0):
    """
    Using an attributed GeoDataFrame containing our target wind toolkit grid ID's,
    attempt to fetch and attribute wind toolkit time series data for a focal toolkit dataset(s).
    """
      
    f = h5.File("/nrel/wtk-us.h5", 'r')
    
    for dataset in datasets:
        if bootstrap_timeseries > 0:
         
         # build an argument spec for scipy.stats.norm
          _kwargs = dict()
          
          _kwargs['variance'] = 12 
          _kwargs['fun'] = round, 
          _kwargs['n_samples']  = bootstrap_timeseries 
          
          all_hours = list(np.linspace(0, 61368, num=timeseries, dtype='int'))
          
          # let y be a list of DataFrames where each element is a bootstrap
          # replicate of our wind toolkit 'hourlies'. We will use the fitted 
          # values from a quadratic regression across replicates for average each row 
          # (x,y coordinate) of the tables -- this may blow the RAM out of 
          # the machine... and NREL may send angry emails about usage.
          
          y_overall = DataFrame(np.empty(shape=(len(gdf['id']),len(all_hours))))
          
          for hour in all_hours:
              _kwargs['mean'] = hour
              
              bs_hourlies = [ int(h) for h in unique(_bootstrap_normal_dist(**_kwargs)) ]
              
              y = DataFrame(np.empty(shape=(len(gdf['id']),len(bs_hourlies))))
              valid_hourlies = []
          
              for i, h in enumerate(bs_hourlies):
                  try:
                    y[i] = f[dataset][h,:].flatten()[gdf['id']]
                    valid_hourlies.append(h)
                  except OSError:
                    y[i] = None
                    valid_hourlies.append(None)
                    
              keep = [h is not None for h in valid_hourlies]
                    
              for index, row in y.iterrows():
                  m_poly = polynomial_regression(polyfit(
                      x=list(np.array(valid_hourlies)[keep]), 
                      y=list(np.array(row)[keep]), deg=2))
                  
                  intercept_m = round(np.mean(np.array(row)[keep]),2)
                  
                  fitted = [ round(m_poly(h),2) for h in 
                    list(np.array(valid_hourlies)[keep]) ]        
                  
                  residuals = np.array(row)[keep] - fitted
                  residuals_intercept = np.array(row)[keep] - intercept_m
                  
                  null_vs_alt_sse = (sum(abs(residuals_intercept)) - sum(abs(residuals)))
                  r_squared = round( null_vs_alt_sse / sum(abs(residuals_intercept)), 2 ) 
                  
                  if r_squared < 0.1:
                      logger.debug("poor regression estimator fit on model for hour="+
                        str(int(hour)))
 
                  y_overall.iloc[index,np.array(all_hours) == hour] = \
                    round(mean(fitted),2)
          
        else:
          gdf = gdf.join(DataFrame({
            dataset : f[dataset][::_HOURS_PER_MONTH, gdf['y'], gdf['x']]
          }, index=target_rows))

    f.close()
    del f # try and cleanly flush our toolkit session 
    
    return(gdf)

