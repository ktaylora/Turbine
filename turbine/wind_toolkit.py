"""
Manage downloading and ingesting tasks for NREL's HSDS Service
"""

import warnings
warnings.filterwarnings('ignore')

import logging
logger = logging.getLogger(__name__)

import os

from shapely.geometry import Point
from shapely.ops import cascaded_union

from geopandas import GeoSeries, GeoDataFrame
from pandas import DataFrame

from scipy.stats import norm

from numpy import poly1d as polynomial_regression
from numpy import unique, polyfit, mean, linspace, memmap, array, zeros, ix_, float64

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
    """
    Wrapper for scipy.stats.norm that will generate n normally distributed
    values about a user-specified mean. This is used for bootstrapping
    hour-periods used in the wind toolkit api. 
    """
    samples = list(norm.rvs(loc=mean, scale=variance, size=n_samples))
    if fun is not None:
        return [fun[0](x) for x in samples]
    return samples

def generate_h5_grid_geodataframe(filter_by_intersection=None, 
    cache_file=_HSDS_CACHE_FILE_PATH):
    """
    Query and subset the latest wind toolkit spatial grid using a user-specified
    geometry. By default, will look for a cached (local) copy of the grid to save
    bandwidth
    """
    if os.path.exists(cache_file):
        logger.debug("Using cached file for project region coordinates for" +
            " the wind toolkit")

        gdf = GeoDataFrame().from_file(cache_file)
        target_rows = gdf['id']

    else:

        logger.debug("Fetching coordinates from wind toolkit HSDS interface")

        f = h5.File("/nrel/wtk-us.h5", 'r')

        n_rows, n_cols = f['coordinates'].shape
        coords = f['coordinates'][:].flatten()

        target_rows = i = list(range(len(coords)))

        target_rows_id = [ math.ceil(x/n_cols) for x in i ]
        target_cols_id = [ math.ceil( n_cols * 
            ( float64(x/n_cols) - math.floor(float64(x/n_cols))) )  for x in i ]

        gdf = GeoDataFrame({
            'geometry' : GeoSeries([Point(reversed(i)) for i in coords]),
            'id' : i,
            'x' :  target_cols_id,
            'y' : target_rows_id
        })

        gdf.crs = _WIND_TOOLKIT_DEFAULT_EPSG
    
    gdf = gdf.iloc[target_rows,:]

    # try and cleanly flush our toolkit session
    f.close()
    del f

    if filter_by_intersection is not None:
        try:
            filter_by_intersection = filter_by_intersection.to_crs(
                _WIND_TOOLKIT_DEFAULT_EPSG)
        except ValueError:
            logger.debug("No CRS found in filter_by_intersection"+
                ". Sometimes PostGIS results returned" +
                " by Vector() won't have a CRS -- assuming EPSG:2163")
            filter_by_intersection.crs = "+init=epsg:2163"
            filter_by_intersection = filter_by_intersection.to_crs(
                _WIND_TOOLKIT_DEFAULT_EPSG)
        target_rows = gdf.loc[gdf.within(cascaded_union(
            filter_by_intersection.geometry))]['id']
        if len(target_rows) is 0:
            raise AttributeError('filter_by_intersection= resulted in no'+
                ' intersecting geometries')
        gdf = gdf.iloc[target_rows,:]

    gdf.to_file(filename=_HSDS_CACHE_FILE_PATH)

    return(gdf)

def _get_wtk_timeslice_by_hour(gdf=None, hour=None, dataset=None):
    f = h5.File("/nrel/wkt-us.h5", 'r')
    
    result = f[dataset][
        hour,
        min(unique(gdf['y'])):max(unique(gdf['y'])),
        min(unique(gdf['x'])):max(unique(gdf['x'])) ]
        
    f.close()
    del f

    return( result )

def _disc_cached_attribute_timeseries(gdf=None, timeseries=_HOURS_PER_MONTH, 
    datasets=None):

    if not isinstance(datasets, Iterable):
        datasets = [datasets]

    for dataset in datasets:
        f = h5.File("/nrel/wtk-us.h5", 'r')

        _WTK_MAX_HOURS=f[dataset].shape[0]

        if len(timeseries) == 1:
            all_hours = linspace(0, _WTK_MAX_HOURS, num=timeseries, dtype='int')
        else:
            all_hours = timeseries

        logger.debug('Building a giant cached array specifying our' +
            ' target WTK hours and sites')

        if os.path.exists('.wtk_result.dat') : os.remove('.wtk_result.dat')

        wtk_result = memmap(
            filename='.wtk_result.dat',
            dtype='float16',
            mode='w+',
            shape=(len(all_hours), len(unique(gdf['y'])), len(unique(gdf['x']))) )

        if os.path.exists('.wtk_selection.dat') : os.remove('.wtk_selection.dat')

        wtk_selection = memmap(
            filename='.wtk_selection.dat',
            dtype='bool',
            mode='w+',
            shape = f[dataset].shape)

        wtk_selection[ix_(
            array([ h in all_hours for h in range(f[dataset].shape[0]) ]),
            array([ i in gdf['y'] for i in range(f[dataset].shape[1]) ]),
            array([ i in gdf['x'] for i in range(f[dataset].shape[2]) ]))] \
                = True

        wtk_result[:] = f[dataset][wtk_selection]

        # join in our full y_overall table with all hourlies for this dataset
        # with our source WTK grid
        wtk_result.set_axis(
            [ dataset + '_' + str(h) for h in all_hours],
            axis=1,
            inplace=True)

        gdf = gdf.join(wtk_result)

        del wtk_result, wtk_selection

    return(gdf)

def _disc_cached_attribute_and_bootstrap_timeseries(gdf=None, timeseries=_HOURS_PER_MONTH,
    datasets=None, n_bootstrap_replicates=30):
    """
    Using an attributed GeoDataFrame containing our target wind toolkit grid
    ID's, attempt to fetch and attribute wind toolkit time series data for a
    focal toolkit dataset(s). This version uses large array discaching and
    array broadcasting to select wind data -- it fails for large project
    geographies
    """
    # build an argument spec for scipy.stats.norm
    _kwargs = dict()

    _kwargs['variance'] = 64
    _kwargs['fun'] = round,
    _kwargs['n_samples']  = n_bootstrap_replicates

    if isinstance(timeseries, Iterable) :
        all_hours = timeseries
    else:
        all_hours = linspace(0, _WTK_MAX_HOURS, num=timeseries, dtype='int')

    y_overall = DataFrame(zeros(shape=(len(gdf['id']),len(all_hours))))

    # pull a number of random hours around our target
    # hour and fit a polynomial regression to the time-series
    for i, hour in enumerate(all_hours):
        _kwargs['mean'] = hour

        logger.debug('Processing hour: ' + str(hour) +
            '; % complete : ' + str(round(i/len(all_hours),2)*100))

        bs_hourlies = [ int(h) for h in
            unique(_bootstrap_normal_dist(**_kwargs)) ]
        bs_hourlies = array(bs_hourlies)[array(bs_hourlies) >= 0]

        y_bs_timeseries = _disc_cached_attribute_timeseries(gdf=gdf, timeseries=bs_hourlies, datasets=dataset)

        logger.debug('Fitting polynomial regression by-site for' +
            ' full time-series')

        for i, site in y_bs_timeseries.iterrows():
            m_poly = polynomial_regression(polyfit(
                x=list(array(bs_hourlies)),
                y=list(array(site)), deg=2))

            intercept_m = round(mean(array(site)),2)

            fitted = [ round(m_poly(h), 2) for h in
                list(array(bs_hourlies)) ]

            residuals = array(site) - fitted
            residuals_intercept = array(site) - intercept_m

            null_vs_alt_sse = (sum(abs(residuals_intercept)) -
                sum(abs(residuals)))
            r_squared = round( null_vs_alt_sse /
                sum(abs(residuals_intercept)), 2 )
            if r_squared < 0.1:
                logger.debug('Poor regression estimator fit on model' +
                    ' for hour='+str(int(hour)))
                y_overall.iloc[i,array(all_hours) == hour] = \
                    round(mean(fitted),2)
            elif r_squared < 0:
                logger.debug('Null model outperformed our regression' +
                    ' estimator hour=' + str(int(hour)) +
                    '; Returning null (mean) time-series estimate:' +
                        str(intercept_m))
                y_overall.iloc[i,array(all_hours) == hour] = \
                    intercept_m
    # join in our full y_overall table with all hourlies for this dataset
    # with our source WTK grid
    y_overall.set_axis(
        [ dataset + '_' + str(h) for h in all_hours],
        axis=1,
        inplace=True)

    return( gdf.join(y_overall) )


def _original_attribute_and_bootstrap_timeseries(gdf=None, timeseries=_HOURS_PER_MONTH,
    datasets=None, n_bootstrap_replicates=30):
    """
    Using an attributed GeoDataFrame containing our target wind toolkit grid
    ID's, attempt to fetch and attribute wind toolkit time series data for a
    focal toolkit dataset(s).
    """
    # build an argument spec for scipy.stats.norm
    _kwargs = dict()

    _kwargs['variance'] = 64
    _kwargs['fun'] = round,
    _kwargs['n_samples']  = n_bootstrap_replicates

    if not isinstance(datasets, Iterable):
        datasets = [datasets]

    for dataset in datasets:

        if os.path.exists('.wtk_selection.dat') : os.remove('.wtk_selection.dat')

        # this is the old (working) implementation
        y_overall = DataFrame(zeros(shape=(len(gdf['id']),len(all_hours))))
        for i, hour in enumerate(all_hours):
            _kwargs['mean'] = hour

            logger.debug('Processing hour: ' + str(hour) +
                '; % complete : ' + str(round(i/len(all_hours),2)*100))

            f = h5.File("/nrel/wtk-us.h5", 'r')

            # otherwise, pull a number of random hours around our target
            # hour and fit a polynomial regression to the time-series

            bs_hourlies = [ int(h) for h in
                unique(_bootstrap_normal_dist(**_kwargs)) ]

            bs_hourlies = array(bs_hourlies)[array(bs_hourlies) >= 0]

            y = DataFrame(
                zeros(shape=(len(gdf['id']),
                len(bs_hourlies))))

            valid_hourlies = []

            for i, h in enumerate(list(bs_hourlies)):
                try:
                    print('Downloading hour slice (' +
                        str(hour) + ') replicate: ' + str(h) + '...')
                    y[i] = f[dataset][h,:].flatten()[gdf['id']]
                    valid_hourlies.append(h)

                except OSError:
                    logger.debug('Caught en error fetching hour slice (' +
                        str(hour) + ') Replicate: ' +
                        str(h) + ' [skipping...]')
                    y[i] = None
                    valid_hourlies.append(None)

            keep = [h is not None for h in valid_hourlies]

            logger.debug('Fitting polynomial regression by-site for' +
                ' full time-series')

            for i, site in y.iterrows():
                m_poly = polynomial_regression(polyfit(
                    x=list(array(valid_hourlies)[keep]),
                    y=list(array(site)[keep]), deg=2))

                intercept_m = round(mean(array(site)[keep]),2)

                fitted = [ round(m_poly(h), 2) for h in
                    list(array(valid_hourlies)[keep]) ]

                residuals = array(site)[keep] - fitted
                residuals_intercept = array(site)[keep] - intercept_m

                null_vs_alt_sse = (sum(abs(residuals_intercept)) -
                    sum(abs(residuals)))
                r_squared = round( null_vs_alt_sse /
                    sum(abs(residuals_intercept)), 2 )

                if r_squared < 0.1:
                    logger.debug('Poor regression estimator fit on model' +
                        ' for hour='+str(int(hour)))
                    y_overall.iloc[i,array(all_hours) == hour] = \
                        round(mean(fitted),2)
                elif r_squared < 0:
                    logger.debug('Null model outperformed our regression' +
                        ' estimator hour=' + str(int(hour)) +
                        '; Returning null (mean) time-series estimate:' +
                        str(intercept_m))
                    y_overall.iloc[i,array(all_hours) == hour] = \
                        intercept_m
            # try and cleanly flush our toolkit session
            f.close()
            del f
        # when finished with all of our bootstrapped hourlies,
        # join in our full y_overall table (all_hours) for this dataset
        # with our source WTK grid and move-on to the next dataset
        y_overall.set_axis(
            [ dataset + '_' + str(h) for h in all_hours],
            axis=1,
            inplace=True)

        gdf = gdf.join(y_overall)

    return(gdf)
