import warnings

warnings.filterwarnings("ignore")

import logging

logger = logging.getLogger(__name__)

import sys
import os

try:
    from geopandas import GeoSeries, GeoDataFrame
except ModuleNotFoundError:
    #!{sys.executable} -m pip install geopandas --upgrade
    from geopandas import GeoSeries, GeoDataFrame

from shapely.geometry import Point
from shapely.ops import cascaded_union

from pandas import DataFrame

from scipy.stats import norm

from numpy import poly1d as polynomial_regression
from numpy import (
    unique,
    polyfit,
    mean,
    linspace,
    memmap,
    array,
    zeros,
    ix_,
    float64,
    prod,
    array_split,
)

import h5pyd as h5

import math

try:
    from tqdm import tqdm
except ModuleNotFoundError:
    #!{sys.executable} -m pip install tqdm --upgrade
    from tqdm import tqdm

_WIND_TOOLKIT_DEFAULT_EPSG = "+init=epsg:4326"
_HOURS_PER_MONTH = 730
_HSDS_CACHE_FILE_PATH = "vector/h5_grid.shp"

N_BOOTSTRAP_REPLICATES = 30
DOY_VARIANCE_PARAMETER = 60  # std. dev. parameter for days used for bootstrapping


def _bootstrap_normal_dist(n_samples=10, mean=0, variance=2, fun=None):
    """
    Wrapper for scipy.stats.norm that will generate n normally distributed
    values about a user-specified mean. This is used for bootstrapping
    hour-periods used in the wind toolkit api.
    :param n_samples:
    :param mean:
    :param variance:
    :param fun:
    :return:
    """
    samples = list(norm.rvs(loc=mean, scale=variance, size=n_samples))
    if fun is not None:
        return [int(fun(x)) for x in samples]
    return samples


def generate_h5_grid_geodataframe(
    filter_by_intersection=None, cache_file=_HSDS_CACHE_FILE_PATH
):
    """
    Query and subset the latest wind toolkit spatial grid using a user-specified
    geometry. By default, will look for a cached (local) copy of the grid to save
    bandwidth
    :param filter_by_intersection:
    :param cache_file:
    :return:
    """
    if os.path.exists(cache_file):
        logger.debug(
            "Using cached file for project region coordinates for the wind toolkit"
        )

        gdf = GeoDataFrame().from_file(cache_file)
        target_rows = gdf["id"]

    else:

        logger.debug("Fetching coordinates from wind toolkit HSDS interface")

        # f = h5.File("/nrel/wtk-us.h5", "r")
        f = h5.File("/nrel/wtk-us.h5", "r", bucket="nrel-pds-hsds")

        n_rows, n_cols = f["coordinates"].shape
        coords = f["coordinates"][:].flatten()

        target_rows = i = list(range(len(coords)))

        target_rows_id = [int(math.ceil(x / n_cols)) for x in i]
        target_cols_id = [
            int(round(n_cols * (float64(x / n_cols) - math.floor(float64(x / n_cols)))))
            for x in i
        ]

        gdf = GeoDataFrame(
            {
                "geometry": GeoSeries([Point(reversed(i)) for i in coords]),
                "id": i,
                "x": target_cols_id,
                "y": target_rows_id,
            }
        )

        gdf.crs = _WIND_TOOLKIT_DEFAULT_EPSG

    gdf = gdf.iloc[target_rows, :]

    # try and cleanly flush our toolkit session
    f.close()
    del f

    if filter_by_intersection is not None:
        try:
            filter_by_intersection = filter_by_intersection.to_crs(
                _WIND_TOOLKIT_DEFAULT_EPSG
            )
        except ValueError:
            logger.debug(
                "No CRS found in filter_by_intersection"
                + ". Sometimes PostGIS results returned"
                + " by Vector() won't have a CRS -- assuming EPSG:2163"
            )
            filter_by_intersection.crs = "+init=epsg:2163"
            filter_by_intersection = filter_by_intersection.to_crs(
                _WIND_TOOLKIT_DEFAULT_EPSG
            )
        target_rows = gdf.loc[
            gdf.within(cascaded_union(filter_by_intersection.geometry))
        ]["id"]
        if len(target_rows) is 0:
            raise AttributeError(
                "filter_by_intersection= resulted in no" + " intersecting geometries"
            )
        gdf = gdf.iloc[target_rows, :]

    gdf.to_file(filename=_HSDS_CACHE_FILE_PATH)

    return gdf


def _polynomial_ts_estimator(y=None, x=None, degree=2):
    m_poly = polynomial_regression(
        polyfit(x=list(array(x)), y=list(array(y)), deg=degree)
    )

    intercept_m = round(mean(array(y)), 2)

    fitted = [round(m_poly(hourly), 2) for hourly in list(array(x))]

    residuals = array(y) - fitted
    residuals_intercept = array(y) - intercept_m

    null_vs_alt_sse = sum(abs(residuals_intercept)) - sum(abs(residuals))
    r_squared = round(null_vs_alt_sse / sum(abs(residuals_intercept)), 2)

    if r_squared < 0.1:
        logger.debug(
            "Warning: Poor regression estimator fit on model hourly ~"
            + str(intercept_m)
            + "; Keeping value, but review output before using."
        )
    elif r_squared < 0:
        logger.debug(
            "Null model outperformed our regression"
            + " estimator hourly ~"
            + str(intercept_m)
            + "; Returning null (mean) time-series estimate:"
            + str(intercept_m)
        )
        return intercept_m

    return mean(round(mean(fitted), 2))


def _query_timeseries(
    f=None, y=None, x=None, hour=None, dataset=None, max_hours=WTK_MAX_HOURS
):
    """
    Accepts a GeoDataFrame of hdf5 grid points attributed with
    (y, x) coordinates and query the hsds interface for an hourly
    time-series values for a dataset specified by the user. This
    implementation will treat the hourly time-series as literal
    and will not do a regression to estimate dataset values
    :param gdf:
    :param timeseries:
    :param datasets:
    :return:
    """
    if not isinstance(timeseries, list):
        logger.debug("Building a sequence from user-specified single scalar value")
        all_hours = _bootstrap_normal_dist(
            n_samples=N_BOOTSTRAP_REPLICATES, mean=hour, variance=30, fun=round
        )
    else:
        logger.debug("Assuming an explicit list of hours from user-specified input")
        all_hours = hour

    _kwargs["mean"] = hour
    _kwargs["variance"] = DOY_VARIANCE_PARAMETER
    _kwargs["fun"] = round
    _kwargs["n_samples"] = N_BOOTSTRAP_REPLICATES

    bs_hourlies = [int(h) for h in unique(_bootstrap_normal_dist(**_kwargs))]
    bs_hourlies = array(bs_hourlies)[array(bs_hourlies) >= 0]
    bs_hourlies = array(bs_hourlies)[array(bs_hourlies) < max_hours]

    logger.debug("Query: x=" + str(x) + "; y=" + str(y) + "; z=" + str(z))

    ret = DataFrame()

    ret[0] = bs_hourlies
    ret[1] = f[dataset][[(z, y, x) for z in bs_hourlies]]

    return ret


def attribute_gdf_w_dataset(
    gdf=None, hour_interval=None, n_bootstrap_replicates=30, dataset=None
):
    """
    """
    logger.debug(
        "Build a target keywords list for our time-series boostrapping" + "procedure"
    )

    _kwargs = dict()

    _kwargs["variance"] = 64
    _kwargs["fun"] = round
    _kwargs["n_samples"] = n_bootstrap_replicates

    logger.debug("Attaching to windtoolkit grid")
    f = h5.File("/nrel/wtk-us.h5", "r", bucket="nrel-pds-hsds")
    WTK_MAX_HOURS = f[dataset].shape[0]

    # assume the user only provided a single scalar value that we should build
    # an hourly time-series with
    all_hours = linspace(
        start=1,
        stop=WTK_MAX_HOURS,
        num=(WTK_MAX_HOURS / hour_interval) + 1,
        dtype="int",
    )

    y_overall = DataFrame(zeros(shape=(len(gdf["id"]), len(all_hours))))

    with tqdm(total=len(y_overall)) as progress:
        for i, row in gdf.iterrows():
            for z in all_hours:
                # print('Processing site :: '+ 'x:' + str(row['x']) + '; y:' + str(row['y']) + '; z:'+ str(z) +';\n')
                site_aggregate = _query_timeseries(f, row["y"], row["x"], z, dataset)
                y_overall.iloc[i, array(all_hours) == z] = _polynomial_ts_estimator(
                    y=site_aggregate[1], x=site_aggregate[0]
                )
                progress.update(i)

    f.close()
    del f

    return y_overall


if __name__ == "__main__":

    datasets = []
    gdf = GeoDataFrame().from_file("vector/h5_grid.shp")

    for dataset in datasets:
        result = attribute_gdf_w_dataset(
            gdf,
            hour_interval=_HOURS_PER_MONTH,
            n_bootstrap_replicates=30,
            dataset=dataset,
        )
