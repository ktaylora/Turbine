import warnings

warnings.filterwarnings("ignore")

import logging

logger = logging.getLogger(__name__)

import sys, os
import time

try:
    from geopandas import GeoSeries, GeoDataFrame
except ModuleNotFoundError:
    raise ModuleNotFoundError(
        "Try running: "
        + "!{sys.executable} -m pip install geopandas --upgrade"
        + "in your python shell and run this again."
    )

from shapely.geometry import Point
from shapely.ops import cascaded_union

from pandas import DataFrame

from scipy.stats import norm

from numpy import poly1d as polynomial_regression
from numpy import unique, polyfit, mean, linspace, array, zeros, float64

try:
    import h5pyd as h5
except ModuleNotFoundError:
    raise ModuleNotFoundError(
        "Try running: "
        + "!{sys.executable} -m pip install git+http://github.com/HDFGroup/h5pyd.git --upgrade "
        + "in your python shell and run this again."
    )

import math

try:
    from tqdm import tqdm
except ModuleNotFoundError:
    raise ModuleNotFoundError(
        "Try running: "
        + "!{sys.executable} -m pip install tqdm --upgrade "
        + "in your python shell and run this again."
    )

_WIND_TOOLKIT_DEFAULT_EPSG = "+init=epsg:4326"
_HOURS_PER_MONTH = 730
_HSDS_CACHE_FILE_PATH = "vector/h5_grid.shp"

N_BOOTSTRAP_REPLICATES = 30
DOY_VARIANCE_PARAMETER = 60  # std. dev. parameter for days used for bootstrapping

f = h5.File("/nrel/wtk-us.h5", "r", bucket="nrel-pds-hsds")

_wtk_datasets = [
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
    "windspeed_80m",
]


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
    global f

    if os.path.exists(cache_file):
        logger.debug(
            "Using cached file for project region coordinates for the wind toolkit"
        )

        gdf = GeoDataFrame().from_file(cache_file)
        target_rows = gdf["id"]

    else:

        logger.debug("Fetching coordinates from wind toolkit HSDS interface")

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

        intersection = []
        with tqdm(total=len(gdf.geometry)) as progress:
            for g in gdf.geometry:
                intersection.append(filter_by_intersection.contains(g))
                progress.update(1)
            progress.close()

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

    if r_squared > 0 and r_squared < 0.1:
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

    return round(mean(fitted), 2)


def _query_timeseries(hour=None, x=None, y=None, dataset=None, max_hours=None):
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
    global f
    if not isinstance(hour, list):
        logger.debug("Building a sequence from user-specified single scalar value")
        all_hours = _bootstrap_normal_dist(
            n_samples=N_BOOTSTRAP_REPLICATES, mean=hour, variance=30, fun=round
        )
    else:
        logger.debug("Assuming an explicit list of hours from user-specified input")
        all_hours = hour

    _kwargs = {}

    _kwargs["mean"] = hour
    _kwargs["variance"] = DOY_VARIANCE_PARAMETER
    _kwargs["fun"] = round
    _kwargs["n_samples"] = N_BOOTSTRAP_REPLICATES

    bs_hourlies = [int(h) for h in unique(_bootstrap_normal_dist(**_kwargs))]
    bs_hourlies = array(bs_hourlies)[array(bs_hourlies) >= 0]
    bs_hourlies = array(bs_hourlies)[array(bs_hourlies) < max_hours]

    logger.debug("Query: z=" + str(hour) + "; x=" + str(x) + "; y=" + str(y))

    ret = DataFrame()

    ret[0] = bs_hourlies
    try:
        ret[1] = f[dataset][[(z, x, y) for z in bs_hourlies]]
    except OSError:
        logger.debug("Dropped our connection -- picking up where we left off")
        f.close()
        time.sleep(10)
        globals()["f"] = h5.File("/nrel/wtk-us.h5", "r", bucket="nrel-pds-hsds")
        return _query_timeseries(hour, x, y, dataset, max_hours)
    return ret


def attribute_gdf_w_dataset(
    gdf=None, hour_interval=None, n_bootstrap_replicates=30, dataset=None
):
    """
    Accepts a user-specified H5 GeoDataFrame and attributes points in the grid
    with fitted values of a polynomial time-series regression for some hourly
    time slice. The hour interval is calculated for a full WTK dataset and each
    hourly is bootstrapped (without replacement) to fit the regression estimator.

    In plain english, this is used to estimate 'monthly' or 'weekly' values of
    windspeed (or other variables) from wind toolkit hourly data in a way that
    is a little more robust than just taking the point values for each time-period.

    :param gdf:
    :param hour_interval:
    :param n_bootstrap_replicates:
    :param dataset:
    :return:
    """
    global f

    logger.debug(
        "Build a target keywords list for our time-series boostrapping" + "procedure"
    )

    _kwargs = dict()

    _kwargs["variance"] = 64
    _kwargs["fun"] = round
    _kwargs["n_samples"] = n_bootstrap_replicates

    WTK_MAX_HOURS, y, x = f[dataset].shape

    if max(gdf["y"]) > y:
        raise ValueError(
            "y coordinate in gdf is outside of the range of our" + "target dataset"
        )
    if max(gdf["x"]) > x:
        raise ValueError(
            "x coordinate in gdf is outside of the range of our" + "target dataset"
        )

    # assume the user only provided a single scalar value that we should build
    # an hourly time-series with
    all_hours = linspace(
        start=1,
        stop=WTK_MAX_HOURS,
        num=(WTK_MAX_HOURS / hour_interval) + 1,
        dtype="int",
    )
    # pre-allocate a destination table with zeros -- if this fails
    # due to memory limitations it's better to discover it now
    # rather than while we are querying the HSDS interface
    y_overall = DataFrame(zeros(shape=(len(gdf), len(all_hours))))
    j = 0
    with tqdm(total=len(gdf)) as progress:
        for i, row in gdf.iterrows():
            for z in all_hours:
                # sample our dataset of interest using time-series boostrapping around
                # hour z
                site_aggregate = _query_timeseries(
                    z, row["x"], row["y"], dataset, WTK_MAX_HOURS
                )
                # fit a quadratic regression for value ~f(hour+hour^2)
                fitted = _polynomial_ts_estimator(
                    y=site_aggregate[1], x=site_aggregate[0]
                )
                y_overall.iloc[j, array(all_hours) == z] = fitted
            j += 1
            # how YOU doin'?
            progress.update(j)
            if j > len(y_overall):
                break

    # clean-up our session -- don't let the socket sit there lurking
    f.close()

    return y_overall


if __name__ == "__main__":

    os.chdir("/home/jovyan")
    # os.chdir("/home/ktaylora/Incoming/Turbine")

    datasets = [
        "windspeed_80m",
        "windspeed_100m",
        "windspeed_200m",
        "winddirection_80m",
        "winddirection_100m",
        "winddirection_200m",
    ]

    gdf = GeoDataFrame().from_file("vector/h5_grid.shp")

    for dataset in datasets:
        logger.debug(
            "Iterating over our target WTK datasets and flush the"
            + "attributed results to disc as we go"
        )
        if os.path.exists("vector/" + dataset + ".shp"):
            print("Existing file found (skipping): " + "vector/" + dataset + ".shp")
        else:
            result = gdf.join(
                attribute_gdf_w_dataset(
                    gdf,
                    hour_interval=_HOURS_PER_MONTH,
                    n_bootstrap_replicates=30,
                    dataset=dataset,
                ),
                rsuffix=dataset,
            )

            result.to_file("vector/" + dataset + ".shp")

            del result
