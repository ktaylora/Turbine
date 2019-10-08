"""
Manage downloading and ingesting tasks for NREL's HSDS Service
"""

import warnings
warnings.filterwarnings('ignore')

from shapely.geometry import Point
from geopandas import GeoSeries, GeoDataFrame

import h5pyd as h5

def extent_to_hsds_bounding_box(extent=None):
    ne = (1,1)
    sw = (2,2)
    return (ne, sw)

f = h5.File("/nrel/wtk-us.h5", 'r')

for d in datasets:
    ds = f[d][tmin:tmax:tskip,bd[1][0]:bd[0][0]:latskip,bd[1][1]:bd[0][1]:lonskip]
    lf[d] = ds

test = f['coordinates'][100:]
o=GeoDataFrame({'geometry':GeoSeries([Point(g) for g in test[1]])})

