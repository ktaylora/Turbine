"""
Performs statistical downscaling (using a 'robust linear estimator') of wind data from NREL's ~2 km grid to a 30 meter grid that we
can use to attribute wind data to the 1 KM US National Grid
"""

# generate a 250 meter sub-grid across the project region extent (using USNG units)
# fetch our 30 meter NED raster data for project region extent
# generate the usual suspects for topographic variables
# use a bilinear interpolation to downscale our attributed NREL grid to our 30 meter NED grid
# use a robust regression to estimate NREL station-level wind variable ~ f(bilinear interpolation + topographic variables)
# write 30 meter downscaled raster surfaces to disc for review
# aggregate mean + variance statistics of the downscaled grid to our US National grid units
# write US National Grid units to PostGIS
