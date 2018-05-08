#!/usr/bin/env R
# Title     : Unix Shell Workflow Interface for 'Turbine'
# Objective : A BASH scriptable interface for 'Turbine' that allows a user to dynamically update and
# evaluate model predictions using the lastest release of the FAA digital obstruction dataset
# Created by: Kyle Taylor (kyle.taylor@pljv.org)
# Created on: 4/19/18

# Default runtime options
TMP_PATH      = "/tmp/r_raster_tmp" # make sure this path has a lot of free space available
WORKSPACE_DIR = "/home/ktaylora/Workspace/turbine" # CWD for our run
N_PASS        = 3 # number of smoother passes to use for feature extraction
DATE_STRING   = tolower(gsub(format(Sys.time(), "%b %d %Y"), pattern = " ", replacement = "_")) # today's date

# Load package defaults
stopifnot(require(Turbine))
setwd(WORKSPACE_DIR)

# MAIN

stopifnot(require(randomForest))
stopifnot(require(raster))

raster::rasterOptions(tmpdir=TMP_PATH)

# sanity check : do we have the most recent version of the FAA dataset?
# -- be noisy about it so that we can script this call from Python
if ( !Turbine:::check_fetch_most_recent_obstruction_file(
        proposed_zip=Turbine:::web_scrape_faa_digital_obstructions(write=F)
      )){
  print("001 : The proposed FAA download isnt newer than what we already have available.")
  stop("001 : The proposed FAA download isnt newer than what we already have available.")
}
# download and read-in the most recent FAA dataset
wind_occurrence_pts <- Turbine:::web_scrape_faa_digital_obstructions()
wind_occurrence_pts <- Turbine:::unpack_faa_zip(wind_occurrence_pts)
# define our project region boundary
boundary <- sp::spTransform(rgdal::readOGR(
  "/gis_data/PLJV/","PLJV_Boundary", verbose=F), sp::CRS(raster::projection(wind_occurrence_pts)
))
# drop turbine locations outside of our project area
wind_occurrence_pts <- wind_occurrence_pts[ !is.na(
  sp::over(wind_occurrence_pts, boundary))[,1], ]
# cache our wind turbine dataset to disk
rgdal::writeOGR(
  wind_occurrence_pts,
  dsn=".",
  layer=paste("wind_turbine_pts_", DATE_STRING, sep=""),
  driver="ESRI Shapefile",
  overwrite=T
)
# generate our pseudo-absences and merge with out input occurrence points
wind_absence_pts <- Turbine:::gen_pseudo_absences(
    pts=wind_occurrence_pts,
    boundary=rgdal::readOGR("/gis_data/PLJV/","PLJV_Boundary", verbose=F)
  )
# if we have an existing product to work with, let's build an evaluation dataset from new turbines
# and let's see how well our last raster prediction does on predicting 'occurrence' of the current dataset
previous_suitability_raster <- list.files(
    WORKSPACE_DIR,
    pattern="^predicted_suitability_",
    full.names=T
  )

if ( length(previous_suitability_raster) > 0 ){
  YEAR <- unlist(strsplit(date(), split=" "))
  YEAR <- as.numeric(YEAR[length(YEAR)])
  previous_suitability_raster <- previous_suitability_raster[length(previous_suitability_raster)]
  previous_suitability_raster <- raster::raster(previous_suitability_raster)
  wind_evaluation_pts <- Turbine:::merge_presences_absences_by_year(
    presences=wind_occurrence_pts,
    absences=wind_absence_pts,
    years=YEAR,
    bag=F
  )
  predicted <- as.numeric(raster::extract(
    previous_suitability_raster, wind_evaluation_pts, df=F, sp=F) > 0.5
  )
  # build a confusion matrix
  keep <- !is.na(predicted)
  evaluation_matrix <- caret::confusionMatrix(
    as.factor(predicted[keep]),
    as.factor(wind_evaluation_pts$response[keep]),
    positive=as.character(1)
  )
}
# build a training dataset of all of our points for model fitting
wind_training_pts <- Turbine:::merge_presences_absences_by_year(
    presences=wind_occurrence_pts,
    absences=wind_absence_pts,
    bag=F
  )
# read-in our previously built raster explanatory data
explanatory_data <- Turbine:::load_explanatory_data()$explanatory_variables
# extract across our predictor dataset
wind_training_pts <- sp::spTransform(
    wind_training_pts,
    sp::CRS(raster::projection(explanatory_data))
  )
training_data <- raster::extract(
    explanatory_data,
    wind_training_pts,
    df=T,
    na.rm=T
  )
# merge our categorical "response" variable into the training data table
training_data <- cbind(
    training_data,
    response=wind_training_pts$response
  )
# fit a random forest model
m_rf <- Turbine:::fit_rf(training_data=training_data)
# generate a 0-to-1 wind suitability raster surface and cache to disk
predicted <- Turbine:::gen_rf_suitability_raster(
    m = m_rf,
    explanatory_vars = explanatory_data,
    write = paste("predicted_suitability_rf_", DATE_STRING, ".tif", sep=""),
    quietly = F
  )
# three-pass gaussian filter to make for cleaner feature extraction
# for our build-out simulation
predicted <- Turbine:::extractDensities(
  predicted,
  use_smoother=N_PASS
)
areas <- sapply(
  X=sp::split(predicted, 1:nrow(predicted)),
  FUN=function(x){
    rgeos::gArea(x)
  })
# attribute a build-out area that satisfies the total area requirement of our WITA projections
# return a positive result for our scripted interface
cat("0\n");
