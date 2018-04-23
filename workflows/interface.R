#!/usr/bin/env R
# Title     : Unix Shell Workflow Interface for 'Turbine'
# Objective : A BASH scriptable interface for 'Turbine' that allows a user to dynamically update and 
# evaluate model predictions using the lastest release of the FAA digital obstruction dataset
# Created by: Kyle Taylor (kyle.taylor@pljv.org)
# Created on: 4/19/18

# Default runtime options
WORKSPACE_DIR = "/gis_data/Wind"

# Load package defaults
stopifnot(require(Turbine))

argv <- commandArgs(trailingOnly=T)

# MAIN

# sanity check : do we have the most recent version of the FAA dataset?
if ( !Turbine:::check_fetch_most_recent_obstruction_file(proposed_zip=Turbine:::web_scrape_faa_digital_obstructions(write=F)) ){
  print("001 : The proposed FAA download isnt newer than what we already have available.")
  stop("001 : The proposed FAA download isnt newer than what we already have available.")
}
# download and read-in the most recent FAA dataset
wind_occurrence_pts <- web_scrape_faa_digital_obstructions()
wind_occurrence_pts <- unpack_faa_zip(wind_occurrence_pts)
# generate our pseudo-absences and merge with out input occurrence points
wind_absence_pts <- gen_pseudo_absences(
    pts=wind_occurrence_pts, 
    boundary=rgdal::readOGR("/gis_data/PLJV/","PLJV_Boundary", verbose=F)
  )
wind_training_pts <- merge_presences_absences_by_year(
    presences=wind_occurrence_pts, 
    absences=wind_absence_pts, 
    bag=F
  )
# read-in our previously built raster explanatory data
explanatory_data <- load_explanatory_data()$explanatory_variables
# extract across our predictor dataset  
training_data <- raster::extract(
    explanatory_data,
    wind_training_pts, 
    df=T
  )
# merge our categorical "response" variable into the training data table 
training_data <- cbind(
    training_data, 
    response=wind_training_data$response
  )
# drop our lurking ID column if it exists
training_data <- training_data[ ,!grepl(tolower(colnames(training_data)), pattern="id") ]
