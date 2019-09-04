#!/usr/bin/Rscript

#
# Tubrine Tabular Asset Ingestion for PostgreSQL/PostGIS
#
# Author : Kyle Taylor (2019)
# License : GPL v.3 (2019)
#

args <- commandArgs(trailingOnly=T)

suppressMessages(require(rjson))
suppressMessages(require(rgdal))
suppressMessages(require(raster))

#
# LOCAL FUNCTIONS
#

debug <- function(x) if(DEBUGGING) cat("DEBUG:", paste(x, collapse=" "), "\n")

#
# RUNTIME ARGUMENTS
#

PROJECT_DIR = "/home/ktaylora/Incoming/turbine_updates_19.8.07/"
    SESSION = rjson::fromJSON(file="config.json")
      TABLE = 'predicted_density_turbines'

#
# MAIN
#

DSN = paste(
  "PG:dbname='", SESSION$database,
  "' host='", SESSION$host,
  "' user='", SESSION$user,
  "' password='", SESSION$password, 
  "' port='5432'",
  sep=""
)

# attempt to commit to our database
result <- try(rgdal::writeOGR(
  transect_centroids,
  dsn=DSN,
  layer=paste("wind", TABLE, sep="."),
  driver="PostgreSQL",
  check_exists=F,
  overwrite_layer=T,
  delete_dsn=F
))

if(!inherits(result, "try-error")){
  quit("no", status=0)
} else {
  quit("no", status=1)
}
