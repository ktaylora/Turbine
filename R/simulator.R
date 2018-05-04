# Title     : Wind energy build-out simulator interface for 'Turbine'
# Objective : Accepts a suitability raster surface and some scalar estimate of build out to simulate wind farm placement
# Created by: Kyle Taylor (kyle.taylor@pljv.org)
# Created on: 4/19/18

#' Downloads a managed CSV file that contains projections of wind
#' energy build-out and extracts 1.) total area of build-out and
#' 2.) number of turbines associated with WITO projections
scrape_turbine_buildout <- function(){
    URL = paste("https://docs.google.com/spreadsheets/d/e/2PACX-1vT1k16zT14aPbe9l",
    "w2SFMKrwdnd83crKe8wsMoyvjb0YI8t0y7_nTJ7vd7IEFwHEZGkQv6YFX15-Ckd/pub?gid=0&si",
    "ngle=true&output=csv", sep="")
    unlink("fetch.csv", force=T)
    download.file(URL, destfile="fetch.csv", quiet=T)
    t <- read.csv("fetch.csv", header=F)
    # the first element is the number of turbines built, the second
    # is the predicted total area of build-out
    y_2018 <- t[(nrow(t)-1):(nrow(t)),2]
    y_2018 <- as.numeric(gsub(matrix(y_2018, ncol=1), pattern=",|[.]", replacement=""))
    y_2020 <- t[(nrow(t)-1):(nrow(t)),3]
    y_2020 <- as.numeric(gsub(matrix(y_2020, ncol=1), pattern=",|[.]", replacement=""))
    y_2030 <- t[(nrow(t)-1):(nrow(t)),4]
    y_2030 <- as.numeric(gsub(matrix(y_2030, ncol=1), pattern=",|[.]", replacement=""))
    y_2050 <- t[(nrow(t)-1):(nrow(t)),5]
    y_2050 <- as.numeric(gsub(matrix(y_2050, ncol=1), pattern=",|[.]", replacement=""))
    t <- cbind(y_2018, y_2020, y_2030, y_2050)
    rownames(t) <- c("turbines built", "total area of buildout")
    unlink("fetch.csv")
    return(t)
}
#' Generates a uniform grid of polygons such that each grid 'cell' is
#' a user-specified size. This follows the original methodology of
#' A. Daniels (2017) used for national wildlife refuge planning (LCD).
generate_uniform_grid <- function(boundary=NULL, cell_size=9397.127){
  # generate a uniform grid of points using the 'sp' package
  grid <- sp::makegrid(boundary, cellsize=cell_size)
  grid <- sp::SpatialPointsDataFrame( grid, data=data.frame(id=1:nrow(grid)) )
  # define our CRS to whatever the boundary object uses
  raster::projection(grid) <- raster::projection(boundary)
  # mask out units outside of our boundary
  grid <- rgeos::gBuffer(grid, width=cell_size/2, capStyle='square', byid=T)
  grid <- grid[ !is.na(sp::over(grid, boundary)[,1]) , ]

  return(grid)
}
#'
#'
calc_attribute_turbines <- function(wind_pts=NULL, turbines_grid=NULL){
   turbines_grid <- sp::spTransform(turbines_grid, sp::CRS(raster::projection(wind_pts)))
   # rows here correspond to each turbine grid unit
   turbines_grid$turbines <- apply(
     rgeos::gIntersects(wind_pts, turbines_grid, byid=T),
     MARGIN=1,
     FUN=sum
   )
  return(turbines_grid)
}
#'
#'
calc_suitability_buildout <- function(total_turbines=NULL, turbines_per_farm=NULL, turbines_grid=NULL, suitability_surface=NULL){
  # calculate mean suitability per
  return(NULL)
}
