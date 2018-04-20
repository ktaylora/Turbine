# Title     : Wind energy build-out simulator interface for 'Turbine'
# Objective : Accepts a suitability raster surface and some scalar estimate of build out to simulate wind farm placement
# Created by: Kyle Taylor (kyle.taylor@pljv.org)
# Created on: 4/19/18

#'
#'
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
