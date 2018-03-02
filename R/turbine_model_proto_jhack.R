# 1. read wind turbine points 
wind_pts <- OpenIMBCR:::readOGRfromPath("/gis_data/Wind/faa_products/wind_turbines_feb_26_2018.shp")

# mask all points outside of the PLJV region
boundary <- OpenIMBCR:::readOGRfromPath("/gis_data/PLJV/PLJV_Boundary.shp")

boundary <- sp::spTransform(boundary, sp::CRS(raster::projection(wind_pts)))
not_overlapping <- is.na(sp::over(wind_pts, boundary)[,1])
wind_pts <- wind_pts[!not_overlapping,]
# subset the turbine dataset into observations from 2008 and beyond.
wind_pts <- wind_pts[wind_pts@data$year>=2008,]

writeOGR(wind_pts,".", "my_wind_pts", driver="ESRI Shapefile", overwrite=T)

# 2. Generate Pseudo-absences

# 2.1 Generate raster of presences
wind_pts_raster <- rasterize(
  wind_pts,
  raster(resolution=(1/111319.9)*30,
         ext=extent(boundary),
         crs=CRS(projection(boundary))),
         field="ors",
         progress='text') > 0 
writeRaster(wind_pts_raster, "wind_pts_raster", format = "GTiff", datatype = "LOG1S", 
            overwrite=T, progress = "text")
# LOG1S not working so R reverts datatype to INT1S. Step takes long time!

# 2.2 Generate/cull pseudo-absences
absences <- sampleRandom(wind_pts_raster,size=round(nrow(wind_pts)*3),sp=T,na.rm=F)
# I multiply my wind_pts by 3 to create some wiggle room.
# The line below takes out points within a 1,000 meter radius of wind_pts (turbines)
absences <- absences[is.na(sp::over(absences,rgeos::gBuffer(wind_pts,byid=T,width=(1/111319.9)*1000))[,1]),]
# This line takes out points that lie within wind_points_raster but outside the PLJV boundary.
# The two extents are not contiguous; wind_points_raster is a square as large as the PLJV boundary
# at its widest and tallest points. We might consider masking out both the non-pljv area and turbine
# area before we generate absences at all so we don't have to allow for wiggle room, because I
# guess there is a small chance a 3x sample still wouldn't be big enough to ensure there are as many
# pseudo-absences as turbines.
absences <- absences[as.vector(!is.na(sp::over(absences,boundary)[,1])),]
if(nrow(absences) > nrow(wind_pts)){
  absences <- absences[sample(1:nrow(absences), size=nrow(wind_pts)),] 
}

# 2.3 combine turbines and pseudo-absences into a dataframe
names(wind_pts) <- names(absences) <- "layer"
training_data <- rbind(wind_pts[,1], absences)
training_data@data <- data.frame(response=c(rep(1,nrow(wind_pts)),rep(0,nrow(absences))))

writeOGR(training_data,"presence_absence_vector", "training data", driver="ESRI Shapefile",overwrite=T)

# 3.1 Extract topographic predictor variables:
topographic_variables <-
  c('aspect.tif',
    'slope.tif',
    'topographic_roughness_3x3.tif',
    'topographic_roughness_11x11.tif',
    'topographic_roughness_33x33.tif'
  )
# line below creates a raster stack of all topographic rasters
stack <- stack(mapply(raster,topographic_variables))

t <- extract(stack,training_data,df=T)
# extract will extract datapoints from each raster in the stack for all ~ 32,000
# training points. How to validate that points were drawn from the same spatial location?
t$response <- training_data@data[,1]
# for some reason response was lost during the extraction process and had to be added back in.

# 3.2 Extract wind speed variables:
# Right now I only have data for windspeed at 100m. It looks like the average turbine height
# is around 300 ft, right between 100 and 80m, which NREL also has a map for. We might want to
# interpolate an 80m raster too, although it does look like the two maps are highly correllated.
setwd("/home/jhack/global_workspace/Joe_Work/turbine_workspace")
wind_speed_100m <- raster("ws_30m")
# Extract wind_speed_100m points. I will rewrite extraction steps later so that all predictor 
# variables in the topographic, windspeed and transmission categories are loaded at once.
f <- extract(wind_speed_100m, training_data, df=T)
# Merge f dataframe with t
t$wind_speed_100m <- f[,2]

# 3.3 Extract transmission data variables:
transmission_variables <-
  c("raster_69",
    "raster_115",
    "raster_138",
    "raster_230",
    "raster_345"
    )
q <- extract(stack(mapply(raster,transmission_variables)),training_data,df=T)
t$raster_69 <- q[,2]
t$raster_115 <- q[,3]
t$raster_138 <- q[,4]
t$raster_230 <- q[,5]
t$raster_345 <- q[,6]

# 4. Run the model!:

# 4.1: variable selection:
# 

t <- na.omit(t)
model_sel <- rf.modelSel(ydata=as.factor(t$response),
                         xdata=t[,!grepl(names(t),pattern="response|ID")],
                         parsimony=0.15,
                         imp.scale="se",
                         r=seq(0,1,0.05),
                         seed=25,
                         final.model=T)

varImpPlot(model_sel$rf.final)

# the command above will run the final model, but I believe Kyle wants to build many different forests.
# to resolve: cross-validation.