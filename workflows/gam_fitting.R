var_response_plot <- function(m=NULL, var=NULL, rug=T){
  dev.new();
  plot(
    m, 
    which=var, 
    type="l", 
    col="red", 
    lwd=1.5,
    rug=T
  );
  grid();
  grid();
}

m_gam <- gamboost(
  as.factor(response)~
    bmono(raster_115)+
    bbs(ws_100m)+
    btree(raster_230)+
    btree(raster_345)+
    bmono(raster_138)+
    btree(slope)+
    btree(rgh_11_11)+
    btree(raster_69), 
  control = boost_control(center=T), 
  family=Binomial(),
  data=training_data
)

test <- predict(
    m_gam, 
    newdata=training_data, 
    type="response"
  )
  

valid <- !is.na(test)
overall_acc <- 
  sum(as.numeric( test[valid] > 0.5 ) == 
  training_data[valid,]$response) / nrow(training_data[valid,])

predicted <- raster::predict(final_stack, type="response")
  predicted <- round(predicted*100)

cat(" -- writing to disk:")

writeRaster(
    predicted, 
    "predicted_wind_suitability.tif",
    datatype="INT1S", 
    overwrite=T, 
    progress='text'
  )
