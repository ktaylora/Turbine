#!/usr/bin/Rscript
#
# Author : Kyle Taylor (PLJV)
# License : GPL v.3 (2019)
#

suppressMessages(require(sf))
suppressMessages(require(rgdal))
suppressMessages(require(raster))
suppressMessages(require(pscl))
suppressMessages(require(rjson))
suppressMessages(require(RPostgreSQL))

#
# RUNTIME ARGUMENTS
#

             VERSION   = "19.5.1"
             DEBUGGING = T
MAGNITUDE_OF_IMBALANCE = 2.5 # how many more zeros should we allow, relative to non-zeros?
    VARIANCE_THRESHOLD = 0.94    # to determine how many components to use for ZIP model
               N_TREES = 800
               SESSION = rjson::fromJSON(file="config.json")
                 TABLE = 'turbine_predicted_densities'

# Database URI for PostGIS / Rgdal
DSN = paste(
  "PG:dbname='", SESSION$database,
  "' host='", SESSION$host,
  "' user='", SESSION$user,
  "' password='", SESSION$password, 
  "' port='", SESSION$port, 
  "'",
  sep=""
)

# Database configuration for tabular data
connection <- RPostgreSQL::dbConnect(
  dbDriver("PostgreSQL"), 
  dbname = SESSION$database,
  host = SESSION$host, 
  port = SESSION$port,
  user = SESSION$user, 
  password = SESSION$password
)

cat(paste("\n## Turbine Wind Suitability Model (v.",VERSION,")\n\n",sep=""))

#
# LOCAL FUNCTIONS
#

#' calculate a deviance statistic from a count (Poisson) model,
#' e.g., : https://goo.gl/KdEUUa
#' @export
est_deviance <- function(m, method="residuals"){
  # by default, use model residual error to estimate deviance
  if (grepl(tolower(method), pattern = "resid")) {
    if ( inherits(m, "unmarkedFit") ) {
      observed <- unmarked::getY(m@data)
      expected <- unmarked::fitted(m)
    } else {
      observed <- m$y
      expected <- fitted(m)
    }
    # Deviance of full model : 2*sum(obs*log(obs/predicted)-(obs-predicted))
    dev.part <- ( observed * log(observed/expected) ) - (observed - expected)
    sum.dev.part <- sum(dev.part,na.rm=T)
    dev.sum <- 2*sum.dev.part
    return(dev.sum)
  # alternatively, use the log-likelihood value returned from our optimization
  # from unmarked. This approach is advocated by B. Bolker, but the likelihood
  # values returned from hierarchical models can look strange. Use this
  # with caution.
  } else if (grepl(tolower(method), pattern = "likelihood")) {
    if ( inherits(m, "unmarkedFit") ) {
      return(2*as.numeric(abs(m@negLogLik)))
    } else if ( inherits(m, "glm") ) {
      return(-2*as.numeric(logLik(m)))
    } else {
      return(NA)
    }
  }
}

#
# MAIN
#

if(DEBUGGING) cat("DEBUG: READING ATTRIBUTED US National Grid\n")

units <- rgdal::readOGR(
  DSN,
  layer="grids.usng_units_1km",
  stringsAsFactors=F,
  verbose=F
)

training_data  <- RPostgreSQL::dbReadTable(
    connection,
    c("ktaylora","turbine_training_data")
)

# These data are MASSIVELY zero-inflated -- finding a wind farm
# is a needle-in-a-haystack problem that is going to make fitting a
# meaningful model difficult. Let's limit the zeros to
# by randomly selecting zeros for modeling from the full extent of
# the PLJV region.
N_ZEROS_TO_USE     <- round(sum(training_data$t_count > 0)* MAGNITUDE_OF_IMBALANCE)
N_ZEROS_TO_USE     <- sample(which(training_data$t_count == 0), size=N_ZEROS_TO_USE)
N_NON_ZEROS_TO_USE <- which(training_data$t_count > 0)

if(DEBUGGING) {
  cat(
    "DEBUG: Using n1=",
    length(N_NON_ZEROS_TO_USE),
    "Non-zero values and n2=",
    length(N_ZEROS_TO_USE),
    "Zero values for modeling\n"
  )
}
# center our explanatory data. Don't center our response data, though
turbine_counts <- as.numeric(as.vector(training_data[,'t_count']))
# force numeric data type for every column
if(DEBUGGING) cat("DEBUG: Forcing numeric values for input training data\n")

for( i in 1:ncol(training_data)){
  training_data[,i] <- as.numeric(as.vector(training_data[,i]))
}

if(DEBUGGING) cat("DEBUG: Mean-variance centering our data\n")
m_scale <- scale(
  training_data[,!grepl(colnames(training_data), pattern="t_count|id")]
)

if(DEBUGGING) cat("DEBUG: Doing a principal components analysis\n")

m_pca <- prcomp(as.data.frame(m_scale))

PCA_MAX_COMPONENT_TO_USE <- min(
  which(summary(m_pca)$importance[3,] > VARIANCE_THRESHOLD)
)

if(DEBUGGING) cat(
  "DEBUG: Needed ",
  PCA_MAX_COMPONENT_TO_USE,
  "components to capture",
  VARIANCE_THRESHOLD*100,
  "% of the variance of the training data\n"
)
if(DEBUGGING) cat("DEBUG: Subseting our original units dataset to deal with zero-inflation\n")

training_data  <- training_data[ c(N_ZEROS_TO_USE, N_NON_ZEROS_TO_USE) , ]
training_data_centered <- as.data.frame(m_scale)[c(N_ZEROS_TO_USE, N_NON_ZEROS_TO_USE),]

# append our un-scaled counts
if(DEBUGGING) cat("DEBUG: Appending our response variable to training data\n")

training_data_centered$t_count <- turbine_counts[c(N_ZEROS_TO_USE, N_NON_ZEROS_TO_USE)]
training_data$t_count <- turbine_counts[c(N_ZEROS_TO_USE, N_NON_ZEROS_TO_USE)]

# fit our full and intercept models
if(DEBUGGING) cat("DEBUG: Fitting a first round of zip models to centered data\n")

m_zip_full <- zeroinfl(t_count ~ . | ., data=training_data)
m_zip_logit_intercept_only <- zeroinfl(t_count ~ . | +1, data=training_data)
m_zip_intercept_only <- update(m_zip_full, . ~ 1)

#pchisq(2 * (logLik(m_zip_full) - logLik(m_zip_intercept_only)), df = 20, lower.tail = FALSE)
# do a PCA to tease apart the variance (and correlation) amoung our predictors
training_data_pca <- as.data.frame(m_pca$x)[c(N_ZEROS_TO_USE, N_NON_ZEROS_TO_USE),]
training_data_pca$t_count <- turbine_counts[c(N_ZEROS_TO_USE, N_NON_ZEROS_TO_USE)]

# fit a zip to an optimal subset of our PCA scores
if(DEBUGGING) cat("DEBUG: Fitting zip models to a subset of the PCA scores matrix\n")

# note that we aren't subsetting components -- throw all of the data at RF and
# let it decide what is useful
m_rf_regression_pca_full <- randomForest::randomForest(
  as.numeric(t_count) ~ .,
  data=training_data_pca,
  do.trace=T,
  ntree=N_TREES, 
  importance=TRUE,
  replace=T,
  corr.bias=T
)

training_data_pca$t_count <- cut(
  as.numeric(as.vector(training_data_pca$t_count)),
  breaks=c(0,1,3,12),
  include.lowest=T,
  labels=c(0,3,10)
)

m_rf_classifier_pca_full <- randomForest::randomForest(
  as.factor(t_count) ~ .,
  data=training_data_pca,
  do.trace=T,
  ntree=N_TREES,
  importance=TRUE
)

training_data_pca$t_count <- as.numeric(turbine_counts[c(N_ZEROS_TO_USE, N_NON_ZEROS_TO_USE)])

m_zip_pca_full <- zeroinfl(
  t_count ~ . | .,
  data= training_data_pca[,c(1:PCA_MAX_COMPONENT_TO_USE, ncol(training_data_pca))]
)

m_zip_pca_logit_intercept_only <- zeroinfl(
  t_count ~ . | +1,
  data= training_data_pca[,c(1:PCA_MAX_COMPONENT_TO_USE, ncol(training_data_pca))]
)

m_zip_pca_intercept_only <- update(m_zip_pca_full, . ~ 1)

# zip model pseudo r-squared values
r_squared <- round( (
  est_deviance(m_zip_pca_intercept_only) - est_deviance(m_zip_pca_full)
) / est_deviance(m_zip_pca_intercept_only) , 3)
if(DEBUGGING) cat("DEBUG: ZIP PCA (Full Model) R-squared:",r_squared,"\n")

# predict across the full PLJV region using competing  models
if(DEBUGGING) cat("DEBUG: Predicting PCA model across the entire study region\n")
turbine_density_map <- pscl:::predict.zeroinfl(
  m_zip_pca_full,
  newdata=as.data.frame(m_pca$x),
  type='response'
)

units_predicted <- units
units_predicted@data <- data.frame(
  density=round(turbine_density_map)
)

if(DEBUGGING) cat("DEBUG: writing ZIP prediction to database:", SESSION$database, "\n")

sf::st_write(
  as(units_predicted, "sf"),
  "local_project_data.gpkg",
  paste("turbine_predicted_densities_zip_v",VERSION,sep=""), 
  update=T
)

turbine_density_map <- predict(
  m_rf_classifier_pca_full,
  newdata=as.data.frame(m_pca$x),
  type='response'
)

units_predicted@data <- data.frame(
  density=as.numeric(as.vector(turbine_density_map))
)

if(DEBUGGING) cat("DEBUG: writing RF prediction to Local Database\n")

sf::st_write(
  as(units_predicted, "sf"),
  "local_project_data.gpkg",
  paste("turbine_predicted_densities_rf_classifier_v",VERSION,sep=""), 
  update=T,
  OVERWRITE=T
)

rf_ensemble <- as.numeric(as.vector(turbine_density_map))

turbine_density_map <- round(as.numeric(predict(
  m_rf_regression_pca_full,
  newdata=as.data.frame(m_pca$x),
  type='response'
)))

units_predicted@data <- data.frame(
  density=round(as.numeric(as.vector(turbine_density_map)))
)

if(DEBUGGING) cat("DEBUG: writing RF prediction to disc as shapefile\n")

sf::st_write(
  as(units_predicted, "sf"),
  "local_project_data.gpkg",
  paste("turbine_predicted_densities_rf_regression_v",VERSION,sep=""), 
  update=T
)

if(DEBUGGING) cat("DEBUG: Building an RF prediction ensemble and commiting to SQL Server\n")

rf_ensemble <- rowMeans(
  cbind(
    as.numeric(as.vector(rf_ensemble)),
    as.numeric(as.vector(turbine_density_map))
  )
)

result <- data.frame(
  id=units$id,
  density=round(as.numeric(as.vector(rf_ensemble)))
)

if( !dbWriteTable(
  connection, 
  c("wind", "turbine_predicted_densities_rf_ensemble"), # note the goofy vector notation for schema/table
  value = result, 
  row.names = FALSE,
  append=T
) ) {
  warning(
    "DEBUG: failed writing to PostgreSQL database",
    SESSION$database,
    "turbine_predicted_densities_rf_ensemble"
  )
} else {
  dbDisconnect(connection);
  rm(connection);
}

if(DEBUGGING) cat("DEBUG: Saving R workspace\n")

save.image(
  paste("turbine_model_fitting_v",VERSION,".RData",sep=""),
  compress=T
)
