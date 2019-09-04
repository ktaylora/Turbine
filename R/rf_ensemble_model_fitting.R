#!/usr/bin/Rscript
#
# Author : Kyle Taylor (PLJV)
# License : GPL v.3 (2019)
#

#suppressMessages(require(sf))
suppressMessages(require(rgdal))
suppressMessages(require(raster))
suppressMessages(require(pscl))
suppressMessages(require(rjson))
suppressMessages(require(RPostgreSQL))

#
# RUNTIME ARGUMENTS
#

             VERSION   = "19.8.07"
             DEBUGGING = T
MAGNITUDE_OF_IMBALANCE = 2.5 # how many more zeros should we allow, relative to non-zeros?
    VARIANCE_THRESHOLD = 0.94    # to determine how many components to use for ZIP model
               N_TREES = 800
	         N_TREES_MAX = 2000
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

debug <- function(x) debug("DEBUG:", paste(x, collapse=" "), "\n")

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

debug("Reading attributed US National Grid from SQL host:", SESSION$host, "\n")

predict_data <- units <- rgdal::readOGR(
  DSN,
  layer="ktaylora.turbine_training_data",
  stringsAsFactors=F,
  verbose=F
)

training_data <- units@data; 
rm(units);

if (DEBUGGING) cat("Variables available in training dataset:", paste(colnames(training_data), collapse=","), "\n") 

# These data are MASSIVELY zero-inflated -- finding a wind farm
# is a needle-in-a-haystack problem that is going to make fitting a
# meaningful model difficult. Let's limit the zeros to
# by randomly selecting zeros for modeling from the full extent of
# the PLJV region

N_ZEROS_TO_USE     <- round(sum(training_data$t_count > 0)* MAGNITUDE_OF_IMBALANCE)
N_ZEROS_TO_USE     <- sample(which(training_data$t_count == 0), size=N_ZEROS_TO_USE)
N_NON_ZEROS_TO_USE <- which(training_data$t_count > 0)

if(DEBUGGING) {
  cat(
    "Using n1=",
    length(N_NON_ZEROS_TO_USE),
    "Non-zero values and n2=",
    length(N_ZEROS_TO_USE),
    "Zero values for modeling\n"
  )
}
# center our explanatory data. Don't center our response data, though
turbine_counts <- as.numeric(as.vector(training_data[,'t_count']))

# force numeric data type for every column
debug("Forcing numeric values for input training data\n")

for( i in 1:ncol(training_data)) training_data[,i] <- as.numeric(as.vector(training_data[,i]))


debug("Mean-variance centering our data\n")
m_scale <- scale(
  training_data[,!grepl(colnames(training_data), pattern="t_count|id")]
)

# The PCA was used for a ZIP (likelihood) model -- not random forests
# But DO take a look at the m_pca object to see which variables are
# actually contributing unique information to the RF model

debug("Building a PCA\n")

m_pca <- prcomp(as.data.frame(m_scale))

predict_data@data <- as.data.frame(m_pca$scale)

PCA_MAX_COMPONENT_TO_USE <- min(
  which(summary(m_pca)$importance[3,] > VARIANCE_THRESHOLD)
)

debug(
  "Needed ",
  PCA_MAX_COMPONENT_TO_USE,
  "components to capture",
  VARIANCE_THRESHOLD*100,
  "% of the variance of the training data\n"
)

debug("Subseting our original units dataset to deal with zero-inflation\n")

training_data  <- training_data[ c(N_ZEROS_TO_USE, N_NON_ZEROS_TO_USE) , ]
training_data_centered <- as.data.frame(m_scale)[c(N_ZEROS_TO_USE, N_NON_ZEROS_TO_USE),]

# append our un-scaled counts
debug("Appending our response variable to training data\n")

training_data_centered$t_count <- turbine_counts[c(N_ZEROS_TO_USE, N_NON_ZEROS_TO_USE)]
training_data$t_count <- turbine_counts[c(N_ZEROS_TO_USE, N_NON_ZEROS_TO_USE)]

# fit our full and intercept models
debug("Fitting a first round of zip models to centered data\n")

m_zip_full <- zeroinfl(t_count ~ . | ., data=training_data_centered)
m_zip_logit_intercept_only <- zeroinfl(t_count ~ . | +1, data=training_data_centered)
m_zip_intercept_only <- update(m_zip_full, . ~ 1)

#pchisq(2 * (logLik(m_zip_full) - logLik(m_zip_intercept_only)), df = 20, lower.tail = FALSE)
# do a PCA to tease apart the variance (and correlation) amoung our predictors
training_data_pca <- as.data.frame(m_pca$x)[c(N_ZEROS_TO_USE, N_NON_ZEROS_TO_USE),]
training_data_pca$t_count <- turbine_counts[c(N_ZEROS_TO_USE, N_NON_ZEROS_TO_USE)]

# note that we aren't subsetting our components -- throw all of the data at RF and
# let it decide what is useful

debug("Over-fitting an initial RF model and determine convergence for our number-of-trees parameter (this may take an hour or two)\n")

m_rf_regression_pca_full <- randomForest::randomForest(
  as.numeric(t_count) ~ .,
  data=training_data_pca,
  do.trace=F,
  ntree=N_TREES_MAX, 
  importance=FALSE,
  replace=T,
  corr.bias=T
)

oob_err_rate_vs_trees <- round(
  diff(diff(m_rf_regression_pca_full$mse)), 
  5 # just-in-case, only accept 5 decimal digits of precision
)

RF_BURNIN_THRESHOLD <- sd(oob_err_rate_vs_trees) 

# take an arbitrary value in the right tail of the 
# err distribution as our cut-off and re-fit our RF

N_TREES <- as.vector(
  round(
    quantile(
      which(abs(oob_err_rate_vs_trees) > RF_BURNIN_THRESHOLD), 
      probs=0.99)
  )
)

debug("Refitting RF using", N_TREES, "trees to control for over-fit\n");

m_rf_regression_pca_full <- randomForest::randomForest(
  as.numeric(t_count) ~ .,
  data=training_data_pca,
  do.trace=F,
  ntree=N_TREES, 
  importance=TRUE,
  replace=T,
  corr.bias=T
)

if( DEBUGGING ){
  cat("R-squared for RF Model: ");
  cat(round(max(m_rf_regression_pca_full$rsq),2));
  cat("\n");
}

# classification trees can efficiently bin discrete counts and represent zero (and other) 
# ordinal intervals nicely. Let's fit some classification trees forest to zero and non-zero
# classes. We will ensemble the results with our regression trees in a full CART that will
# (hopefully) do a pretty good job of representing zero and non-zero values

if( DEBUGGING ) cat("Over-fitting RF classification trees to determine burn-in\n")

training_data_pca$t_count <- ( as.vector(training_data_pca$t_count) > 0 )

m_rf_classifier_pca_full <- randomForest::randomForest(
  as.factor(t_count) ~ .,
  data=training_data_pca,
  do.trace=F,
  ntree=N_TREES_MAX,
  importance=FALSE
)

oob_err_rate_vs_trees <- round(
  diff(diff(m_rf_classifier_pca_full$err.rate[,1])), 
  5 # just-in-case, only accept 5 decimal digits of precision
)

RF_BURNIN_THRESHOLD <- sd(oob_err_rate_vs_trees) 

# take an arbitrary value in the right tail of the 
# err distribution as our cut-off and re-fit our RF

N_TREES <- as.vector(
  round(
    quantile(
      which(abs(oob_err_rate_vs_trees) > RF_BURNIN_THRESHOLD), 
      probs=0.81)
  )
)

debug("Refitting RF (classifier) using", N_TREES, "trees to control for over-fit\n");

m_rf_classifier_pca_full <- randomForest::randomForest(
  as.factor(t_count) ~ .,
  data=training_data_pca,
  do.trace=F,
  ntree=N_TREES,
  importance=FALSE
)


#
# Fit a ZIP model
#

debug("Fitting a ZIP model for reference using PCA subsetting");

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
debug("ZIP PCA (Full Model) R-squared:",r_squared,"\n")

cat("Staging for region-wide prediction\n");

# df <- predict_data@data[,!grepl(colnames(predict_data@data), pattern="t_count|id")]
# for( i in 1:ncol(df)) df[,i] <- as.numeric(as.vector(df[,i]))
# df <- as.data.frame(scale(df))
# df <- as.data.frame(prcomp(df)$x)
# predict_data@data <- cbind(data.frame(t_count=predict_data$t_count), df); 
# rm(df);

# predict across the full PLJV region using competing models

debug("Predicting RF classifier across the entire study region\n")

rf_turbine_presence_absence <- predict(
  m_rf_classifier_pca_full, 
  newdata=predict_data@data,
  type='response'
)

debug("Predicting PCA model across the entire study region\n")
zip_turbine_density_map <- pscl:::predict.zeroinfl(
  m_zip_pca_full,
  newdata=predict_data@data,
  type='response'
)

debug("DEBUG : Predicting RF regression model across the entire study region\n")
rf_regression_density_map <- predict(
  m_rf_regression_pca_full,
  newdata=predict_data@data,
  type='response'
)

result <- data.frame(
  full_ens_density=round(rowMeans(matrix(c(rf_regression_density_map, zip_turbine_density_map), ncol=2)),2)*as.numeric(rf_turbine_presence_absence),
  zip_ens_density=zip_turbine_density_map*as.numeric(rf_turbine_presence_absence),
  rf_regression_ens_density=rf_regression_density_map*as.numeric(rf_turbine_presence_absence),
  rf_classifer_filter=as.numeric(rf_turbine_presence_absence)
)

result <- round(result, 2); result[result < 0] <- 0

predict_data@data <- result;

debug("Writing ensemble predictions to database:", SESSION$database, "\n")

#sf::st_write(
#  as(units_predicted, "sf"),
#  "local_project_data.gpkg",
#  paste("turbine_predicted_densities_zip_v",VERSION,sep=""), 
#  update=T
#)

if( sum(grepl(dbListTables(connection), pattern="turbine_predicted_densities_rf_ensemble")) > 0 ){ 
  suppressMessages(dbSendStatement(connection, "DROP TABLE wind.turbine_predicted_densities_rf_ensemble CASCADE"))
}

rgdal::writeOGR(
  predict_data,
  dsn=DSN,
  layer="wind.turbine_predicted_densities_rf_ensemble",
  driver="PostgreSQL"
)

dbDisconnect(connection); rm(connection);

debug("Saving R workspace\n")

save.image(
  paste("workflows/turbine_model_fitting_v",VERSION,".rdata",sep=""),
  compress=T
)

quit("no", 0)
