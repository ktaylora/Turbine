#!/usr/bin/env R
# Title     : Model builder interface for 'Turbine'
# Objective : Accept input FAA 'occurrence' data and generate pseudo-absences -- optionally fit a number of models
# to input data, including Random Forests, GAM's, and Maxent
# Created by: Kyle Taylor (kyle.taylor@pljv.org)
# Created on: 4/19/18

#' function that will merge presences and absences SpatialPoints data.frame using a 'year' attribute
#' @export
merge_presences_absences_by_year <- function(presences=NULL, absences=NULL, years=NULL, bag=T){
    if(is.null(years)){
        years <- unique(presences@data[,1])
    }
    presences  <- presences[presences$year %in% years, ]
    if(bag){
      presences <- presences[
          sample(1:nrow(presences), replace=T, size=nrow(presences)),
        ]
      absences <- absences[
        sample(1:nrow(absences), size=nrow(presences), replace=T),
      ]
    } else {
      absences <- absences[
        sample(1:nrow(absences), size=nrow(presences), replace=F),
      ]
    }
    # merge and return to user
    presences@data <- data.frame(response=rep(1, nrow(presences)))
    absences@data <- data.frame(response=rep(0, nrow(absences)))
    return(sp::rbind.SpatialPointsDataFrame(presences,absences))
}
#' generate a large pool of random pseudo-absences within a geographic boundary
#' using a set of user-specified wind turbine point locations. The pool will be
#' arbitrarily large (nrow(pts)*iter) with the intention that the dataset will
#' be downsampled through a bagging procedure down-the-line.
#' @export
gen_pseudo_absences <- function(
  pts=NULL,
  boundary=NULL,
  iter=5,
  buffer_width=1000,
  quietly=T
){
  # let's reproject our data to a planar geography for rgeos
  original_crs <- sp::CRS(raster::projection(pts))
  pts <- sp::spTransform(
      pts,
      sp::CRS(raster::projection("+init=epsg:2163"))
    )
  boundary <- sp::spTransform(
      boundary,
      sp::CRS(raster::projection("+init=epsg:2163"))
    )
  # buffer our turbine occurrences by some user-specified distance
  buffered_turbine_occurrences <- rgeos::gBuffer(
      pts, byid=T, width=buffer_width
    )
  if(!quietly) cat(" -- generating a large pool of pseudo-absences\n")
  wind_pts_raster <- raster::raster(
      resolution=30,
      ext=raster::extent(boundary),
      crs=sp::CRS(raster::projection(boundary))
    )
  # set our pool to zero
  raster::values(wind_pts_raster) <- 0
  # randomly generate a pool of pseudo-absences across n=iter steps
  if(!quietly) cat(" -- iteritively re-sampling absence space:")
  pseudo_absences <- do.call(sp::rbind.SpatialPointsDataFrame, lapply(
    X=1:iter,
    FUN=function(r){
      if(!quietly) cat(".")
      absences <- raster::sampleRandom(
          wind_pts_raster,
          size=round(nrow(pts)*1.2),
          sp=T,
          na.rm=F
        )
      absences <- absences[
        is.na(sp::over(absences, buffered_turbine_occurrences)[,1]),]
      return(absences)
    }
  ))
  if(!quietly) cat("\n")
  pseudo_absences@data <- data.frame(response=rep(0, nrow(pseudo_absences)))
  # reproject back to our native CRS and return to user
  return(sp::spTransform(pseudo_absences, original_crs))
}
#' load an input R data file containing a raster stack and (optionally) our
#' input training data
#' @export
load_explanatory_data <- function(path="."){
  model_fitting_data <- new.env()
  # are we a directory?
  if(dir.exists(path)){
     dir <- path
    path <- list.files(path, pattern="[.]rdata$", full.names=T)
  } else if(file.exists(path)){
     dir <- paste(dir[1:(length(dir)-1)], collapse="/")
    path <- path
  } else {
    stop("couldn't find or use any rdata files specified by path= argument")
  }
  # path could be a single rdata file, or a number of rdata files
  if(length(path) > 0){
    catch <- sapply(
        X=1:length(path),
        FUN=function(i)(
          return(class(try(load(path[i], envir=model_fitting_data))))
        )
    )
    if("try-error" %in% catch) stop("failed to load an input rdata file -- this shouldn't happen")
    # make a list out of our data and return to user
    ret <- list()
    if ( sum(grepl(ls(envir=model_fitting_data), pattern="wind_pts")) > 0 ){
      ret$training_data <- get("wind_pts", envir=model_fitting_data)
    }
    ret$explanatory_variables <- raster::stack(paste(dir,"explanatory_vars.tif",sep="/"))
    names(ret$explanatory_variables) <- get("n", envir=model_fitting_data)
    return(ret)
  }
}
#' fit a generalized additive model to a user-specified dataset using an optional functional specification
#' @export
fit_boosted_gam <- function(
  formula=NULL,
  training_data=NULL,
  vars=NULL,
  control=mboost::boost_control(center=T)
  ){
  stopifnot(require(mboost))
  if (is.null(vars)) {
      vars <- unique(colnames(as.data.frame(training_data)))
  }
  if (is.null(formula)) {
      # if the user didn't specify a formula, we will use the default b-spline base learners from 'mboost'
      m <- mboost::gamboost(
        formula=as.factor(response)~.,
        data=as.data.frame(training_data)[,vars],
        family=mboost::Binomial(),
        control = control,
      )
  } else {
      m <- mboost::gamboost(
        formula=as.formula(formula),
        data=as.data.frame(training_data)[,vars],
        family=mboost::Binomial(),
        control = control,
      )
  }
  # return fitted model object to user
  return(m)
}
#' fit a random forests model to the input dataset
#' @export
fit_rf <- function(training_data=NULL, max_trees=2000){
  training_data <- na.omit(training_data)

  N_TREES <- which.min(c(max_trees, ceiling(nrow(training_data)/2)))

  if(!require(randomForest)){
    install.packages('randomForest', repos="https://cran.revolutionanalytics.com")
    require(randomForest)
  }

  m <- randomForest::randomForest(
    y=as.factor(training_data$response),
    x=training_data[,!grepl(tolower(colnames(training_data)), pattern="response|id")],
    data=training_data,
    na.rm=T,
    ntree=N_TREES
  )

  return(m)
}
#' fit a maxent model to the input dataset
#' @export
#' training_data is only two columns: 1: longitude and 2: latitude, which correspond to the locations of turbines.
#' predictor_stack is a stack of all the predictors from the study area that will be sampled by the maxent function to create absences.
fit_maxent <- function(formula=NULL, training_data=NULL, predictor_stack=NULL){
  # install and load dismo, the ecology maxent package
  if(!require(dismo)){
    install.packages('dismo', repos="https://cran.revolutionanalytics.com")
    require(dismo)
  }
  # section on making sure java package is installed and the maxent.jar is in the right folder.
  file.copy(
  from=paste(paste(system.file(package="Turbine"), "data", sep="/"), "/maxent.jar", sep=""),
  to=system.file("java", package="dismo", overwrite=T)
  )
  # Model fitting:
  m <- maxent(x = predictor_stack, p = training_data, nbg = nrow(as.data.frame(training_data)), progress = 'text')
  return(m)
}
