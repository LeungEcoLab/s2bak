#' @title Make predictions using fitted SO, S2, BaK or S2BaK class models
#'
#' #' @description S3 methods for making models predictions using the output
#' from \link[s2bak]{fit.s2bak}, \link[s2bak]{fit.s2bak.so},
#' \link[s2bak]{fit.s2bak.s2} and \link[s2bak]{fit.s2bak.bak}.
#'
#' If the complete model of class `s2bak` is provided, the user can specify
#' what type of model predictions to use to allow for specific outputs from
#' the sub-components of the model. When applying `predict.s2bak`,
#' the function assumes that `predict.fun` is the same for all models.
#'
#' @param model Outputted model from fit.s2bak.so, fit.s2bak.s2, fit.s2bak or
#' s2bak.combine.
#' @param newdata Environmental data to predict
#' @param trait Trait data, optional if output = "s2" or output = "so"
#' @param predict.fun Prediction function for SDM, which must match the model
#' function used for s2bak.s2 and s2bak.so models).
#' @param output The choice of how predictions are made using the `s2bak` model.
#'
#' Only one type of output can be selected: sightings-only (output = "so"),
#' sightings-survey (S2; output = "s2") and bias-adjusted sightings-only
#' (output = "sobak").
#'
#' If interested in using the `s2bak.bak` model on its own with provided model
#' predictions, see \link[s2bak]{predict.s2bak.bak}.
#'
#' The default method is `s2bak`, which will return a single set
#' of predictions for each species, prioritizing
#' S2 models over bias-adjusted sightings-only models and .
#' If the full S2BaK model is provided, the function will
#' make predictions using \link[s2bak]{s2bak.predict.s2} for species with
#' sightings and survey data, and adjusted sightings-only predictions using BaK
#' for species without survey data.
#' If only `s2bak.so` and `s2bak.bak` models are provided
#' (that is, an S2BaK class model with S2 is NA), the functions returns only
#' adjusted sightings-only predictions, requiring trait data for the species.
#' Conversely, if only `s2bak.s2` model is provided, then predictions are made
#' only with s2bak.predict.s2.
#'
#' #' If the `all` is selected, all sets of predictions are returned as a
#' list: `so`, `s2`, `sobak` and `s2bak`.
#' If models are missing, then only available predictions are made.
#' For example, if the `s2bak` model contains only
#' `s2bak.s2` and `s2bak.so` models, then only those model
#' predictions will be returns without BaK adjustment.
#'
#' @param useReadout Whether to force useReadout in SOS2
#' @param ncores Paralellize SO/S2 predictions, if NA then auto-pick it
#' @param ... Any additional arguments for the predict.fun
#' @return Model predictions as a data.frame with columns for each species and
#' rows for each location. If output = "all" is selected, a lot of
#' data.frame predictions are returned.
#' @rdname predict.s2bak
#' @export predict.s2bak
#' @export
predict.s2bak <- function(model,
                          newdata,
                          trait = NA,
                          predict.fun,
                          predict.bak.fun,
                          output = c("s2bak", "all", "so", "s2", "sobak")[1],
                          ncores = 1,
                          useReadout = FALSE,
                          ...) {

  if (!(output %in% c("so", "s2", "all", "sobak", "s2bak"))){
    stop(paste("Invalid output selected. Specify output as",
          "`s2bak`, `all`, `so`, `s2` or `sobak`"))
  }

  # Check if empty = TRUE
  check_mods <- unlist(lapply(model, FUN = function(x){
                                                        return(x@empty)
                                                      }))
  names(check_mods) <- names(model)

  if (all(check_mods)) {
    stop("No valid models provided.")
  }

  # We output a list if we use `all`
  predictions <- list()

  # We need to make SO predictions in every case except s2 output
  if (output %in% c("so", "sobak", "all", "s2bak")) {
    # Do we have a model?
    # `check_mods` will be TRUE if NA
    if (check_mods["s2bak.so"]) {
      if (output == "all") {
        warning(paste("s2bak.so sub-model missing:",
                      "sightings-only predictions excluded from `all`."))

      } else {
        stop("s2bak.so sub-model missing from `s2bak` model.")

      }

    } else {
      predictions[["so"]] <- predict.s2bak.so(model@s2bak.so,
                                              newdata, predict.fun,
                                              useReadout = useReadout,
                                              ncores = ncores, ...)

    }
  }

  if (output %in% c("s2", "all", "s2bak")) {
    # Do we have a model?
    # `check_mods` will be TRUE if NA
    if (check_mods["s2bak.s2"]) {
      if (output == "all") {
        warning(paste("s2bak.s2 sub-model missing:",
                      "survey-sightings predictions excluded from `all`."))

      } else {
        stop("s2bak.s2 sub-model missing from `s2bak` model.")

      }

    } else {
      predictions[["s2"]] <- predict.s2bak.s2(model@s2bak.s2,
                                              newdata, predict.fun,
                                              useReadout = useReadout,
                                              ncores = ncores, ...)

    }
  }

  if (output %in% c("sobak", "all", "s2bak")) {
    # We need s2bak.SO and s2bak.BaK
    if (check_mods["s2bak.so"] | check_mods["s2bak.bak"]) {
      if (output == "all") {
        warning(paste("s2bak.so, s2bak.bak sub-model missing:",
                      "bias-adjusted sightings-only predictions",
                      "excluded from `all`."))

      } else {
        stop("s2bak.so, s2bak.bak sub-model missing from `s2bak` model.")

      }

    } else {
      predictions[["sobak"]] <- predict.s2bak.bak(model@s2bak.bak,
                                              predictions = predictions[["so"]],
                                              trait, newdata,
                                              predict.bak.fun = predict.bak.fun)

    }

  }

  if (output %in% c("s2bak", "all")) {
    # Can't
    if (any(check_mods)) {
      if (output == "all") {
        warning(paste("Missing one or more sub-models: `s2bak` predictions",
                      "excluded from `all`."))

      } else {
        stop("Missing sub-models: cannot make output = `s2bak` predictions.")

      }

    } else {
      # get full species list
      speciesList <- colnames(predictions[["sobak"]])
      # get S2 species
      spl_s2 <- colnames(predictions[["s2"]])
      # get SO species
      spl_so <- speciesList[!(speciesList %in% spl_s2)]

      predictions[["s2bak"]] <- as.data.frame(
                                  cbind(predictions[["s2"]][, spl_s2],
                                        predictions[["sobak"]][, spl_so])
                                )

      colnames(predictions[["s2bak"]]) <- c(spl_s2, spl_so)

      # Re-order columns
      predictions[["s2bak"]] <- predictions[["s2bak"]][, speciesList]
    }
  }

  # if output is all, return the full list
  # otherwise return only what the user specified
  if (output == "all") {
    return(predictions)
  } else {
    return(predictions[[output]])
  }

}

#' @description The function predict.s2bak.s2 is used for `s2bak.so` and
#' `s2bak.s2` models. It automatically detects which
#' model class is provided and makes the appropriate adjustments.
#' \link[s2bak]{predict.s2bak.so} is a wrapper function
#' that detects the class of the inputed model and makes the appropriate
#' prediction.
#'
#' For `s2bak.s2` models, the binary variable denoting sightings-only or survey
#' is stored and will be checked by the function.
#'
#' @param model Models of class `s2bak.so` or `s2bak.s2` can
#' be used to make predictions for each species fitted.
#' If the object does not have stored SDMs, it will check to see if there
#' is readout (alternatively, a readout can be forced with useReadout = T).
#'
#' @param newdata A data.frame containing the values . All variables needed for
#' prediction should be included.
#' @param predict.fun Predict function linked to the SDM used. Functions have
#' the structure of model and newdata as the first and second arguments,
#' respectively.
#' @param survey_var Character name for the predictor variable determining
#' a site is sightings-only (1) or survey data (0).
#' The column is automatically within the function, and is used to define with
#' formula.
#' @param useReadout logical; if TRUE will do readout over stored SDMs.
#' If there are no SDMs then it will automatically check for readout
#' @param ncores Number of cores to fit the SDMs, default is 1 core but can be
#' automatically set if ncores=NA. If ncores > number of available cores - 1,
#' set to the latter.
#' @param ... Additional arguments passed into function for
#' prediction (predict.fun).
#' @return Generates a matrix of predictions with rows being indices in the
#' data.frame, and columns representing each species.
#' @rdname predict.s2bak
#' @export predict.s2bak.s2
#' @export
predict.s2bak.s2 <- function(model,
                             newdata,
                             predict.fun,
                             useReadout = FALSE,
                             ncores = 1, ...) {
  cat("Predictions using", class(model), "model\n")

  # Get survey_var
  survey_var <- model@survey_var

  # Set cores
  registerDoParallel()
  if (is.numeric(ncores) | is.na(ncores)) {
    NC <- max(detectCores() - 1, 1)
    options(cores = min(NC, ncores, na.rm = TRUE))
  } else {
    stop("Invalid number of cores.")
  }

  # Add SO = 0, regardless of model class
  # Throw a warning if they have 'so' already in the data
  if (survey_var %in% colnames(newdata)) {
    warning(paste(survey_var, "column present in newdata, with same name as",
    "`survey_var`. All values were converted to 0 for prediction."))
  }
  newdata <- as.data.frame(newdata)
  newdata[, survey_var] <- 0

  # Get whether we are dealing with a readout or model
  if (!is.null(model@sdm) & !useReadout) {
    cat("Reading models from argument\n")
    dir <- FALSE
  } else if (useReadout | !is.null(model@options$readout)) {
    cat("Reading models from directory\n")
    dir <- TRUE

    lf <- model@options$readout

    if (any(!(file.exists(lf)))) stop("Missing models in directory.")
  } else {
    stop("Invalid 'model' type, no readout or SDM available")
  }

  # Get length of species
  speciesList <- model@speciesList

  tmp <- as.data.frame(matrix(
    NA,
    nrow = nrow(newdata),
    ncol = length(speciesList)
  ))
  names(tmp) <- speciesList
  tmp[] <- foreach(i = speciesList) %dopar% {
    if (dir) {
      # Read in appropriate model, whose filename
      # should be found in l$options$readout
      tmp_sdm <- readRDS(lf[which(names(lf) == i)])
      # Verify whether this is actually true
      if (tmp_sdm$species == i) {
        tmp_sdm <- tmp_sdm$sdm
      } else {
        stop("Invalid model, incorrect species-filename association")
      }
    } else {
      tmp_sdm <- model@sdm[[i]]
    }

    # NULL model
    if (is.null(tmp_sdm)) {
      tmp1 <- NA
    } else {
      tmp1 <- predict.fun(tmp_sdm, newdata, ...)
    }

    return(tmp1)
  }

  # Check failures
  failed <- apply(tmp, 2, FUN = function(x) {
    all(is.na(x))
  })
  if (any(failed)) {
    warning(paste(
      "Failed to predict for the following species:",
      paste(speciesList[failed], collapse = ", ")
    ))
  }

  return(tmp)
}

#' @rdname predict.s2bak
#' @export predict.s2bak.so
#' @export
predict.s2bak.so <- function(model,
                             newdata,
                             predict.fun,
                             useReadout = FALSE,
                             ncores = 1, ...) {
  return(predict.s2bak.s2(
    model = model,
    newdata = newdata,
    predict.fun = predict.fun,
    useReadout = useReadout,
    ncores = ncores, ...
  ))
}

#' @description predict.s2bak.bak adjusts sightings-only predictions using a
#' fitted s2bak.bak model, using trait, environmental data and
#' sightings-only (\link[s2bak]{predict.s2bak.so}) predictions.
#'
#' @param model Output from fit.s2bak.bak(), with fitted BaK model
#' @param predictions Sightings-only predictions as a matrix or data.frame with
#' rows as sites and columns as species. Assumes as type="response", and rows of
#' data.frame correspond to newdata rows.
#' @param data Environmental data, with rows corresponding to rows of
#' predictions
#' @param trait Trait data, with column 'species' matching those in predictions.
#' @param predict.bak.fun Model function for predicting bias adjustment model
#' (e.g., \link[stats]{predict.glm}). Needs to match `bak.fun`
#' @param truncate Numeric minimum and maximum range of predicted values. Values
#' very close to zero or one cannot be meaningfully distinguished, however
#' these extreme values may have disproportionally large consequences on
#' likelihoods due to logit transformation.
#' @return Model predictions but with adjustments made by the BaK model.
#' Note the default right now is type="response"
#' @rdname predict.s2bak
#' @export predict.s2bak.bak
#' @export
predict.s2bak.bak <- function(model, predictions, trait, data,
                              predict.bak.fun,
                              truncate = c(0.0001, 0.9999)) {
  cat("Predictions using s2bak.bak model\n")
  predictions <- as.matrix(predictions)
  rownames(predictions) <- 1:nrow(predictions)
  predictions2 <- melt(predictions)
  colnames(predictions2) <- c("loc", "species", "pred")

  # scale_so_sp
  trait$pred <- predict.bak.fun(model@bak$bias_sp, trait)
  # scale_so_l
  data$pred <- predict.bak.fun(model@bak$bias_loc, data)

  predictions2$scale_so_l <- data$pred[predictions2$loc]
  predictions2$scale_so_sp <- trait$pred[match(predictions2$species,
                                               trait$species)]
  predictions2$zso_only <- s2bak.truncate(
    as.vector(as.matrix(-log((1 - predictions2$pred) / predictions2$pred))),
    logit(truncate[1]),
    logit(truncate[2])
  )

  predictions2$pred_out <- predict(model@bak$bias_adj,
                                   predictions2,
                                   type = "response")

  # Convert back to wide format
  predictions2 <- dcast(predictions2, loc ~ species, value.var = "pred_out")
  predictions2 <- predictions2[, -which(colnames(predictions2) == "loc")]

  return(predictions2)

}