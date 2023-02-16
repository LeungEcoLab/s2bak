#' @title Adjust SO model predictions using BaK
#'
#' @description Make model adjustments using the output from the BaK output,
#' requiring trait, environmental data and SO predictions.
#'
#' @param predictions Sightings-only predictions as a matrix or data.frame with
#' rows as sites and columns as species. Assumes as type="response", and rows of
#' data.frame correspond to newdata rows.
#' @param bak Output from s2bak.BaK(), with fitted BaK model
#' @param data Environmental data, with rows corresponding to rows of
#' predictions
#' @param trait Trait data, with column 'species' matching those in predictions.
#' @return Model predictions but with adjustments made by the BaK model.
#' Note the default right now is type="response"
#' @rdname s2bak.predict
#' @export
s2bak.predict.BaK <- function(predictions, bak, trait, data) {
  predictions <- as.matrix(predictions)
  rownames(predictions) <- 1:nrow(predictions)
  predictions2 <- melt(predictions)
  colnames(predictions2) <- c("loc", "species", "pred")

  # scale_so_sp
  trait$pred <- predict.glm(bak$bak$bias_sp, trait)
  # scale_so_l
  data$pred <- predict.glm(bak$bak$bias_loc, data)

  predictions2$scale_so_l <- data$pred[predictions2$loc]
  predictions2$scale_so_sp <- trait$pred[match(predictions2$species,
                                                trait$species)]
  predictions2$zso_only <- s2bak.truncate(
        as.vector(as.matrix(-log((1 - predictions2$pred) / predictions2$pred))),
        -15,
        15
      )

  predictions2$pred_out <- predict(bak$bak$bias_adj,
                                    predictions2,
                                    type = "response")

  # Convert back to wide format
  predictions2 <- dcast(predictions2, loc ~ species, value.var = "pred_out")
  predictions2 <- predictions2[, -which(colnames(predictions2) == "loc")]

  return(predictions2)
}

#' @title Make predictions using fitted SO, S2, BaK or S2BaK class models
#'
#' @description The function automatically detects which model class is used,
#' which can be either the output from s2bak.S2 or s2bak.SO.
#'
#' @param model Fitted SO or S2 models of class `s2bak.SO` or `s2bak.S2` to
#' use for prediction. If the object does
#' not have stored SDMs, it will check to see if there is readout
#' (alternatively, you could force readout with useReadout = T).
#' @param newdata A data.frame containing the values . All variables needed for
#' prediction should be included.
#' @param predict.fun Predict function linked to the SDM used. The default used
#' is predict.gam from the mgcv package. Functions have the structure of model
#' and newdata as the first and second arguments, respectively.
#' @param class Whether we're dealing with SO or S2 models, which is only used
#' when providing a directory so we can identify which class
#' we are interested in.
#' @param useReadout logical; if TRUE will do readout over stored SDMs.
#' If there are no SDMs then it will automatically check for readout
#' @param ncores Number of cores to fit the SDMs, default is 1 core but can be
#' automatically set if ncores=NA. If ncores > number of available cores - 1,
#' set to the latter.
#' @param ... Additional arguments passed into function for
#' prediction (predict.fun).
#' @return Generates a matrix of predictions with rows being indices in the
#' data.frame, and columns representing each species.
#' @rdname s2bak.predict
#' @export
s2bak.predict.SOS2 <- function(model,
                                newdata,
                                predict.fun = predict.gam,
                                useReadout = FALSE,
                                ncores = 1, ...) {
  cat("Predicting of class", class(model), "\n")

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
  if ("so" %in% colnames(newdata)) {
    warning("'so' column present in newdata. All values were converted to 0.")
  }
  newdata <- as.data.frame(newdata)
  newdata$so <- 0

  # Get whether we are dealing with a readout or model
  if (!is.null(model$sdm) & !useReadout) {
    cat("Reading models from argument\n")
    dir <- FALSE
  } else if (useReadout | !is.null(model$options$readout)) {
    cat("Reading models from directory\n")
    dir <- TRUE

    lf <- model$options$readout

    if (any(!(file.exists(lf)))) stop("Missing models in directory.")
  } else {
    stop("Invalid 'model' type, no readout or SDM available")
  }

  # Get length of species
  speciesList <- model$speciesList

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
      tmp_sdm <- model$sdm[[i]]
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

#' @title Predict with S2BaK objects
#'
#' @description Wrapper function that detects the class of the inputed model
#' and makes the appropriate prediction.
#'
#' If the provided model is s2bak.S2 or s2bak.SO, predictions will be made
#' using s2bak.predict.SOS2. If an SO model with BaK is provided (that is, an
#' S2BaK class model with S2 = NULL, for example through combine without a
#' provided S2), the functions returns adjusted predictions requiring trait
#' data for the species. If the full S2BaK model (SO, S2 and BaK are provided),
#' the function will make predictions with S2 for species with sightings and
#' survey data, and adjusted predictions using SO and BaK for species
#' without survey data.
#'
#' @param model Outputted model from s2bak.SO, s2bak.S2, s2bak.S2BaK or
#' s2bak.combine.
#' @param newdata Environmental data to predict
#' @param trait Trait data, optional if model is S2 or SO
#' @param predict.fun prediction function for SDM (for S2 and SO)
#' @param useReadout Whether to force useReadout in SOS2
#' @param ncores Paralellize SO/S2 predictions, if NA then auto-pick it
#' @param ... Any additional arguments for the predict.fun
#' @return Model predictions as a data.frame with columns for each species and
#' rows for each location
#' @rdname s2bak.predict
#' @export
s2bak.predict <- function(model,
                          newdata,
                          trait = NA,
                          predict.fun = predict.gam,
                          ncores = 1,
                          useReadout = FALSE,
                          ...) {
  # Simplest cases, which require a call to s2bak.predict.SOS2

  if (class(model) == "s2bak.SO" | class(model) == "s2bak.S2") {
    return(s2bak.predict.SOS2(model, newdata, predict.fun,
                              useReadout = useReadout,
                              ncores = ncores, ...))
  } else if (class(model) == "s2bak.S2BaK") {
    if (is.null(model$s2bak.SO) | is.null(model$s2bak.BaK)) {
      stop("Missing SO or BaK model(s).")
    }
    if (all(is.na(trait))) stop("Missing trait data for BaK prediction.")
    # Get species list, and identify which ones are S2 and with ones are SO-BaK
    # SO should contain all species
    speciesList <- model$s2bak.SO$speciesList

    # Check if NULL
    if (!is.null(model$s2bak.S2)) {
      speciesList.s2 <- model$s2bak.S2$speciesList
      speciesList.so <- speciesList[!(speciesList %in% speciesList.s2)]

      # Make predictions for S2 species
      out.S2 <- s2bak.predict.SOS2(model = model$s2bak.S2, newdata = newdata,
                                    predict.fun = predict.fun,
                                    useReadout = useReadout, ncores = ncores, ...)
    } else {
      # No S2, fit for all species
      speciesList.so <- speciesList
    }

    # Make predictions for SO species
    out.SO <- s2bak.predict.SOS2(model = model$s2bak.SO, newdata = newdata,
                                  predict.fun = predict.fun,
                                  useReadout = useReadout, ncores = ncores, ...)

    # Make adjustment for SO species
    out.SO <- s2bak.predict.BaK(out.SO, model$s2bak.BaK, trait, newdata)

    # Combine and output
    out <- cbind(out.SO, out.S2)
    out <- out[, speciesList]
    return(out)
  } else {
    stop("Invalid class: requires s2bak.SO, s2bak.S2 or s2bak.S2BaK.")
  }
}
