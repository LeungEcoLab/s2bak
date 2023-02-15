#' @name s2bak
#' @title Build sightings-only or S2 species distribution models for
#' multiple species.
#'
#' @description s2bak.SO function fits SDMs for each provided species within
#' the same system, using a specified SDM approach (or the default which are
#' GAMs from the mgcv package). Parallelisation is possible when processing each
#' SDM, with the default being 1 core.
#'
#' The s2bak.S2 function fits SDMs using species sightings, background sites and
#' survey sites, differentiating between them using a binary 'so' predictor,
#' denoting sightings-only (1) or survey (0).
#'
#' Saving SDMs to the output may be computationally intensive, particularly with
#' large datasets and many species. To reduce the possibility of crashing,
#' readout and the version = "short" may be used, which does not output the
#' fitted models but instead saves it to the directory specified in readout.
#'
#' @param formula Formula for the SDMs, which serves as input for the given SDM.
#' Assumes the structure follows "Y ~ X". Alternatively, a list of formulas can
#' be provided with index names corresponding to species. In this case, species
#' will be fit using their corresponding formula. The response variable can
#' have any name, as the function will detect it and name the columns
#' accordingly.
#' @param data Full environmental data used for fitting.
#' @param obs Species observations as a data.frame, with a column for species
#' name (must be labelled 'species') and column of index of observations to
#' reflect presences. If the index column name is not found in 'data', it
#' assumes row number.
#' @param surv Survey data, given as a data.frame of row/index as the first
#' column for the data following the same name as obs, and each column
#' representing a given species presence (1) or absence (0). It will add the an
#' additional binary predictor to the formula(s), "so", which denotes
#' sightings-only or not. If 'so' is already in formula (e.g. if modifying
#' the variable in any way, then set surv.formula = FALSE). If left as NA,
#' it will fit the SDMs as presence-only models with the function of choice.
#' NOTE: Assumes all species are provided in the survey data, and that they are
#' surveyed over the same sites (i.e. matrix-type).
#' @param background Background sites (pseudo-absences) used to fit the
#' presence-only model, provided as a vector of indices of data (following the
#' same column name as observations). If the index column name is not found in
#' 'data', it assumes row number within 'data'. If left as NA, it will randomly
#' sample 'nbackground' sites, with or without overlap ('overlapbackground').
#' Currently, only one set of background sites can be used.
#' @param sdm.fun Model (as function) used for fitting (the default is
#' \link[mgcv]{gam} from the mgcv package). The function must have the formula
#' as their first argument, and 'data' as the parameter for the dataset
#' (including presences and background sites within the data.frame).
#' A wrapper function for \link[maxnet]{maxnet} is also available as
#' \link[s2bak]{s2bak.maxnet}.
#' @param overlapbackground Whether sampled background sites that overlap with
#' observations should be included. By default, it allows overlap. If FALSE,
#' number of background sites may be less than specified or provided.
#' @param nbackground Number of background sites to sample. Only applies
#' if background = FALSE.
#' @param surv.formula Whether the binary variable 'so', denoting whether sites
#' are sightings-only or survey, should be added to the formula(s). If survey
#' data is not provided, then 'so' will not be added to the formula(s)
#' regardless of whether it is TRUE or FALSE. If there is survey data and
#' surv.formula = FALSE, then 'so' will not be added, and must be included in
#' the initial function call (or it will throw an error).
#' @param ncores Number of cores to fit the SDMs, default is 1 core but can be
#' automatically set if ncores=NA. If ncores > number of available cores - 1,
#' set to the latter.
#' @param readout Directory to save fitted SDMs and background sites. If NA, it
#' will not save any SDMs. Provides an additional output that shows where the
#' SDM is saved (with file name). The output in this directory can later be
#' used in other functions such as \link[s2bak]{s2bak.predict.SOS2}.
#' @param version Whether the SDMs should be included in the output. With
#' "short", no the fitted SDMs are not provided. Setting to "full" (default)
#' will output the list with all SDMs. Setting to "short" and combined with
#' readout, can considerably reduce RAM usage while saving the progress so far,
#' which is useful when dealing with many species or large datasets.
#' @param ... Other arguments that are passed to the SDM function (sdm.fun).
#' @return An object of class "s2bak.S2", providing fitted SDMs for each species
#' based on the provided SDM modelling approach. The primary difference between
#' SO and S2 models are the additional data points from the survey data, and an
#' additional binary predictor 'so' which denotes whether the data is from
#' presence-background (1) or presence-absence data (0).
#' @rdname s2bak
#' @export
s2bak.S2 <- function(formula, data, obs, surv = NA, background = NA,
                     sdm.fun = gam, overlapbackground = TRUE, 
                     nbackground = 10000,
                     surv.formula = TRUE,
                     ncores = 1,
                     readout = NA, version = c("full", "short")[1], ...) {
  # Set cores
  registerDoParallel()
  if (is.numeric(ncores) | is.na(ncores)) {
    NC <- max(detectCores() - 1, 1)
    options(cores = min(NC, ncores, na.rm = TRUE))
  } else {
    stop("Invalid number of cores.")
  }

  # Get the list of species in the data
  specieslist <- as.character(unique(obs$species))
  cat("Fitting", ifelse(all(is.na(surv)), "sightings-only", "S2"),
      "model for", length(specieslist), "species\n")
  cat(ncores, "core(s)\n")

  # Add / if not there
  if (!is.na(readout)) {
    readout <- paste0(readout, ifelse(substr(readout, nchar(readout),
                      nchar(readout)) == "/", "", "/"))
  }

  # Index name, if it's row #s or a specific column
  if (ncol(obs) != 2 | !("species" %in% colnames(obs))) {
    stop("Invalid columns in 'obs'.")
  } else {
    # Get column name of obs that isn't species
    ind <- colnames(obs)[which(colnames(obs) != "species")]
    if (!(ind %in% colnames(data))) {
      ind <- NA # Row number instead of column name
    }
  }

  # Check formula type
  if (typeof(formula) == "language") {
    # Only one formula
    flong <- FALSE
  } else if (typeof(formula) == "list") {
    # Multiple formulas, but we still need to check if all species are there
    flong <- TRUE
    if (!all(specieslist %in% names(formula))) {
      stop(cat("Incorrect formula type.",
      "Please provide a single formula or a list of formulas with each",
      "index named after the corresponding species."))
    }
  } else {
    stop(cat("Invalid formula type. Please provide a single formula or a",
    "list of formulas with each index named after the corresponding species."))
  }

  # Output for the model
  # Saves the species list and functions used for fitting/predicting
  # Also save the arguments
  l <- list(
    species.list = specieslist,
    options = list(
      call = formula,
      functions = list(
        sdm = sdm.fun
      ),
      version = version,
      readout = readout,
      overlapbackground = overlapbackground
    )
  )

  class(l) <- ifelse(all(is.na(surv)), "s2bak.SO", "s2bak.S2")

  # Add background options to the output
  if (is.na(background)) {
    l$options$background <- "Randomly sampled background sites"
    l$options$nbackground <- nbackground
    # Sample background sites, if NA
    # Also provide background indices
    if (is.na(background)) {
      # If there is a specific index name, use that, otherwise use rownumber
      if (is.na(ind)) {
        background <- sample(1:nrow(data), nbackground)
      } else {
        background <- sample(data[, ind], nbackground)
      }
    }
  } else {
    l$options$background <- "User-provided background sites"
    l$options$nbackground <- length(background)
  }
  l$options$sitesbackground <- background

  # Save our options, and also add a filepath
  if (!is.na(readout)) {
    fopts <- paste0(readout,
                    ifelse(all(is.na(surv)),
                        "opts.SO.rds",
                        "opts.S2.rds"))
    saveRDS(l, fopts)
    cat("Output saved to ", fopts, "\n")
  }

  # get response variable name
  yy <- as.character(ifelse(flong, formula[[1]][[2]], formula[[2]]))

  # Add 'so' if we have survey data and surv.formula = TRUE
  if (!all(is.na(surv)) & surv.formula) {
    formula <- s2bak.addPred(formula)
  }

  # Fit SDM for each species
  # Saves output of each one as well
  cat("Fitting SDMs for", length(specieslist), "species\n")
  l$sdm <- foreach(i = specieslist) %dopar% {
    # Generate fitting data using 'data', 'obs' and 'background'
    # If overlapbackground == FALSE, remove overlaps
    if (!overlapbackground) {
      background <- background[!(background %in% obs[obs$species == i, 
                                  which(colnames(obs) != "species")])]
    }
    # Get the list of indices
    ind2 <- c(obs[obs$species == i, which(colnames(obs) != "species")],
              background)
    # Row number
    if (is.na(ind)) {
      tmp_dat <- data[ind2, ]
    } else {
      # Obs index name
      tmp_dat <- data[match(ind2, data[, ind]), ]
    }

    # Add `yy` as column
    tmp_dat <- cbind(tmp_dat, 0)
    colnames(tmp_dat)[ncol(tmp_dat)] <- yy
    # Set first N rows as sightings
    tmp_dat[1:nrow(obs[obs$species == i, ]), yy] <- 1
    tmp_dat <- as.data.frame(tmp_dat)

    # If we have a list of formulas, find corresponding index
    if (flong) {
      ff <- formula[[i]]
    } else {
      ff <- formula
    }

    # Add survey data, if it exists
    if (!all(is.na(surv))) {
      # Generate survey data
      if (is.na(ind)) {
        # Note: assumes that the first column of surv is the index
        #### MAYBE CHANGE THIS ####
        tmp_surv <- data[surv[, 1], ]
      } else {
        tmp_surv <- data[match(surv[, ind], data[, ind]), ]
      }
      tmp_surv <- as.data.frame(tmp_surv)
      tmp_surv$pa <- surv[, i]

      # Add 'so': sightings only for `tmp_dat` (1) or
      # survey data for `tmp_surv` (0)
      tmp_dat$so <- 1
      tmp_surv$so <- 0

      tmp_dat <- rbind(tmp_dat, tmp_surv)

      # Check if 'so' is in our formula.. if it isn't throw an error
      if (!("so" %in% labels(terms(ff)))) {
        stop("Provided survey data but 'so' is missing as a predictor")
      }
    }

    tmp_sdm <- tryCatch(sdm.fun(ff, data = tmp_dat, ...),
      error = function(e) {
        return(NULL)
      }
    )

    if (!is.na(readout)) {
      fname <- paste0(readout, gsub("\\s+", "_", i), "-", class(l), ".rds")
      saveRDS(list(species = i, sdm = tmp_sdm), fname)
      cat("\t", fname, "\n")

      # Also save l in case of failure
      saveRDS(l, fopts)
    }

    if (version == "full") {
      return(tmp_sdm)
    } else if (version == "short") {
      return(NA) # This will be NULL afterwards
    } else {
      stop("Invalid version. Please specify either 'long' or 'short'")
    }
  }
  names(l$sdm) <- specieslist

  wh <- which(unlist(lapply(l$sdm, is.null)))
  if (length(wh) > 0) {
    l$failure <- specieslist[wh]
  } else {
    l$failure <- NULL
  }


  if (version == "short") {
    l$sdm <- NULL
  } else if (!(version %in% c("full", "short"))) {
    stop("Invalid version. Please specify either 'long' or 'short'")
  }

  # Save options one last time, logging the filenames
  if (!is.na(readout)) {
    lf <- paste0(readout, gsub("\\s+", "_", specieslist), "-", class(l), ".rds")
    names(lf) <- specieslist
    # Don't include options in the output list of files
    l$options$readout <- lf
    # Save opts one more time
    saveRDS(l, fopts)
  }

  if (length(l$failure) > 0) {
    warning(paste0(
      "Models failed to fit for the following species: ",
      paste0(l$failure, collapse = ", ")
    ))
  }

  return(l)
}

#' @rdname s2bak
#' @export
s2bak.SO <- function(formula, data, obs, background = NA,
                     sdm.fun = gam, overlapbackground = TRUE,
                     nbackground = 10000,
                     ncores = 1,
                     readout = NA, version = c("full", "short")[1], ...) {
  return(s2bak.S2(formula, data, obs,
    surv = NA, background = background,
    sdm.fun = sdm.fun,
    overlapbackground = overlapbackground, nbackground = nbackground,
    ncores = ncores, surv.formula = TRUE,
    readout = readout, version = version, ...
  ))
}

#' @title Generate bias adjustment kernel (BaK) for a fitted sightings-only SDM.
#'
#' @description Builds bias adjustment kernel (BaK) for a fitted sightings-only
#' SDM. Provides three models, as generalized linear models (GLMs): modelling
#' location bias, modelling species bias and an adjustment model that combines
#' model predictions with the output from the other two models.
#'
#' The dataset assumes the rows for the predictions, survey and data arguments
#' all align and match.
#'
#' @param predictions Sightings-only (SO) model predictions over the survey
#' sites for all species, beyond those found in the survey data, as a matrix
#' with columns for each species and rows for each site.
#' @param surv Survey data as a matrix, with rows corresponding to each site.
#' @param data Survey environmental data with columns for environmental data
#' and rows for each site. Assumes all columns correspond to relevant
#' environmental data and match with predictions/survey data. Also assumes
#' that inputted columns are all relevant in modelling location bias.
#' @param trait Full trait data for the species predictions, as a data.frame
#' with 'species' as a column and relevant traits for the remainder.
#' Like with the predictions, the species in the dataset do not necessarily
#' have to possess survey data, but will be used in the final adjustment model
#' as final output.
#' @return Bias adjustment models, the kernels (location and species),
#' as a second-order GLM.
#' @rdname s2bak
#' @export
s2bak.BaK <- function(predictions, surv, data, trait) {
  wh <- which(!colnames(surv) %in% trait$species)
  if (length(wh) > 0) {
    warning(paste("Columns from survey data missing from trait data and",
                  "not excluded from the fitting:", colnames(surv)[wh]))
    surv <- surv[, -wh]
  }
  specieslist <- colnames(surv)

  # Model outputs
  out <- list(
    species = specieslist,
    fullspecies = unique(trait$species)
  )
  class(out) <- "s2bak.BaK"

  cat("Fitting bias adjustment kernel on", length(specieslist),
  "survey species and", length(out$fullspecies), "species total.\n")

  # Throw error if predictions, surv and data are not aligned
  if (length(unique(c(nrow(predictions), nrow(data), nrow(surv)))) > 1) {
    stop(cat("Differing rows for predictions, survey and environment data:",
            "sites/rows must correspond with each other."))
  }

  # Throw error if species missing from traits data.frame
  if (!all(specieslist %in% trait$species) | !all(specieslist %in% colnames(predictions))) {
    stop("Species missing from trait data.")
  }

  # Generate polynomial predictor variable names
  # Can replace with more complex models, if desired
  # (but right now it's GLM with this)
  tn <- colnames(trait)[colnames(trait) != "species"]
  en <- colnames(data)
  numer_env <- unlist(lapply(data[1, en], is.numeric))
  numer_tr <- unlist(lapply(trait[1, tn], is.numeric))
  names_xbe <- c(en[numer_env], paste("I(", en[numer_env], "^2)", sep = ""))
  names_xbt <- c(tn[numer_tr], paste("I(", tn[numer_tr], "^2)", sep = ""))

  # Get mean for replacing missing values (for the full dataset)
  msd <- rep(NA, time = length(tn))
  names(msd) <- tn
  for (vv in tn) {
    if (numer_tr[vv]) {
      msd[vv] <- mean(trait[, vv], na.rm = TRUE)
    } else {
      msd[vv] <- 0
    }
  }

  # Location biases, based on environment
  fit_l <- data
  # Get summed predictions
  fit_l$so_only <- apply(predictions, 1, sum, na.rm = TRUE)
  # Get summed occurrences
  fit_l$pa <- apply(surv, 1, sum, na.rm = TRUE)
  # Get log-ratio
  fit_l$lr <- log((fit_l$pa + 1) / (fit_l$so_only + 1))

  # Species biases, based on traits for survey species
  trait2 <- trait[match(specieslist, trait$species), ]
  fit_sp <- data.frame(species = specieslist, stringsAsFactors = FALSE)
  fit_sp <- cbind(fit_sp, trait2[, -which(colnames(trait2) == "species")])
  for (i in 1:nrow(fit_sp)) {
    ## Does this work? Should we just assume that this is provided
    #### TBD ####
    wh <- which(is.na(fit_sp[i, ]))
    if (length(wh) > 0) {
      warning("Missing trait data, which were imputated using mean of the column.")
      fit_sp[i, wh] <- msd[colnames(wh)]
    }
  }
  # Get summed predictions
  fit_sp$so_only <- apply(predictions[, specieslist], 2, sum, na.rm = TRUE)
  # Get summed occurrences
  fit_sp$pa <- apply(surv[, specieslist], 2, sum, na.rm = TRUE)
  # Get log-ratio
  fit_sp$lr <- log((fit_sp$pa + 1) / (fit_sp$so_only + 1))

  out$bak <- list()
  # Fit bias kernels
  out$bak$bias_loc <- glm(formula(paste("lr ~", paste(names_xbe, collapse = "+"))), data = fit_l)
  out$bak$bias_sp <- glm(formula(paste("lr ~", paste(names_xbt, collapse = "+"))), data = fit_sp)

  # Generate the final data.frame to model bias adjustment (based on mk_bak)
  # Include all species, beyond just the ones found in the survey
  trait$pred <- predict(out$bak$bias_sp, trait)
  fit_adj <- lapply(specieslist, FUN = function(x) {
    ddf <- data.frame(
      loc = 1:nrow(data), species = x, pa = 0,
      zso_only = NA, scale_so_l = NA, scale_so_sp = NA,
      stringsAsFactors = FALSE
    )
    if (x %in% colnames(surv)) ddf$pa <- surv[, x]
    ddf$zso_only <- s2bak.truncate(as.vector(as.matrix(-log((1 - predictions[, x]) / predictions[, x]))), -15, 15)
    ddf$scale_so_l <- predict(out$bak$bias_loc, data)
    ddf$scale_so_sp <- trait$pred[trait$species == x]
    return(ddf)
  })
  fit_adj <- do.call(rbind, fit_adj)

  # Fit bias adjustment model
  out$bak$bias_adj <- glm(pa ~ zso_only + scale_so_sp + scale_so_l, family = "binomial", data = fit_adj)

  return(out)
}

#'  @title Fit entire S2Bak model
#'
#' @description Build S2BaK from top to bottom. Has functionality for
#' parallelization, but the default is 1 core.
#'
#' Fits SO models for all species, S2 models for species with survey data and
#' a BaK model for adjusted predictions.
#'
#' Assumes that all columns/variables in 'data' are relevant for the
#' location bias model.
#'
#' @return An S2BaK class object containing S2, SO and BaK.
#' @rdname s2bak
#' @export
s2bak.S2BaK <- function(formula, data, obs, surv, trait,
                        background = NA, sdm.fun = gam, predict.fun = predict.gam, overlapbackground = TRUE, nbackground = 10000,
                        surv.formula = TRUE,
                        ncores = 1,
                        readout = NA, version = c("full", "short")[1], ...) {
  # Output for s2bak.S2BaK
  out <- list()
  class(out) <- "s2bak.S2BaK"

  # Index name, if it's row #s or a specific column
  if (ncol(obs) != 2 | !("species" %in% colnames(obs))) {
    stop("Invalid columns in 'obs'.")
  } else {
    # Get column name of obs that isn't species
    ind <- colnames(obs)[which(colnames(obs) != "species")]
    if (!(ind %in% colnames(data))) {
      ind <- NA # Row number instead of column name
    }
  }

  ## First, fit SO model
  out$s2bak.SO <- s2bak.SO(formula = formula, data = data, obs = obs,
                            background = background, sdm.fun = sdm.fun,
                            overlapbackground = overlapbackground,
                            nbackground = nbackground, ncores = ncores,
                            readout = readout, version = version, ...)

  ## Next, fit S2 model
  ## Only for species with survey data
  if (all(is.na(surv))) stop("Missing survey data.")

  out$s2bak.S2 <- s2bak.S2(formula = formula, data = data,
                            obs = obs[obs$species %in% colnames(surv), ],
                            surv = surv,
                            background = background,
                            sdm.fun = sdm.fun,
                            overlapbackground = overlapbackground,
                            nbackground = nbackground,
                            surv.formula = surv.formula,
                            ncores = ncores, readout = readout,
                            version = version, ...)

  ## Finally, fit BaK
  if (is.na(ind)) {
    #### ASSUMES 1ST COLUMN OF SURV IS ROWNUM ####
    ## Should maybe just detect it
    surv_dat <- data[surv[, 1], ]
    surv2 <- surv[, -1]
  } else {
    surv_dat <- data[match(surv[, ind], data[, ind]),
                            -which(colnames(data) == ind)]
    surv2 <- surv[, -which(colnames(surv) == ind)]
  }

  # Make predictions using out$SO (always type = "response")
  predictions <- s2bak.predict.SOS2(out$s2bak.SO, surv_dat,
                                      doReadout = !is.na(readout),
                                      ncores = ncores, type = "response")

  out$s2bak.BaK <- s2bak.BaK(predictions, surv2, surv_dat, trait)

  return(out)
}

#' @title Combine SO, S2 and BaK objects
#'
#' @description Combines separately fitted SO, S2 and BaK models into a single
#' S2BaK object. The output will be identical to running the entire process in
#' s2bak.S2BaK function, and can be used in the same situations.
#'
#' Models can be partially provided, for instance only SO and BaK, in which case
#' the use of s2bak predict functions will not apply S2 but instead only run
#' predictions with sightings-only and BaK adjustment.
#'
#' @param SO output from s2bak.SO or s2bak.S2 without survey data.
#' @param S2 output from s2bak.S2 function
#' @param BaK output from s2back.BaK function
#' @return Object of class s2bak.S2Bak, equivalent to having run
#' \link[s2bak]{s2bak.S2BaK} for the entire dataset
#' @export
s2bak.combine <- function(SO = NULL, S2 = NULL, BaK = NULL) {
  out <- list(
    s2bak.SO = SO,
    s2bak.S2 = S2,
    s2bak.BaK = BaK
  )
  class(out) <- "s2bak.S2BaK"
  return(out)
}

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

  predictions2$pred.out <- predict(bak$bak$bias_adj,
                                    predictions2,
                                    type = "response")

  # Convert back to wide format
  predictions2 <- dcast(predictions2, loc ~ species, value.var = "pred.out")
  predictions2 <- predictions2[, -which(colnames(predictions2) == "loc")]

  return(predictions2)
}

#' @title Make predictions using fitted SO, S2, BaK or S2BaK class models
#'
#' @description The function automatically detects which model class is used,
#' which can be either the output from s2bak.S2 or s2bak.SO.
#'
#' @param model Fitted SO or S2 models to use for prediction. If the object does
#' not have stored SDMs, it will check to see if there is readout
#' (alternatively, you could force readout with doReadout = T).
#' @param newdata A data.frame containing the values . All variables needed for
#' prediction should be included.
#' @param predict.fun Predict function linked to the SDM used. The default used
#' is predict.gam from the mgcv package. Functions have the structure of model
#' and newdata as the first and second arguments, respectively.
#' @param class Whether we're dealing with SO or S2 models, which is only used
#' when providing a directory so we can identify which class
#' we are interested in.
#' @param doReadout logical; if TRUE will do readout over stored SDMs.
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
                                doReadout = FALSE,
                                ncores = 1, ...) {
  cat("Predicting of class", class(model), "\n")

  # Set cores
  registerDoParallel()
  if (is.numeric(ncores) | is.na(ncores)) {
    NC <- max(detectCores() - 1, 1)
    options(cores = min(NC, ncores, na.rm = RUET))
  } else {
    stop("Invalid number of cores.")
  }

  # Add SO = 0, regardless of model class
  # Throw a warning if they have 'so' already in the data
  if ("so" %in% colnames(newdata)) {
    warning("'so' column present in newdata. All values were converted to 0.")
  }
  newdata$so <- 0

  # Get whether we are dealing with a readout or model
  if (!is.null(model$sdm) & !doReadout) {
    cat("Reading models from argument\n")
    dir <- FALSE
  } else if (doReadout | !is.null(model$options$readout)) {
    cat("Reading models from directory\n")
    dir <- TRUE

    lf <- model$options$readout

    if (any(!(file.exists(lf)))) stop("Missing models in directory.")
  } else {
    stop("Invalid 'model' type, no readout or SDM available")
  }

  # Get length of species
  specieslist <- model$species.list

  tmp <- as.data.frame(matrix(NA, nrow = nrow(newdata), ncol = length(specieslist)))
  names(tmp) <- specieslist
  tmp[] <- foreach(i = specieslist) %dopar% {
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
      paste(specieslist[failed], collapse = ", ")
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
#' @param doReadout Whether to force doReadout in SOS2
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
                          doReadout = FALSE,
                          ...) {
  # Simplest cases, which require a call to s2bak.predict.SOS2

  if (class(model) == "s2bak.SO" | class(model) == "s2bak.S2") {
    return(s2bak.predict.SOS2(model, newdata, predict.fun,
                              doReadout = doReadout,
                              ncores = ncores, ...))
  } else if (class(model) == "s2bak.S2BaK") {
    if (is.null(model$s2bak.SO) | is.null(model$s2bak.BaK)) {
      stop("Missing SO or BaK model(s).")
    }
    if (all(is.na(trait))) stop("Missing trait data for BaK predictin.")
    # Get species list, and identify which ones are S2 and with ones are SO-BaK
    # SO should contain all species
    specieslist <- model$s2bak.SO$species.list

    # Check if NA
    if (!is.null(model$s2bak.S2)) {
      specieslist.s2 <- model$s2bak.S2$species.list
      specieslist.so <- specieslist[!(specieslist %in% specieslist.s2)]

      # Make predictions for S2 species
      out.S2 <- s2bak.predict.SOS2(model = model$s2bak.S2, newdata = newdata,
                                    predict.fun = predict.fun,
                                    doReadout = doReadout, ncores = ncores, ...)
    } else {
      # No S2, fit for all species
      specieslist.so <- specieslist
    }

    # Make predictions for SO species
    out.SO <- s2bak.predict.SOS2(model = model$s2bak.SO, newdata = newdata,
                                  predict.fun = predict.fun,
                                  doReadout = doReadout, ncores = ncores, ...)

    # Make adjustment for SO species
    out.SO <- s2bak.predict.BaK(out.SO, model$s2bak.BaK, trait, newdata)

    # Combine and output
    out <- cbind(out.SO, out.S2)
    out <- out[, specieslist]
    return(out)
  } else {
    stop("Invalid class: requires s2bak.SO, s2bak.S2 or s2bak.S2BaK.")
  }
}

#' @title Modified maxnet function for s2bak
#'
#' @description Wrapper function for maxnet::maxnet function. The function
#' converts the arguments for maxnet to the necessary structure for fitting an
#' SDM within S2BaK, where the first and second arguments are the formula and
#' dataset, respectively.
#'
#' @param formula Formula for fitting the maxnet model
#' @param data Fitting dataset including response variable
#' @return Output of the maxnet model
#' @export
s2bak.maxnet <- function(formula, data, ...) {
  require(maxnet)
  # Get response variable
  yy <- as.character(formula[[2]])
  p <- data[, which(colnames(data) == yy)]
  data <- data[, which(colnames(data) != yy)]
  # If there's a Y variable, remove it from formula
  formula[[2]] <- NULL
  return(maxnet::maxnet(p, data, f = formula, ...))
}