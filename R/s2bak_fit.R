#' @name s2bak
#' @title Build sightings-only or S2 species distribution models for
#' multiple species.
#'
#' @description fit.so function fits SDMs for each provided species within
#' the same system, using a specified SDM approach (or the default which are
#' GAMs from the mgcv package). Parallelisation is possible when processing each
#' SDM, with the default being 1 core.
#'
#' The fit.s2 function fits SDMs using species sightings, background sites and
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
#' be provided with names corresponding to species. In this case, species
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
#' the variable in any way, then set addSurvey = FALSE). If left as NA,
#' it will fit the SDMs as presence-only models with the function of choice.
#' NOTE: Assumes all species are provided in the survey data, and that they are
#' surveyed over the same sites (i.e. matrix-type).
#' @param background Background sites (pseudo-absences) used to fit the
#' presence-only model, provided as a vector of indices of data (following the
#' same column name as observations). If the index column name is not found in
#' 'data', it assumes row number within 'data'. If left as NA, it will randomly
#' sample 'nbackground' sites, with or without overlap ('overlapBackground').
#' Currently, only one set of background sites can be used.
#' @param sdm.fun Model (as function) used for fitting (the default is
#' \link[mgcv]{gam} from the mgcv package). The function must have the formula
#' as their first argument, and 'data' as the parameter for the dataset
#' (including presences and background sites within the data.frame).
#' @param overlapBackground Whether sampled background sites that overlap with
#' observations should be included. By default, it allows overlap. If FALSE,
#' number of background sites may be less than specified or provided.
#' @param nbackground Number of background sites to sample. Only applies
#' if background = NA.
#' @param addSurvey Whether the binary variable 'so' should be added to the
#' formula(s), denoting whether sites are sightings-only or survey. If survey
#' data is not provided, then 'so' will not be added to the formula(s)
#' regardless of whether it is TRUE or FALSE. If there is survey data and
#' addSurvey = FALSE, then 'so' will not be added, and must be included in
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
fit.s2 <- function(formula, data, obs, surv = NA, background = NA,
                     sdm.fun = gam, overlapBackground = TRUE,
                     nbackground = 10000,
                     addSurvey = TRUE,
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

  # Are we fitting SO or S2?
  mode <- ifelse(all(is.na(surv)), "s2bak.SO", "s2bak.S2")

  # Get the list of species in the data
  speciesListFull <- as.character(unique(obs$species))
  if (mode == "s2bak.SO") {
    speciesList <- speciesListFull
  } else if(mode == "s2bak.S2") {
    speciesList <- colnames(surv)
    ## Assumes that if ind is NA, remove first entry (index)
    ## Otherwise remove whatever has the value of `ind`
    if (is.na(ind)) {
      speciesList <- speciesList[-1]
    }else{
      speciesList <- speciesList[-which(speciesList == ind)]
    }
  }
  cat("Fitting", mode, "model for",
      length(speciesList), "out of",
      length(speciesListFull), "species\n")
  cat(ncores, "core(s)\n")

  # Add / if not there
  if (!is.na(readout)) {
    readout <- paste0(readout, ifelse(substr(readout, nchar(readout),
                      nchar(readout)) == "/", "", "/"))
  }

  # Check formula type
  if (typeof(formula) == "language") {
    # Only one formula
    flong <- FALSE
  } else if (typeof(formula) == "list") {
    # Multiple formulas, but we still need to check if all species are there
    flong <- TRUE
    if (!all(speciesList %in% names(formula))) {
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
    speciesList = speciesList,
    options = list(
      call = formula,
      functions = list(
        sdm = sdm.fun
      ),
      version = version,
      readout = readout,
      overlapBackground = overlapBackground
    )
  )

  class(l) <- mode

  # Add background options to the output
  if (all(is.na(background))) {
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

  # Add 'so' if we have survey data and addSurvey = TRUE
  if (!all(is.na(surv)) & addSurvey) {
    formula <- s2bak.addPred(formula)
  }

  # Fit SDM for each species
  # Saves output of each one as well
  l$sdm <- foreach(i = speciesList) %dopar% {
    # Generate fitting data using 'data', 'obs' and 'background'
    # If overlapBackground == FALSE, remove overlaps
    if (!overlapBackground) {
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
  names(l$sdm) <- speciesList

  wh <- which(unlist(lapply(l$sdm, is.null)))
  if (length(wh) > 0) {
    l$failure <- speciesList[wh]
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
    lf <- paste0(readout, gsub("\\s+", "_", speciesList), "-", class(l), ".rds")
    names(lf) <- speciesList
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
fit.so <- function(formula, data, obs, background = NA,
                     sdm.fun = gam, overlapBackground = TRUE,
                     nbackground = 10000,
                     ncores = 1,
                     readout = NA, version = c("full", "short")[1], ...) {
  return(fit.s2(formula, data, obs,
    surv = NA, background = background,
    sdm.fun = sdm.fun,
    overlapBackground = overlapBackground, nbackground = nbackground,
    ncores = ncores, addSurvey = TRUE,
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
#' Assumes the first column is the index matching the environmental data.
#' @param data Environmental data with columns for environmental data
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
fit.bak <- function(predictions, surv, data, trait) {
  wh <- which(!colnames(surv)[-1] %in% trait$species)
  if (length(wh) > 0) {
    warning(paste("Columns from survey data missing from trait data and",
                  "not excluded from the fitting:", colnames(surv)[wh]))
    surv <- surv[, -wh]
  }
  speciesList <- colnames(surv)[-1]

  # Model outputs
  out <- list(
    speciesList = speciesList,
    speciesListFull = unique(trait$species)
  )
  class(out) <- "s2bak.BaK"

  cat("Fitting bias adjustment kernel on", length(speciesList),
  "survey species and", length(out$speciesListFull), "species total.\n")

  # Throw error if predictions, surv and data are not aligned
  if (length(unique(c(nrow(predictions), nrow(data), nrow(surv)))) > 1) {
    stop(cat("Differing rows for predictions, survey and environment data:",
            "sites/rows must correspond with each other."))
  }

  # Throw error if species missing from traits data.frame
  if (!all(speciesList %in% c(trait$species, colnames(predictions)))) {
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
  fit_l <- as.data.frame(data)
  # Get summed predictions
  fit_l$so_only <- apply(predictions, 1, sum, na.rm = TRUE)
  # Get summed occurrences
  fit_l$pa <- apply(surv, 1, sum, na.rm = TRUE)
  # Get log-ratio
  fit_l$lr <- log((fit_l$pa + 1) / (fit_l$so_only + 1))

  # Species biases, based on traits for survey species
  trait2 <- trait[match(speciesList, trait$species), ]
  fit_sp <- data.frame(species = speciesList, stringsAsFactors = FALSE)
  fit_sp <- cbind(fit_sp, trait2[, -which(colnames(trait2) == "species")])
  for (i in 1:nrow(fit_sp)) {
    ## Does this work? Should we just assume that this is provided
    #### TBD ####
    wh <- which(is.na(fit_sp[i, ]))
    if (length(wh) > 0) {
      warning(cat("Missing trait data, which were imputed",
                  "using mean of the column."))
      fit_sp[i, wh] <- msd[colnames(wh)]
    }
  }

  # Get summed predictions
  fit_sp$so_only <- apply(predictions[, speciesList], 2, sum, na.rm = TRUE)
  # Get summed occurrences
  fit_sp$pa <- apply(surv[, speciesList], 2, sum, na.rm = TRUE)
  # Get log-ratio
  fit_sp$lr <- log((fit_sp$pa + 1) / (fit_sp$so_only + 1))

  out$bak <- list()
  # Fit bias kernels
  out$bak$bias_loc <- glm(formula(
                            paste("lr ~", paste(names_xbe, collapse = "+"))),
                          data = fit_l)
  out$bak$bias_sp <- glm(formula(
                            paste("lr ~", paste(names_xbt, collapse = "+"))),
                          data = fit_sp)

  # Generate the final data.frame to model bias adjustment (based on mk_bak)
  # Include all species, beyond just the ones found in the survey
  trait$pred <- predict(out$bak$bias_sp, trait)
  fit_adj <- lapply(speciesList, FUN = function(x) {
    ddf <- data.frame(
      loc = 1:nrow(data), species = x, pa = 0,
      zso_only = NA, scale_so_l = NA, scale_so_sp = NA,
      stringsAsFactors = FALSE
    )
    if (x %in% colnames(surv)) ddf$pa <- surv[, x]
    ddf$zso_only <- s2bak.truncate(
                            as.vector(
                                as.matrix(-log(
                                    (1 - predictions[, x]) /
                                    predictions[, x]
                                  ))
                            ),
                            -15,
                            15
                          )
    ddf$scale_so_l <- predict(out$bak$bias_loc, data)
    ddf$scale_so_sp <- trait$pred[trait$species == x]
    return(ddf)
  })
  fit_adj <- do.call(rbind, fit_adj)

  # Fit bias adjustment model
  out$bak$bias_adj <- glm(pa ~ zso_only + scale_so_sp + scale_so_l,
                          family = "binomial", data = fit_adj)

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
fit.s2bak <- function(formula, data, obs, surv, trait,
                        background = NA, sdm.fun = gam,
                        predict.fun = predict.gam,
                        overlapBackground = TRUE, nbackground = 10000,
                        addSurvey = TRUE,
                        ncores = 1,
                        readout = NA, version = c("full", "short")[1], ...) {
  # Output for fit.s2bak
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
  out$s2bak.SO <- fit.so(formula = formula, data = data, obs = obs,
                            background = background, sdm.fun = sdm.fun,
                            overlapBackground = overlapBackground,
                            nbackground = nbackground, ncores = ncores,
                            readout = readout, version = version, ...)

  ## Next, fit S2 model
  ## Only for species with survey data
  if (all(is.na(surv))) stop("Missing survey data.")

  out$s2bak.S2 <- fit.s2(formula = formula, data = data,
                            obs = obs,
                            surv = surv,
                            background = background,
                            sdm.fun = sdm.fun,
                            overlapBackground = overlapBackground,
                            nbackground = nbackground,
                            addSurvey = addSurvey,
                            ncores = ncores, readout = readout,
                            version = version, ...)

  ## Finally, fit BaK
  if (is.na(ind)) {
    #### ASSUMES 1ST COLUMN OF SURV IS ROWNUM ####
    surv_dat <- as.data.frame(data[surv[, 1], ])
  } else {
    wh <- which(colnames(data) == ind)
    surv_dat <- as.data.frame(data[match(surv[, ind], data[, ind]), ])
    # Move `ind` it to the first column
    surv_dat <- surv_dat[, c(ind, colnames(surv_dat)[-wh])]
  }

  # Make predictions using out$SO (always type = "response")
  ### THIS MIGHT BE A PROBLEM WITH OTHER PREDICT FUNCTIONS!! ####
  predictions <- s2bak.predict.SOS2(out$s2bak.SO, surv_dat,
                                      predict.fun = predict.fun,
                                      useReadout = !is.na(readout),
                                      ncores = ncores, type = "response")

  out$s2bak.BaK <- fit.bak(predictions, surv, surv_dat, trait)

  return(out)
}

#' @title Combine SO, S2 and BaK objects
#'
#' @description Combines separately fitted SO, S2 and BaK models into a single
#' S2BaK object. The output will be identical to running the entire process in
#' fit.s2bak function, and can be used in the same situations.
#'
#' Models can be partially provided, for instance only SO and BaK, in which case
#' the use of s2bak predict functions will not apply S2 but instead only run
#' predictions with sightings-only and BaK adjustment.
#'
#' @param so output from fit.so or fit.s2 without survey data.
#' @param s2 output from fit.s2 function
#' @param bak output from s2back.BaK function
#' @return Object of class s2bak.S2Bak, equivalent to having run
#' \link[s2bak]{fit.s2bak} for the entire dataset
#' @export
combine.s2bak <- function(so = NULL, s2 = NULL, bak = NULL) {
  out <- list(
    s2bak.SO = so,
    s2bak.S2 = s2,
    s2bak.BaK = bak
  )
  class(out) <- "s2bak.S2BaK"
  return(out)
}
