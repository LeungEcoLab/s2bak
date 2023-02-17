#' @name s2bak
#' @title Build sightings-only or S2 species distribution models for
#' multiple species.
#'
#' @description fit.s2bak.so function fits SDMs for each provided species within
#' the same system, using a specified SDM approach (or the default which are
#' GAMs from the mgcv package). Parallelization is possible when processing each
#' SDM, with the default being 1 core.
#'
#' The fit.s2bak.s2 function fits SDMs using species sightings,
#' background sites and survey sites, differentiating between them using a
#' binary 'so' predictor, denoting sightings-only (1) or survey (0).
#'
#' Saving SDMs to the output may be computationally intensive, particularly with
#' large datasets and many species. To reduce issues with memory,
#' readout and the version = "short" may be used, which does not output the
#' fitted models but instead saves it to the directory specified in readout.
#'
#' @param formula Formula for the model functions.
#' Assumes the structure follows "Y ~ X". Alternatively, a named list of
#' formulas can be provided corresponding to species names.
#' In this case, species will be fit using their corresponding formula.
#' The response variable can have any name, as the function name the column
#' accordingly.
#' @param data_obs A data.frame containing the covariates used for fitting
#' `s2bak.SO` and `s2bak.S2` models with sightings data.
#' The index of the data.frame linking sites to observations should
#' correspond to the indices in `obs`.
#' @param data_surv A data.frame containing the covariates used for fitting
#' `s2bak.S2` models with survey data. The index of the data.frame linking sites
#' to survey presences should correspond to the indices in `surv`. Default is
#' NA, as survey data is not necessary to fit `s2bak.SO` models.
#' @param obs A data.frame of species observations, with a column for species
#' name (must be labelled 'species') and column of index of observations to
#' reflect presences. If the index column name is not found in 'data', it
#' assumes row number.
#' @param surv A data.frame of species presences for the survey data used to
#' fit `s2bak.S2` models (optional otherwise), with a column for species
#' name (must be labelled 'species') and column of index of observations to
#' reflect presences. If the index column name is not found in 'data', it
#' assumes row number.
#' It will add the an additional binary predictor to the formula(s), `so`,
#' denoting whether a sites is sightings-only (1) or survey data (0).
#' If 'so' is already in formula (e.g. if modifying
#' the variable in any way, then set addSurvey = FALSE). If left as NA,
#' it will fit the SDMs as presence-only models with the function of choice.
#' NOTE: Assumes all species are provided in the survey data, and that they are
#' surveyed over the same sites (i.e. matrix-type).
#' @param sdm.fun Model (as function) used for fitting (the default is
#' \link[stats]{glm}). The function must have the formula
#' as their first argument, and 'data' as the parameter for the dataset
#' (including presences and background sites within the data.frame).
#' @param background Background sites (pseudo-absences) used to fit the
#' presence-only model, provided as a vector of indices of data (following the
#' same column name as observations). If the index column name is not found in
#' 'data', it assumes row number within 'data'. If left as NA, it will randomly
#' sample 'nbackground' sites, with or without overlap ('overlapBackground').
#' Currently, only one set of background sites can be used.
#' @param nbackground Number of background sites to sample. Only applies
#' if background = NA.
#' @param overlapBackground Whether sampled background sites that overlap with
#' observations should be included. By default, it allows overlap. If FALSE,
#' number of background sites may be less than specified or provided.
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
#' used in other functions such as \link[s2bak]{predict.s2bak.s2}.
#' @param version Whether the SDMs should be included in the output. With
#' "short", no the fitted SDMs are not provided. Setting to "full" (default)
#' will output the list with all SDMs. Setting to "short" and combined with
#' readout, can considerably reduce RAM usage while saving the progress so far,
#' which is useful when dealing with many species or large datasets.
#' @param ... Other arguments that are passed to the SDM function (sdm.fun).
#' @return An object of class "s2bak.s2", providing fitted SDMs for each species
#' based on the provided SDM modelling approach. The primary difference between
#' SO and S2 models are the additional data points from the survey data, and an
#' additional binary predictor 'so' which denotes whether the data is from
#' presence-background (1) or presence-absence data (0).
#' @rdname s2bak
#' @export
fit.s2bak.s2 <- function(formula,
                         data_obs, data_surv = NA,
                         obs, surv = NA,
                         sdm.fun = glm,
                         background = NA,
                         nbackground = 10000,
                         overlapBackground = TRUE,
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
    if (!(ind %in% colnames(data_obs))) {
      ind <- NA # Row number instead of column name
    }
  }

  # Are we fitting SO or S2?
  mode <- ifelse(all(is.na(surv)), "s2bak.so", "s2bak.s2")

  # Get the list of species in the data
  # (S2 models may have fewer species than SO)
  speciesListFull <- as.character(unique(obs$species))
  if (mode == "s2bak.so") {
    speciesList <- speciesListFull
  } else if (mode == "s2bak.s2") {
    speciesList <- as.character(unique(surv$species))
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
      stop(paste("Incorrect formula type.",
                 "Please provide a single formula or a list of",
                 "formulas with each index named after",
                 "the corresponding species.\n"))
    }
  } else {
    stop(paste("Invalid formula type. Please provide a single formula or a",
               "list of formulas with each index named after the",
               "corresponding species.\n"))
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
        background <- sample(1:nrow(data_obs), nbackground)
      } else {
        background <- sample(data_obs[, ind], nbackground)
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
    # Generate fitting data using 'data_obs', 'obs' and 'background'
    # If overlapBackground == FALSE, remove overlaps
    wh_i <- which(colnames(obs) != "species")
    if (!overlapBackground) {
      background <- background[!(background %in% obs[obs$species == i, wh_i])]
    }
    # Get the list of indices
    ind2 <- c(obs[obs$species == i, wh_i], background)
    # Row number
    if (is.na(ind)) {
      tmp_dat <- data_obs[ind2, ]
    } else {
      # Obs index name
      tmp_dat <- data_obs[match(ind2, data_obs[, ind]), ]
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
    # This is easier: if species present at site, then 1, otherwise 0
    if (!all(is.na(surv))) {
      # If overlapBackground == FALSE, remove overlaps
      wh_i <- which(colnames(surv) != "species")
      # Generate survey data
      tmp_surv <- as.data.frame(data_surv)
      tmp_surv$pa <- 0
      if (is.na(ind)) {
        tmp_surv$pa[surv[surv$species == i, wh_i]] <- 1
      } else {
        tmp_surv$pa[match(surv[surv$species == i, ind], tmp_surv[, ind])] <- 1
      }

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
fit.s2bak.so <- function(formula, data_obs, obs,
                         sdm.fun = glm, background = NA,
                         nbackground = 10000,
                         overlapBackground = TRUE,
                         ncores = 1,
                         readout = NA, version = c("full", "short")[1], ...) {
  return(fit.s2bak.s2(formula, data_obs, data_surv = NA, obs,
                      surv = NA, background = background,
                      sdm.fun = sdm.fun,
                      overlapBackground = overlapBackground,
                      nbackground = nbackground,
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
#' @param trait Full trait data for the species predictions, as a data.frame
#' with 'species' as a column and relevant traits for the remainder.
#' Like with the predictions, the species in the dataset do not necessarily
#' have to possess survey data, but will be used in the final adjustment model
#' as final output.
#' @return Bias adjustment models, the kernels (location and species),
#' as a second-order GLM.
#' @rdname s2bak
#' @export
fit.s2bak.bak <- function(predictions, data_surv, surv, trait) {
  # Check if we're missing trait species
  wh <- which(!(surv$species %in% trait$species))
  if (length(wh) > 0) {
    warning(paste("Species from survey data missing trait/predictions data and",
                  "excluded:", unique(surv$species[wh])))
    surv <- surv[-wh,]
  }
  speciesList <- unique(surv$species)

  ind <- colnames(surv)[which(colnames(surv) != "species")]

  # Model outputs
  out <- list(
    speciesList = speciesList,
    speciesListFull = unique(c(surv$species, trait$species))
  )
  class(out) <- "s2bak.bak"

  cat("Fitting bias adjustment kernel on", length(speciesList),
      "survey species and", length(out$speciesListFull), "species total.\n")

  # Throw error if predictions, surv and data are not aligned
  if (nrow(predictions) != nrow(data_surv)) {
    stop(paste("Differing rows for predictions and environment data:",
               "sites/rows must correspond with each other.\n"))
  }

  # Generate polynomial predictor variable names
  # Can replace with more complex models, if desired
  # (but right now it's GLM with this)
  tn <- colnames(trait)[colnames(trait) != "species"]
  en <- colnames(data_surv)[colnames(data_surv) != ind]
  numer_env <- unlist(lapply(data_surv[1, en], is.numeric))
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
  fit_l <- as.data.frame(data_surv)
  # Get summed predictions by site
  fit_l$so_only <- apply(predictions, 1, sum, na.rm = TRUE)
  # Get summed occurrences by site
  fit_l$pa <- 0
  sum_pres <- table(surv[, ind])
  ### ind ### ASUMES ROW NUMBER HERE.
  fit_l$pa[as.numeric(names(sum_pres))] <- sum_pres
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
      warning(paste("Missing trait data, which were imputed",
                    "using mean of the column.\n"))
      fit_sp[i, wh] <- msd[colnames(wh)]
    }
  }

  # Get summed predictions by species
  fit_sp$so_only <- apply(predictions[, speciesList], 2, sum, na.rm = TRUE)
  # Get summed occurrences by species
  fit_sp$pa <- 0
  sum_pres <- table(surv$species)
  ### ind ### ASUMES ROW NUMBER HERE.
  fit_sp$pa[match(names(sum_pres), fit_sp$species)] <- sum_pres

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
  # `x` = species name
  fit_adj <- lapply(speciesList, FUN = function(x) {
    ddf <- data.frame(
      loc = 1:nrow(data_surv), species = x, pa = 0,
      zso_only = NA, scale_so_l = NA, scale_so_sp = NA,
      stringsAsFactors = FALSE
    )
    ## ind ASSUMES ROW NUMBER HERE
    ddf$pa[surv[surv$species == x, ind]] <- 1

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
    ddf$scale_so_l <- predict(out$bak$bias_loc, data_surv)
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
fit.s2bak <- function(formula,
                      data_obs, data_surv = NA,
                      obs, surv = NA, trait,
                      sdm.fun = glm,
                      predict.fun = predict.glm,
                      background = NA,
                      nbackground = 10000,
                      overlapBackground = TRUE,
                      addSurvey = TRUE,
                      ncores = 1,
                      readout = NA, version = c("full", "short")[1], ...) {
  # Output for fit.s2bak
  out <- list()
  class(out) <- "s2bak"

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
  out$s2bak.SO <- fit.s2bak.so(formula = formula,
                               data_obs = data_obs,
                               obs = obs,
                               sdm.fun = sdm.fun,
                               background = background,
                               nbackground = nbackground,
                               overlapBackground = overlapBackground,
                               ncores = ncores,
                               readout = readout, version = version, ...)

  ## Next, fit S2 model
  ## Only for species with survey data
  if (all(is.na(surv))) stop("Missing survey data.")

  out$s2bak.S2 <- fit.s2bak.s2(formula = formula,
                               data_obs = data_obs,
                               data_surv = data_surv,
                               obs = obs,
                               surv = surv,
                               sdm.fun = sdm.fun,
                               background = background,
                               nbackground = nbackground,
                               overlapBackground = overlapBackground,
                               addSurvey = addSurvey,
                               ncores = ncores, readout = readout,
                               version = version, ...)

  # Make predictions using out$SO (always type = "response")
  ### THIS MIGHT BE A PROBLEM WITH OTHER PREDICT FUNCTIONS!! ####
  predictions <- predict.s2bak.s2(out$s2bak.SO, data_surv,
                                  predict.fun = predict.fun,
                                  useReadout = !is.na(readout),
                                  ncores = ncores, type = "response")

  out$s2bak.BaK <- fit.s2bak.bak(predictions, data_surv, surv, trait)

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
#' @param so output from fit.s2bak.so or fit.s2bak.s2 without survey data.
#' @param s2 output from fit.s2bak.s2 function
#' @param bak output from fit.s2bak.bak function
#' @return Object of class s2bak, equivalent to having run
#' \link[s2bak]{fit.s2bak}
#' @export
combine.s2bak <- function(so = NULL, s2 = NULL, bak = NULL) {
  out <- list(
    s2bak.SO = so,
    s2bak.S2 = s2,
    s2bak.BaK = bak
  )
  class(out) <- "s2bak"
  return(out)
}
