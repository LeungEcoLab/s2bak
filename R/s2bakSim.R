#' @title Simulate example data for fitting SDMs and S2BaK
#'
#' @description This function simulates the data needed for fitting SDMs for
#' multiple species including generated trait and environments. The function is
#' to primarily to illustrate the data structure, as well as demonstration
#' of the provided functions.
#'
#' @param species numeric. The number of species used for fitting
#' @param sites numeric. The number of sites in the dataset
#' @param surv numeric. The number of survey sites
#' @param minimum.sightings numeric. The function will re-generate a new species
#' until the minimum number of sightings (and survey presences) is greater
#' than or equal to the provided value.
#' @param env.prob numeric vector of size two for the range of probabilities for
#' establishment given environment
#' @param trait.prob numeric vector for the range of probabilities for
#' establishment given species trait
#' @param max.iter numeric. Number of attempts to generate a given species if
#' the number of presences is less than the minimum.
#'
#' @examples
#' ## See ?fit.s2bak
#'
#' @return a list with a data.frame of the environment for each site, a
#' data.frame of species traits, a data.frame of species sightings, and
#' the survey data as a matrix.
#' @export
s2bakSim <- function(species = 30, sites = 3000, surv = 300,
                     minimum.sightings = 30,
                     env.prob = c(0.5, 1), trait.prob = c(0.5, 1),
                     max.iter = 100) {
  if (length(env.prob) != 2 || length(trait.prob) != 2) {
    stop(paste("Invalid range for probabilities.",
               "Ensure that two numeric values are provided."))
  }
  if (surv > sites) {
    stop(paste("Number of survey sites exceeds total number of sites:",
               "surv must be less than sites"))
  }

  # Ensure that minimum is the first index and maximum is the second
  env.prob <- sort(env.prob)
  trait.prob <- sort(trait.prob)

  l <- 1:sites

  # Create four environments
  envs <- lapply(1:4, rnorm, n = sites)

  # Generate the output of the function
  # We already have the information for the 1/4 of it (environment)
  output <- list(Environment = cbind(
    l, do.call(cbind, envs)
  ))
  colnames(output$Environment) <- c("index", paste0("Environment", 1:4))

  vtraits <- list()

  # Create species-specific data (sightings and survey)
  # Survey matrix that we will cbind later
  survey_matrix <- list()
  # Sightings data.frame that we will rbind later
  sightings_data <- list()

  for (i in 1:species) {
    # Generate species repeatedly until obs >= minimum.sightings
    count <- 1
    go <- TRUE
    while (go && count <= max.iter) {
      # defines where the species is best suited
      m <- runif(2, -2, 2)
      traits <- rnorm(2)
      # trait that affects detection, unknown to the user
      unk_traits <- rnorm(2)
      # standard deviation (i.e., how much of a specialist the species it)
      sdv <- runif(1, .5, 2)

      # Probability of detection based on environments 1 and 3
      p_d_e <- (env.prob[2] - env.prob[1]) *
        dnorm(dst(envs[[1]], envs[[3]])) /
        dnorm(0) + env.prob[1]

      # Probability of detection based on traits
      p_d_t <- (trait.prob[2] - trait.prob[1]) *
        dnorm(dst(traits[1], unk_traits[2])) /
        dnorm(0) + trait.prob[1]

      # Probability of establishment based on environment 1 and 2
      # p_e maximum value is always 1, may want to relax that later
      p_e <- dnorm(dst(envs[[1]] - m[1], envs[[2]] - m[2]), sd = sdv) /
        dnorm(0, sd = sdv)

      # Generate binary based on probability of establishment
      # Therefore the 'true' presence-absence
      t_inv <- rbinom(sites, 1, p_e)

      ## Probability of sighting
      # given probabilities of detection and probability of establishment
      p <- p_d_t * p_d_e * t_inv
      # This is to generate sightings
      obs <- rbinom(sites, 1, p)

      if (sum(obs) >= minimum.sightings && sum(t_inv) >= minimum.sightings) {
        vtraits[[i]] <- traits
        go <- FALSE
      }
    }

    if (count >= max.iter) {
      stop("Maximum iteration reached without producing valid species.")
    }

    # Format it according to our specs
    survey_matrix[[i]] <- t_inv[1:surv]

    # Observations
    ### SHOULD THE FIRST `surv` ROWS BE EXCLUDED THOUGH?
    sightings_data[[i]] <- data.frame(
      species = paste0("Species", i),
      index = which(obs == 1)
    )
  }

  # We can make the data.frame for traits now
  output$Trait <- data.frame(
    species = paste0("Species", 1:species),
    Trait1 = sapply(vtraits, FUN = function(x) {
      x[1]
    }),
    Trait2 = sapply(vtraits, FUN = function(x) {
      x[2]
    })
  )

  output$Sightings <- do.call(rbind, sightings_data)
  survey_matrix <- do.call(cbind, survey_matrix)
  colnames(survey_matrix) <- paste0("Species", 1:species)
  output$Survey <- cbind(index = 1:surv, survey_matrix)

  return(output)
}
