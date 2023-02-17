#' @title Simulate example data for fitting SDMs and S2BaK
#'
#' @description This function simulates the data needed for fitting SDMs for
#' multiple species including generated trait and environments. The function is
#' to primarily to illustrate the data structure, as well as demonstration
#' of the package functions.
#'
#' @param species numeric. The number of species used for fitting
#' @param sites numeric. The number of environment sites in the dataset
#' for species sightings
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
#' @return a list of data.frames for the environment of the sightings,
#' the environment of the surveys, species traits, species sightings, and
#' the survey presences.
#' @export
s2bakSim <- function(species = 30, sites = 3000, surv = 300,
                     minimum.sightings = 30,
                     env.prob = c(0.5, 1), trait.prob = c(0.5, 1),
                     max.iter = 100) {

  if (length(env.prob) != 2 || length(trait.prob) != 2) {
    stop(paste("Invalid range for probabilities.",
               "Ensure that two numeric values are provided."))
  }

  # Ensure that minimum is the first index and maximum is the second
  env.prob <- sort(env.prob)
  trait.prob <- sort(trait.prob)

  # Create four environments
  # for observations and survey sites
  envs <- lapply(1:4, rnorm, n = sites + surv)

  # To store (known) traits
  vtraits <- list()
  # Create species-specific data (sightings and survey)
  # Survey matrix that we will rbind later
  # Note that this means that if a species is not presence,
  # then it is an absence
  surv_presences <- list()
  # Sightings data.frame that we will rbind later
  obs_presences <- list()

  for (i in 1:species) {
    # Generate species repeatedly until obs >= minimum.sightings
    count <- 1
    go <- TRUE
    while (go && count <= max.iter) {
      # defines where the species is best suited
      m <- runif(2, -2, 2)
      traits <- rnorm(2)
      # standard deviation (i.e., how much of a specialist the species it)
      sdv <- runif(1, .5, 2)

      # Probability of detection based on environments 1 and 3
      p_d_e <- (env.prob[2] - env.prob[1]) *
        dnorm(dst(envs[[1]], envs[[3]])) /
        dnorm(0) + env.prob[1]

      # Probability of detection based on traits
      p_d_t <- (trait.prob[2] - trait.prob[1]) *
        dnorm(dst(traits[1], traits[2])) /
        dnorm(0) + trait.prob[1]

      # Probability of establishment based on environment 1 and 2
      # p_e maximum value is always 1, may want to relax that later
      p_e <- dnorm(dst(envs[[1]] - m[1], envs[[2]] - m[2]), sd = sdv) /
        dnorm(0, sd = sdv)

      # Generate binary based on probability of establishment
      # Therefore the 'true' presence-absence
      t_inv <- rbinom(sites + surv, 1, p_e)

      ## Probability of sighting
      # given probabilities of detection and probability of establishment
      p <- p_d_t * p_d_e * t_inv
      # This is to generate sightings
      obs <- rbinom(sites + surv, 1, p)

      # Subset obs to sites and t_inv to surv
      obs <- obs[1:sites]
      t_inv <- t_inv[(sites + 1):(sites + surv)]

      if (sum(obs) >= minimum.sightings && sum(t_inv) >= minimum.sightings) {
        vtraits[[i]] <- traits
        go <- FALSE
      }
    }

    if (count >= max.iter &&
        sum(obs) < minimum.sightings &&
        sum(t_inv) < minimum.sightings) {
      stop("Maximum iteration reached without producing valid species.")
    }

    # Format it according to our specs
    surv_presences[[i]] <- data.frame(
      species = paste0("Species", i),
      index = which(t_inv == 1)
    )

    # Observations
    obs_presences[[i]] <- data.frame(
      species = paste0("Species", i),
      index = which(obs == 1)
    )
  }

  # Generate the output of the function
  # Environments are separated by sites and surv
  output <- list(
    Environment_Sightings =  as.data.frame(cbind(
      1:sites, do.call(cbind, envs)[1:sites, ]
    )),
    Environment_Survey =  as.data.frame(cbind(
      1:surv, do.call(cbind, envs)[(sites + 1):(sites + surv), ]
    ))
  )
  colnames(output$Environment_Sightings) <- c("index",
                                                paste0("Environment", 1:4))
  colnames(output$Environment_Survey) <- c("index",
                                                paste0("Environment", 1:4))

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

  output$Sightings <- as.data.frame(do.call(rbind, obs_presences))
  output$Survey <-  as.data.frame(do.call(rbind, surv_presences))

  return(output)
}
