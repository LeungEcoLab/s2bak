#' @title Sample target-group background sites
#'
#' @description This function samples background sites used in species
#' distribution modelling (SDM), using the target-group background (TGB)
#' sampling method.
#'
#' Sites within the species observations, weighted by the number of sightins
#' in `obs`.
#'
#' without sightings will be excluded from sampling.
#' @param obs Species observations as a data.frame, with a column for species
#' name (must be labelled 'species') and column of index of observations to
#' reflect presences. If the index column name is not found in 'data', it
#' assumes row number.
#' @param nbackground The number of background sites to sample. If the number
#' of sites specified is less than the available sites in `data`, then the
#' maximum will be sampled.
#'
#' @examples
#' ## WIP ##
#' ## But here we would generate data, sample background sites then fit model
#'
#' @return A vector of indices corresponding to the `data` argument, which can
#' be used to fit SDMs through link[s2bak]{fit.s2bak},
#' \link[s2bak]{fit.s2bak.so} and \link[s2bak]{fit.s2bak.s2}
#' @export
tgb_sample <- function(obs, nbackground = 10000) {
  # Figure out index name. If not found then use first column of data
  # And rename obs' index name
  ind <- colnames(obs)[which(colnames(obs) != "species")]

  # Otherwise, we'll sample according to sightings
  # As data.frame, first column is index and second column is counts
  obs_count <- as.data.frame(table(obs[, ind]), stringsAsFactors = FALSE)

  # Convert to numeric if the original was numeric
  if (all(is.numeric(obs[, ind]))) {
    obs_count[, 1] <- as.numeric(obs_count[, 1])
  }

  if (nrow(obs_count) <= nbackground) {
    warning(paste("'nbackground' exceeds number of available sites; returning",
      nrow(obs_count), "sites\n"))
    return(obs_count[, 1])
  }

  # Otherwise, do weighted sampling
  obs_count$prob <- obs_count[, 2] / sum(obs_count[, 2])

  # Now sample
  bg_sample <- sample(x = obs_count[, 1],
                      size = nbackground,
                      prob = obs_count$prob)

  return(bg_sample)
}
