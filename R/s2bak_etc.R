
#' @title Distance function
#' @description Calculate distance between two values
#' @param x Numerical first value, order does not matter
#' @param y Numerical second value, order does not matter
#' @return Square root of the sum of the values squared
dst <- function(x, y) {
  return((x^2 + y^2)^.5)
}

#' @title Truncate values based on min and max input
#'
#' @description Truncate provided numerical values based on a defined minimum and maximum cutoff. Used in \link[s2bak]{s2bak.BaK} and \link[s2bak]{s2bak.predict.BaK}.
#'
#' @param x Numerical values to truncate
#' @param min Minimum value to truncate, where \code{x <= min} with be \code{min}. If left as NA it will use the minimum value of x.
#' @param max Maximum value to truncate, where \code{x >= max} with be \code{max}. If left as NA it will use the maximum value of x.
#' @return Truncated x values
s2bak.truncate <- function(x, min = NA, max = NA) {
  if (is.na(min)) min <- min(x)
  if (is.na(max)) max <- max(x)

  x[x <= min] <- min
  x[x >= max] <- max

  return(x)
}

#' @title Add predictor to the formula.
#'
#' @description Adds 'so' (or other specified predictor) to the formula as an additional predictor. The arugment 'formula' assumes the formula follows a "Y ~ X" format. The functino will not add the predictor if it is already in the formula.
#'
#' @param formula Formula or list of formulas to add 'pred' to
#' @param pred Predictor name to add to formula
#' @return A formula with added predictor
s2bak.addPred <- function(formula, pred = "so") {
  # Check if formula or list of formulas
  if (typeof(formula) == "language") {
    flong <- FALSE
    formula <- list(formula)
  } else if (typeof(formula) == "list") {
    flong <- TRUE
  } else {
    # Throw error in all other cases than list and language
    stop("Invalid formula. Please provide a formula or list of formulas.")
  }

  out <- lapply(formula, FUN = function(x) {
    if (!(pred %in% labels(terms(x)))) {
      ff <- as.character(x)
      ff[3] <- paste(ff[3], "+", pred)
      ff <- formula(paste(ff[2], ff[1], ff[3]))
      return(ff)
    } else {
      return(x)
    }
  })

  if (flong) {
    return(out)
  } else {
    return(out[[1]])
  }
}

#' @title Scale continuous values within data.frame.
#'
#' @description Scales provided columns to a mean of 0 and standard deviation of 1. Does not scale binary and categorical predictors. If providing projected data as well, then it will scale them using the mean and standard deviation of the inputted environment and output both as a list.
#'
#' @param data Values to scale, with column names indicating the variable
#' @param project Projected values to scale, using the mean and standard deviations of 'env'. If NA, then the function will only scale the environmental values provided
#' @param getScaleValues Return mean and standard deviations of each factor
#' @return Returns a data.frame of scaled columns for non-categorical and non-binary variables. If getScaleValues or project data is provided, returns as a list.
s2bak.scale <- function(data, project = NA, getScaleValues = FALSE) {
  numer_env <- unlist(lapply(data[1, ], is.numeric))

  if (!all(is.na(project))) {
    # Throw an error if there is a mismatch in environment names
    diff <- c(setdiff(colnames(data), colnames(project)), setdiff(colnames(project), colnames(data)))
    if (length(diff) > 0) {
      stop(paste("Column names do not align:", paste(diff, collapse = ", ")))
    }
  }

  if (getScaleValues) {
    sv <- matrix(nrow = 2, ncol = length(colnames(data)[numer_env]))
    colnames(sv) <- colnames(data[numer_env])
  }

  for (cname in colnames(data)[numer_env]) {
    m <- mean(data[, cname])
    s <- sd(data[, cname])
    if (getScaleValues) {
      sv[1:2, cname] <- c(m, s)
    }
    data[, cname] <- (data[, cname] - m) / s
    if (all(!is.na(project))) {
      project[, cname] <- (project[, cname] - m) / s
    }
  }

  if (all(is.na(project))) {
    return(data)
  } else {
    out <- list(Scaled.values = data)

    if (any(!is.na(project))) out$Projected.values <- project

    if (getScaleValues) out$ScaleMSD <- sv
  }
  return(out)
}
