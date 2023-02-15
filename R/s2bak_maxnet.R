## Lacks predict function right now! So let's not export it for now...
## We should consider this after everything is in place...

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
#'
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