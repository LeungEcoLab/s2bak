## Functions for summary wrapper

#' @name summary.s2bak
#'
#' @title Generate summary statistics for s2bak objects
#'
#' @description This function will allow the user to call the `summary`
#' function for `s2bak`, `s2bak.so`, `s2bak.s2` and `s2bak.bak` class objects.
#'
#' This model assumes that the model used in the fitting process has an
#' associated `summary` function associated with it.
#'
#' As these models may be fitted for multiple species or sub-models, the user
#' can specify which index or species to view.
#'
#' @param object Object of class s2bak`, `s2bak.so`, `s2bak.s2` or `s2bak.bak`.
#' @param species For `s2bak.so`, `s2bak.s2` and `s2bak` models.
#' The character name of the species to access the model,
#' or the index if an integer is provided.
#' @param ... Additional arguments provided by the object's `summary` function.
#'
#' @return The relevant summary output for a given species.
#'
#' @rdname summary.s2bak
#' @export
summary.s2bak.s2 <- function(object, species, ...) {
  return(summary(object@sdm[[species]], ...))

}

## SO is just a wrapper for s2, we call it so we only need to make changes once
## If necessary
#' @rdname summary.s2bak
#' @export
summary.s2bak.so <- function(object, species, ...) {
  return(summary.s2bak.s2(object, species, ...))

}

## BaK and S2BaK are slightly different, since BaK are three objects
## And S2BaK are collections of the three sub-models..
#' @param bak For `s2bak.bak` and `s2bak` models. Select which of the three
#' submodels to display, or all of them (default bak = "all").
#' @rdname summary.s2bak
#' @export
summary.s2bak.bak <- function(object,
                              bak = c("all", "bias_site",
                                        "bias_species", "bias_adj")[1],
                              ...) {
  if (bak == "all") {
    return(lapply(object@bak, summary, ...))
  } else {
    return(summary(object@bak[[bak]], ...))
  }
}

#' @param submodel Which submodel (`s2bak.so`, `s2bak.s2` or `s2bak.bak`)
#' to display (default submodel = "all", which displays all of them).
#' Note that species must be specified if attempting a summary of `s2bak.so`
#' or `s2bak.s2` models.
#' @rdname summary.s2bak
#' @export
summary.s2bak <- function(object,
                          submodel = c("all", "so", "s2", "bak")[1],
                            species = NULL,
                            bak = c("all", "bias_site",
                                      "bias_species", "bias_adj")[1],
                          ...) {

  if (submodel == "all") {
    return(list(
      summary.s2bak.so(object@s2bak.so, species, ...),
      summary.s2bak.s2(object@s2bak.s2, species, ...),
      summary.s2bak.bak(object@s2bak.bak, bak, ...)
    ))

  } else if (submodel == "so") {
    return(summary.s2bak.so(object@s2bak.so, species, ...))

  } else if (submodel == "s2") {
    return(summary.s2bak.s2(object@s2bak.s2, species, ...))

  } else {
    # BaK
    return(summary.s2bak.bak(object@s2bak.bak, bak, ...))

  }
}