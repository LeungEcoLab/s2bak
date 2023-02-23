## Defining an S4 class for s2bak.so, s2bak.s2 and s2bak.bak

#' @title Display information for s2bak.so and s2bak.s2 class object
#'
#' @description Display information for s2bak.so and s2bak.s2 class object
#'
#' @param object object of class `s2bak.so` or `s2bak.s2`
#' @return Printout of information
display_sos2 <- function(object) {
  if (object@empty) {
    cat("Empty", class(object), "class object\n")
  } else {
    if (typeof(object@options$call) == "language") {
      cat("Formula:", paste(deparse(object@options$call,
                                    width.cutoff = 500),
                            collapse = ""),
          "\n\n")
    } else {
      cat("Formula: list of", length(object@options$call),
          "formulas\n\n")
    }
    cat(length(object@speciesList), "species\n")
    cat(object@options$nbackground,
        tolower(object@options$background),
        ifelse(object@options$overlapBackground, "with", "no"),
        "overlap\n")
    if (length(object@failure) > 0) {
      cat("\t", length(object@failure), "failed species\n")
    }
    cat("Version:", object@options$version, "model object\n")
    cat("File output:", object@options$readout, "\n")
  }
}

#' @title Display information for s2bak.bak class object
#'
#' @description Display information for s2bak.bak class object
#'
#' @param object object of class `s2bak.bak`
#' @return Printout of information
display_bak <- function(object) {
  if (object@empty) {
    cat("Empty", class(object), "class object\n")
  } else {
    cat("Species bias formula:",
        paste(deparse(object@options$call$bias_species,
            width.cutoff = 500),
          collapse = ""),
        "\n\n")
    cat("Site bias formula:",
        paste(deparse(object@options$call$bias_site,
            width.cutoff = 500),
          collapse = ""),
        "\n\n")

    cat(length(object@speciesList), "survey species,",
    length(object@speciesListFull), "total species\n")
    cat(object@options$nsites, "survey sites\n\n")
  }
}

## Defining classes
s2bak.so <- setClass(Class = "s2bak.so",
                    slots = c(speciesList = "character",
                              options = "list",
                              sdm = "list",
                              failure = "character",
                              empty = "logical"
                    ),
                    prototype = prototype(empty = TRUE))

s2bak.s2 <- setClass(Class = "s2bak.s2",
                    slots = c(speciesList = "character",
                              options = "list",
                              sdm = "list",
                              failure = "character",
                              empty = "logical"
                    ),
                    prototype = prototype(empty = TRUE))

s2bak.bak <- setClass(Class = "s2bak.bak",
                      slots = c(speciesList = "character",
                                speciesListFull = "character",
                                options = "list",
                                bak = "list",
                                empty = "logical"
                    ),
                    prototype = prototype(empty = TRUE))

s2bak <- setClass(Class = "s2bak",
                  slots = c(s2bak.so = "s2bak.so",
                            s2bak.s2 = "s2bak.s2",
                            s2bak.bak = "s2bak.bak"))

## Show methods: what happens when you just call an object
setMethod(f = "show",
          signature = "s2bak.so",
          definition = function(object){
            cat("Sightings-only (s2bak.so) model(s)\n\n")
            display_sos2(object)
          })

setMethod(f = "show",
          signature = "s2bak.s2",
          definition = function(object){
            cat("Sightings-Survey (s2bak.s2) model(s)\n\n")
            display_sos2(object)
          })

setMethod(f = "show",
          signature = "s2bak.bak",
          definition = function(object){
            cat("Bias-adjustment kernel models\n\n")
            display_bak(object)
          })

setMethod(f = "show",
          signature = "s2bak",
          definition = function(object){
            cat("S2BaK fitted model\n")
            cat("----------\n\n")

            cat("Sightings-only (s2bak.so) model(s)\n\n")
            display_sos2(object@s2bak.so)
            cat("----------\n\n")

            cat("Sightings-Survey (s2bak.s2) model(s)\n\n")
            display_sos2(object@s2bak.s2)
            cat("----------\n\n")

            cat("Bias-adjustment kernel models\n\n")
            display_bak(object@s2bak.bak)
          })
