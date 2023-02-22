## Defining an S4 class for s2bak.so, s2bak.s2 and s2bak.bak

#' @title Display information for s2bak.so and s2bak.s2 class object
#'
#' @description Display information for s2bak.so and s2bak.s2 class object
#'
#' @param object object of class `s2bak.so` or `s2bak.s2`
#' @return Printout of information
display_sos2 <- function(object) {
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

if (0) {
  ## Defining classes
  s2bak.so <- setClass(Class = "s2bak.so",
                      slots = c(speciesList = "character",
                                options = "list",
                                sdm = "list",
                                failure = "character"
                      ))

  s2bak.s2 <- setClass(Class = "s2bak.s2",
                      slots = c(speciesList = "character",
                                options = "list",
                                sdm = "list",
                                failure = "character"
                      ))

  s2bak.bak <- setClass(Class = "s2bak.bak",
                        slots = c(speciesList = "character",
                                  speciesListFull = "character",
                                  options = "list",
                                  bak = "list"
                        ))

  s2bak <- setClass(Class = "s2bak",
                    slots = c(s2bak.so = "s2bak.so",
                              s2bak.s2 = "s2bak.s2",
                              s2bak.bak = "s2bak.bak"
                    ))

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
              cat("Sightings-Survey (s2bak.s2) model(s)\n")
              display_sos2(object)
            })

  setMethod(f = "show",
            signature = "s2bak.bak",
            definition = function(object){
              cat("Bias-adjustment kernel models\n")
              cat(length(object@speciesList), "survey species\n")
              cat(length(object@speciesListFull), "total species\n")

              cat("WORK IN PROGRESS\n")
              cat("SPECIES FUNCTION CALL\n")
              cat("LOCATION FUNCITON CALL\n")

            })

  setMethod(f = "show",
            signature = "s2bak",
            definition = function(object){
              cat("S2BaK fitted model\n")

              cat("WORK IN PROGRESS\n")
              cat("DISPLAY EVERYTHING?\n")

            })
}
