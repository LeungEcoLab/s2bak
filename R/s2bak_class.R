## Defining an S4 class for s2bak.so, s2bak.s2 and s2bak.bak
#### WORK IN PROGRESS ####
## THESE SHOULD NOT HAVE THE SAME NAME AS THE FITTING FUNCTIONS!!!

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
                                  options = "list",
                                  sdm = "list",
                                  failure = "character"
                        ))

  # Note sure if this one works the way I think... WIP
  s2bak <- setClass(Class = "s2bak",
                    slots = c(s2bak.SO = "s2bak.so",
                              s2bak.S2 = "s2bak.s2",
                              s2bak.BaK = "s2bak.bak"
                    ))

  ## Show methods: what happens when you just call an object
  setMethod(f = "show",
            signature = "s2bak.so",
            definition = function(object){
              cat("Sightings-only (SO) model\n")
              cat("\tSpecies list:", object$speciesList, "\n")
              cat("WORK IN PROGRESS!\n")
            })

  setMethod(f = "show",
            signature = "s2bak.s2",
            definition = function(object){
              cat("S2 fitted model\n")
              cat("\tSpecies list:", object$speciesList, "\n")
              cat("WORK IN PROGRESS!\n")
            })

  setMethod(f = "show",
            signature = "s2bak.bak",
            definition = function(object){
              cat("BaK fitted model\n")
              cat("\tSpecies list:", object$speciesList, "\n")
              cat("WORK IN PROGRESS!\n")
            })

  setMethod(f = "show",
            signature = "s2bak",
            definition = function(object){
              cat("S2BaK fitted model\n")
              cat("\tSpecies list:", object$speciesList, "\n")
              cat("WORK IN PROGRESS!\n")
            })
}
