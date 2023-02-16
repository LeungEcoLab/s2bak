## Defining an S4 class for s2bak.SO, s2bak.S2 and s2bak.BaK
#### WORK IN PROGRESS ####
## THESE SHOULD NOT HAVE THE SAME NAME AS THE FITTING FUNCTIONS!!!

if(0){
    ## Defining classes
    s2bak.SO <- setClass(Class = "s2bak.SO",
                        slots = c(speciesList = "character",
                                    options = "list",
                                    sdm = "list",
                                    failure = "character"
                                ))

    s2bak.S2 <- setClass(Class = "s2bak.S2",
                        slots = c(speciesList = "character",
                                    options = "list",
                                    sdm = "list",
                                    failure = "character"
                                ))

    s2bak.BaK <- setClass(Class = "s2bak.BaK",
                        slots = c(speciesList = "character",
                                    options = "list",
                                    sdm = "list",
                                    failure = "character"
                                ))

    # Note sure if this one works the way I think... WIP
    s2bak.S2BaK <- setClass(Class = "s2bak.S2BaK",
                        slots = c(s2bak.SO = "s2bak.SO",
                                    s2bak.S2 = "s2bak.S2",
                                    s2bak.BaK = "s2bak.BaK"
                                ))

    ## Show methods: what happens when you just call an object
    setMethod(f = "show",
            signature = "s2bak.SO",
            definition = function(object){
                cat("Sightings-only (SO) model\n")
                cat("\tSpecies list:", object$speciesList, "\n")
                cat("WORK IN PROGRESS!\n")
            })

    setMethod(f = "show",
            signature = "s2bak.S2",
            definition = function(object){
                cat("S2 fitted model\n")
                cat("\tSpecies list:", object$speciesList, "\n")
                cat("WORK IN PROGRESS!\n")
            })

    setMethod(f = "show",
            signature = "s2bak.BaK",
            definition = function(object){
                cat("BaK fitted model\n")
                cat("\tSpecies list:", object$speciesList, "\n")
                cat("WORK IN PROGRESS!\n")
            })

    setMethod(f = "show",
            signature = "s2bak.S2BaK",
            definition = function(object){
                cat("S2BaK fitted model\n")
                cat("\tSpecies list:", object$speciesList, "\n")
                cat("WORK IN PROGRESS!\n")
            })
}
