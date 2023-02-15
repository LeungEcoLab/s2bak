## Custom output for defined classes within `s2bak`
#### TESTING STUFF! ### DEAL WITH THIS LATER
## THESE SHOULD NOT HAVE THE SAME NAME AS THE FITTING FUNCTIONS!!!
if(0){
s2bak.SO <- setClass(Class = "s2bak.SO",
                    slots = c(species.list = "character",
                                options = "list",
                                sdm = "list",
                                failure = "character"
                                ))

s2bak.S2 <- setClass(Class = "s2bak.S2",
                    slots = c(species.list = "character",
                                options = "list",
                                sdm = "list",
                                failure = "character"
                                ))

setMethod(f = "show",
          signature = "s2bak.SO",
          definition = function(object){
              cat("Sightings-only (SO) model\n")
              cat("\tSpecies list:", object$species.list, "\n")
              cat("WORK IN PROGRESS!\n")
          })

setMethod(f = "show",
          signature = "s2bak.S2",
          definition = function(object){
              cat("S2 fitted model\n")
              cat("\tSpecies list:", object$species.list, "\n")
              cat("WORK IN PROGRESS!\n")
          })

setMethod(f = "show",
          signature = "s2bak.BaK",
          definition = function(object){
              cat("BaK fitted model\n")
              cat("\tSpecies list:", object$species.list, "\n")
              cat("WORK IN PROGRESS!\n")
          })

setMethod(f = "show",
          signature = "s2bak.S2BaK",
          definition = function(object){
              cat("S2BaK fitted model\n")
              cat("\tSpecies list:", object$species.list, "\n")
              cat("WORK IN PROGRESS!\n")
          })
}
