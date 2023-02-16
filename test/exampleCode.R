### Script to test examples for `s2bak` package
# If working, they will be imported into the `s2bak` documentation

rm(list = ls()); gc()
graphics.off()

setwd("../test/")

#### INSTALL `s2bak` ####

## Should re-install every time we want to test a new version
library(devtools)
library(roxygen2)
if (0) {
    tryCatch(detach("package:s2bak", unload = TRUE))
    install_github("https://github.com/LeungEcoLab/s2bak", force = TRUE)
}

#### TEST CODE ####

library(s2bak)

## Create dataset for 20 species (default is 30).
dat <- s2bakSim(species = 20)
## Subset survey data species to simulate an incomplete dataset
## This doesn't work the way we want to so it will likely be re-done
dat$Survey <- dat$Survey[, c(1, sample(2:11, 10))]

## First sightings-only model (`s2bak.SO`) using GAM and single core
## `s2bakSim` produces a smaller dataset, so we reduce `nbackground`
## Note that in this instance, background site indices are not provided,
## so they are randomly sampled from `data`
model_so <- fit.s2bak.so(pa ~ Environment1 + Environment2 +
                            Environment3 + Environment4,
                        data = dat$Environment,
                        obs = dat$Sightings,
                        ncores = 1,
                        sdm.fun = glm,
                        nbackground = 2000,
                        family = binomial
                    )

## We can also fit models using GAM
library(mgcv)
model_so_gam <- fit.s2bak.so(pa ~ s(Environment1) + s(Environment2) +
                            s(Environment3) + s(Environment4),
                        data = dat$Environment,
                        obs = dat$Sightings,
                        ncores = 1,
                        sdm.fun = mgcv::gam,
                        nbackground = 2000,
                        method = "GCV.Cp",
                        select = TRUE,
                        family = binomial
                    )

## Fitting using S2 (that is, including survey data into the SDM)
## This will limit it the number of models fit to the species with survey data!
model_s2 <- fit.s2bak.s2(pa ~ Environment1 + Environment2 +
                            Environment3 + Environment4,
                        data = dat$Environment,
                        obs = dat$Sightings,
                        surv = dat$Survey,
                        ncores = 1,
                        sdm.fun = glm,
                        nbackground = 2000,
                        family = binomial
                    )

## If we wanted to provide our own background data, we would provide the indices
bgsites <- sample(as.data.frame(dat$Environment)$index, 2000)
model_s2 <- fit.s2bak.s2(pa ~ Environment1 + Environment2 +
                            Environment3 + Environment4,
                        data = dat$Environment,
                        obs = dat$Sightings,
                        surv = dat$Survey,
                        ncores = 1,
                        sdm.fun = glm,
                        background = bgsites,
                        family = binomial
                    )

## Fit BaK using the presence-only model (`model.so`)
## Generate survey data (this seems a bit out of sorts..
## re-format to match SOS2?)
surv_env <- as.data.frame(dat$Environment[dat$Survey[, 1], ])
## Generate predictions on survey sites
## `predict` here is equivalent to `predict.s2bak.so`
so_preds <- predict(model_so, surv_env,
                    predict.fun = predict.glm,
                    ncores = 1, type = "response")

## Fit BaK using survey predictions generated from sightings-only model
## This is super unintuitive...
model_bak <- fit.s2bak.bak(so_preds, dat$Survey, surv_env, dat$Trait)

## We can combine all the three objects into a single s2bak object
## This produces the same output as using `s2bak.S2BaK` (see below)
model_s2bak <- combine.s2bak(model_so, model_s2, model_bak)

## All of the above is done step-by-step, but we can do all of it using
## `s2bak.S2BaK`, which calls fits the sightings-only and S2 SDMs as well as BaK
model_s2bak2 <- fit.s2bak(pa ~ Environment1 + Environment2 +
                                Environment3 + Environment4,
                            data = dat$Environment,
                            obs = dat$Sightings,
                            surv = dat$Survey,
                            trait = dat$Trait,
                            background = NA,
                            sdm.fun = glm,
                            predict.fun = predict.glm,
                            overlapBackground = TRUE,
                            nbackground = 2000,
                            addSurvey = TRUE,
                            ncores = 1,
                            readout = NA,
                            version = "full",
                            family = binomial
                        )

## We can compare the predictive results of BaK vs S2 vs SO for survey species
all_predictions <- predict()

## If memory issues arise, we can output the fitted models to file instead
## We can also specify "short" output for this reason
## Fitted models will not be stored in memory
newdir <- tempfile(pattern = "so", tmpdir = ".")
dir.create(newdir)
print(newdir)
model_so <- fit.so(pa ~ Environment1 + Environment2 +
                        Environment3 + Environment4,
                        data = dat$Environment,
                        obs = dat$Sightings,
                        ncores = 1,
                        sdm.fun = glm,
                        nbackground = 2000,
                        readout = newdir,
                        version = "short",
                        family = binomial
                    )

## Predictions can still be made in this way, using the same object
## useReadout forces files to be read
so_preds <- predict.s2bak.SOS2(model_so, surv_env,
                                predict.fun = predict.glm,
                                ncores = 1, useReadout = TRUE)

#### ETC ####

# For within-function testing purposes
if (0) {
    formula <- pa ~ Environment1 + Environment2 +
                    Environment3 + Environment4
    data <- dat$Environment
    obs <- dat$Sightings
    surv <- dat$Survey
    trait <- dat$Trait
    background <- NA
    sdm.fun <- glm
    predict.fun <- predict.glm
    overlapBackground <- TRUE
    nbackground <- 2000
    addSurvey <- TRUE
    ncores <- 1
    readout <- NA
    version <- "full"
    ### Although this won't work with `...`
    family <- binomial
}
