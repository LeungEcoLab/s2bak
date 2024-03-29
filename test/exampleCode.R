### Script to test examples for `s2bak` package
# If working, they will be imported into the `s2bak` documentation

rm(list = ls())
gc()
graphics.off()

setwd("../test/")

#### INSTALL `s2bak` ####

if (0) {
  ## Should re-install every time we want to test a new version
  library(devtools)
  library(roxygen2)
  tryCatch(detach("package:s2bak", unload = TRUE))
  install_github("https://github.com/LeungEcoLab/s2bak", force = TRUE)
}

#### TEST CODE ####

library(s2bak)

## Create dataset for 20 species (default is 30).
dat <- s2bakSim(species = 20)
## Subset survey data species to simulate an incomplete dataset
survey_species <- sample(unique(dat$Survey$species), 10)
dat$Survey <- dat$Survey[dat$Survey$species %in% survey_species, ]

## First sightings-only model (`s2bak.SO`) using GAM and single core
## `s2bakSim` produces a smaller dataset, so we reduce `nbackground`
## Note that in this instance, background site indices are not provided,
## so they are randomly sampled from `data`
model_so <- fit.s2bak.so(pa ~ Environment1 + Environment2 +
                           Environment3 + Environment4,
                         data_obs = dat$Environment_Sightings,
                         obs = dat$Sightings,
                         sdm.fun = glm,
                         ncores = 1,
                         nbackground = 2000,
                         family = binomial
)

## We can also fit using other models, such as mgcv::gam
library(mgcv)
model_so_gam <- fit.s2bak.so(pa ~ s(Environment1) + s(Environment2) +
                               s(Environment3) + s(Environment4),
                             data_obs = dat$Environment_Sightings,
                             obs = dat$Sightings,
                             sdm.fun = mgcv::gam,
                             ncores = 1,
                             nbackground = 2000,
                             method = "GCV.Cp",
                             select = TRUE,
                             family = binomial
)

## Fitting using S2 (that is, including survey data into the SDM)
## `survey_var` should be included in the formula, or set addSurvey = TRUE
## Which will simply add so as an additional predictor
model_s2 <- fit.s2bak.s2(formula = pa ~ Environment1 + Environment2 +
                           Environment3 + Environment4 + so,
                         data_obs = dat$Environment_Sightings,
                         data_surv = dat$Environment_Survey,
                         obs = dat$Sightings,
                         surv = dat$Survey,
                         ncores = 1,
                         sdm.fun = glm,
                         nbackground = 2000,
                         survey_var = "so",
                         addSurvey = FALSE,
                         family = binomial
)

# We can specify GAM with SO as an interaction, meaning that survey-sightings
# will have a different effect
model_s2_gam <- fit.s2bak.s2(formula = pa ~ s(Environment1, so) +
                               s(Environment2, so) +
                               s(Environment3, so) + s(Environment4, so),
                         data_obs = dat$Environment_Sightings,
                         data_surv = dat$Environment_Survey,
                         obs = dat$Sightings,
                         surv = dat$Survey,
                         ncores = 1,
                         sdm.fun = mgcv::gam,
                         nbackground = 2000,
                         survey_var = "so",
                         addSurvey = FALSE,
                         family = binomial
)

## If we wanted to provide our own background data, we would provide the indices
## For instance with target-group background sampling
bgsites <- tgb_sample(obs = dat$Sightings, nbackground = 2000)
model_s2_tgb <- fit.s2bak.s2(pa ~ Environment1 + Environment2 +
                           Environment3 + Environment4 + so,
                         data_obs = dat$Environment_Sightings,
                         data_surv = dat$Environment_Survey,
                         obs = dat$Sightings,
                         surv = dat$Survey,
                         ncores = 1,
                         sdm.fun = glm,
                         background = bgsites,
                         family = binomial
)

## Fit BaK using the presence-only model (`model.so`)

## Generate predictions on survey sites
## `predict` here is equivalent to `predict.s2bak.so`
so_preds <- predict(model_so,
                    newdata = dat$Environment_Survey,
                    predict.fun = predict.glm,
                    ncores = 1, type = "response")

## Fit BaK using survey predictions generated from sightings-only model
model_bak <- fit.s2bak.bak(bias_site ~ Environment1 + Environment2 +
                                        Environment3 + Environment4 +
                                        I(Environment1^2) + I(Environment2^2) +
                                        I(Environment3^2) + I(Environment4^2),
                          bias_species ~ Trait1 + Trait2 +
                                          I(Trait1^2) + I(Trait2^2),
                      so_preds, dat$Environment_Survey,
                      dat$Survey, dat$Trait,
                      bak.fun = glm,
                      predict.bak.fun = predict.glm,
                      bak.arg = list(family = gaussian))

## We can combine all the three objects into a single s2bak object
## This produces the same output as using `s2bak.S2BaK` (see below)
model_s2bak <- combine.s2bak(model_so, model_s2, model_bak)

# Values can be partially combined, and also we can also replace models using
# replace.s2bak
model_s2bak <- combine.s2bak(model_so, model_s2)
model_s2bak <- replace.s2bak(model_bak, model_s2bak)

## All of the above is done step-by-step, but we can do all of it using
## `s2bak.S2BaK`, which calls fits the sightings-only and S2 SDMs as well as BaK
model_s2bak2 <- fit.s2bak(formula = pa ~ Environment1 + Environment2 +
                            Environment3 + Environment4,
                          formula_survey = pa ~ Environment1 + Environment2 +
                            Environment3 + Environment4 + so,
                          formula_site = bias_site ~ Environment1 +
                                        Environment2 + Environment3 +
                                        Environment4 +
                                        I(Environment1^2) + I(Environment2^2) +
                                        I(Environment3^2) + I(Environment4^2),
                          formula_species = bias_species ~ Trait1 + Trait2 +
                                          I(Trait1^2) + I(Trait2^2),
                          data_obs = dat$Environment_Sightings,
                          data_surv = dat$Environment_Sightings,
                          obs = dat$Sightings,
                          surv = dat$Survey,
                          trait = dat$Trait,
                          background = NA,
                          sdm.fun = glm,
                          predict.fun = predict.glm,
                          bak.fun = glm,
                          predict.bak.fun = predict.glm,
                          overlapBackground = TRUE,
                          nbackground = 2000,
                          survey_var = "so",
                          addSurvey = FALSE,
                          ncores = 1,
                          readout = NA,
                          version = "full",
                          family = binomial
)

## We can obtain predictions for all modelling approaches
## by specifying output = "all"
all_predictions <- predict(model_s2bak2,
                                  dat$Environment_Sightings,
                                  trait = dat$Trait,
                                  predict.fun = predict.glm,
                                  predict.bak.fun = predict.glm,
                                  output = "all", type = "response"
)

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
so_preds <- predict(model_so, surv_env,
                               predict.fun = predict.glm,
                               ncores = 1, useReadout = TRUE)

#### ETC ####

# For within-function testing purposes
if (0) {
  formula <- pa ~ Environment1 + Environment2 +
    Environment3 + Environment4
  data_obs <- dat$Environment_Sightings
  data_surv <- dat$Environment_Survey
  obs <- dat$Sightings
  surv <- dat$Survey
  trait <- dat$Trait
  sdm.fun <- glm
  predict.fun <- predict.glm
  background <- NA
  nbackground <- 2000
  overlapBackground <- TRUE
  addSurvey <- TRUE
  ncores <- 1
  readout <- NA
  version <- "full"
  ### Although this won't work with `...`
  family <- binomial

  model <- model_s2bak2

}
