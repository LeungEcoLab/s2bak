# `s2bak` R package

>Framework and initial code created by Brian Leung, with code adapted for package by Dat Nguyen and Brian Leung.

**THIS PACKAGE IS IN ITS EARLY STAGES AND A WORK IN PROGRESS. USE AT YOUR OWN RISK.**

This R package adapts the code used for the Sightings, Survey and Bias-adjustment Kernel (S2BaK) framework within Leung et al., 2019, a flexible, efficient and integrative bias-adjustment model to fit species distribution models (SDMs) for multiple species. The framework allows the user to apply any SDM approach to the framework, such as Generalized Additive Models (GAMs) or MaxEnt (Phillips et al., 2006).

**An explanation of S2BaK will be added to this README. Please refer to Leung et al., 2019 in the meanwhile.**

## Installation

The `s2bak` package, currently hosted on GitHub, can be installed using `devtools` R package:

```R
# install.packages("devtools")
library(devtools)
install_github("https://github.com/LeungLab2/s2bak")
library(s2bak)
```

## Package information

Four data objects are primarily required to fit the S2BaK model. **Making the data structures more flexible is a current goal**:

1. An environment data.frame with each row corresponding to a site, including spatial predictors of bias.
2. A species-trait data.frame with each row being a species, and each trait being a column.
3. A two-column data.frame of species sightings linked to the row-index of the environment data.
4. A matrix of binary presence (1) and absence (0) data, with site as rows, species as columns. The first column is assumed to be row-index.

Simulated sample data can be generated using the function `s2bakSim`, which is used for demonstration and to illustrate the data structure.

The functions `s2bak.S2`, `s2bak.SO` and `s2bak.S2BaK` fit SDMs for multiple species, but the user may specify a different SDM model than the default (`glm`). The package currently supports fitting SDM model functions that support formulae (e.g., `mgcv::gam` and `glm`). **Supporting the most popular SDM methods is a goal, such as MaxEnt**.

Given that the framework is intended to be used with multiple species, the package has integrated parallelization using the `doParallel` package. Within the functions `s2bak.S2`, `s2bak.SO` and `s2bak.S2BaK` functions, which are all used to fit SDMs, the user can specify the number of cores using the `ncores` parameter.

## References

Leung, B., Hudgins, E. J., Potapova, Anna, & Ruiz-Jaen, Maria. (2019). A new baseline for countrywide α-diversity and species distributions: illustration using >6000 plant species in Panama. Ecological Applications.

Microsoft Corporation and Steve Weston (2022). doParallel: Foreach Parallel Adaptor for the 'parallel' Package. R package version 1.0.17. https://CRAN.R-project.org/package=doParallel

Phillips, S. J., Anderson, R. P., & Schapire, R. E. (2006). Maximum entropy modeling of species geographic distributions. Ecological modelling, 190(3-4), 231-259.

Wood, S.N. (2011) Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. Journal of the Royal Statistical Society (B) 73(1):3-36