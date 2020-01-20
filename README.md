# SPACancer
R package for Poisson regression of background mutation rate.

# Installation
Install directly from github with the R package devtools:
```r
install.packages("devtools")
library(devtools)
install_github("huangweiqing1993/SPACancer")
```

# Small usage example
```r
library(SPACancer)

# Load the example dataset:
data(mutation)

# Fit a mixed-effects Poisson regression model
fit <- ProBMR(m, M, x)

# Get background mutation rates
bmr <- fit$bkgd2$bmr

# Single gene analysis
convolution_test(bmr, n)
```
