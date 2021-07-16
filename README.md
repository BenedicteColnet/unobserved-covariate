This repository contains the code to reproduce the results from the paper __Generalizing a causal effect:  sensitivity analysis and missing covariates__.

In a few words: it illustrates a sensitivity method to assess the impact of a missing covariate when you want to generalize a treatment effect from a randomized controlled trial to a target population.

Typical data you need for those applications are:

- A randomized controlled trial (RCT)
- An observational sample

# How to use this repository

- The simulation section can be launched from scratch using the Rmarkdown `simulation-unobserved-covariates.Rmd`. 
- The semi-synthetic benchmark can also be launched from scratch using the `STAR.Rmd`. 
- The code for the medical data analysis is available in the folder `code_medical_application`, but with no data for privacy reasons.


The libraries used are:

```
library(ggplot2) 
library(MASS) 
library(tidyverse) 
library(matlib) 
library(stringi)
library(Rmisc) 
library(caret)
library(rpart)
library(ggrepel) 
library(metR)
library(lubridate)
library(grf)
```

The code runs on R version `3.6.1`.

Approximate time needed to reproduce the analyses on a standard desktop machine: one hour.


# Want to use the functions?
Feel free to open the `functions.R` document and use function for your own purpose. The function use a concatenated data frame of the RCT and the observational sample. So you need to have your covariates (with any columns names), the outcome `Y`, the treatment `A`, and an indicatrix of being in the RCT (`S=1`) or not.


# Any questions?

Feel free to submit an issue! (I'll make my best to answer quickly)
Note that I launched non parametric version of these estimators, and we could see to improve this code together if you are interested to.

Cheers,
Bénédicte

