# remotes::install_github("nt-williams/lcmmtp")
# remotes::install_github("nt-williams/mlr3superlearner")
library(tidyverse)
library(lcmmtp)
library(mlr3superlearner)

# Run test data

dat_t1_covars <- sample(t1_covars)

# set d prime and d star functions
d_ap <- function(data, trt) rep(1, length(data[[trt]]))
d_as <- function(data, trt) rep(0, length(data[[trt]]))

# An `lcmmtp_variables` object mapping observed variables to the assumed variable structure.
vars <- lcmmtp:::lcmmtp_variables$new(
  L = list(c("L1_potassium", "miss_L1_potassium","L1_ddimer", "miss_L1_ddimer")),
  A = c("A1"),
  Z = list(c("Z1_potassium", "miss_Z1_potassium","Z1_ddimer", "miss_Z1_ddimer")),
  M = c("M1"),
  Y = c("Y1"),
  cens = c("C1")
)


t1_covars <- as.data.frame(t1_covars)
fit <- lcmmtp(t1_covars,
              vars,
              d_ap,
              d_as,
              "glm",
              folds=2
)

fit

nrow(t1_covars["M1"][complete.cases(t1_covars["M1"])])



