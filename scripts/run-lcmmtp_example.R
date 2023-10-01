# remotes::install_github("nt-williams/lcmmtp")
# remotes::install_github("nt-williams/mlr3superlearner")
library(tidyverse)
library(lcmmtp)
library(mlr3superlearner)

# Run test data

?lcmmtp::lcmmtp() # note this shows the package description, not the examples that show up in Github

data("lcmmtp_foo")

# set d prime and d star functions
d_ap <- function(data, trt) rep(1, length(data[[trt]]))
d_as <- function(data, trt) rep(0, length(data[[trt]]))

# An `lcmmtp_variables` object mapping observed variables to the assumed variable structure.
vars <- lcmmtp:::lcmmtp_variables$new(
    L = list(c("L_1"), c("L_2")),
    A = c("A_1", "A_2"),
    Z = list(c("Z_1"), c("Z_2")),
    M = c("M_1", "M_2"),
    Y = "Y",
    cens = c("c1", "c2")
)

fit <- lcmmtp(lcmmtp_foo,
       vars,
       d_ap,
       d_as
       )

fit

# modify data set so that if the observation is observed at c1, they have an outcome
lcmmtp_foo_mod <- lcmmtp_foo |>
  mutate(Y_mod = case_when(is.na(Y) & c1 == 1 ~ 0,
                       TRUE ~ Y
                       ))

# variables for 1 time point
tp1_vars <- lcmmtp:::lcmmtp_variables$new(
  L = list(c("L_1")),
  A = c("A_1"),
  Z = list(c("Z_1")),
  M = c("M_1"),
  Y = "Y_mod",
  cens = c("c1")
)

tp1_fit <- lcmmtp(lcmmtp_foo_mod,
              tp1_vars,
              d_ap,
              d_as,
              "glm",
              folds=2
)

#Error in `[.default`(M, complete.cases(M), ) : 
#incorrect number of dimensions


