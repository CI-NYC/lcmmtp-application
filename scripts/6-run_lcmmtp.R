### Run lcmmtp software on final dataset
### Kat Hoffman

library(tidyverse) 
library(lcmmtp) # remotes::install_github("nt-williams/lcmmtp")
library(mlr3superlearner) # remotes::install_github("nt-williams/mlr3superlearner")

all_wide <- read_rds("data/derived/all_wide.rds")

# set d prime and d star functions
d_ap <- function(data, trt) rep(1, length(data[[trt]]))
d_as <- function(data, trt) rep(0, length(data[[trt]]))

test <-
  all_wide |>
    select(L_value_1_glucose,
           L_value_2_glucose,
           A_1,
           A_2,
           Z_value_1_glucose,
           Z_value_2_glucose,
           M_1,
           M_2,
           Y_1,
           Y_2,
           C_1,
           C_2
           )
  

# An `lcmmtp_variables` object mapping observed variables to the assumed variable structure.
vars <- lcmmtp:::lcmmtp_variables$new(
  L = list(c("L_value_1_glucose"), c("L_value_2_glucose")),
  A = c("A_1", "A_2"),
  Z = list(c("Z_value_1_glucose"), c("Z_value_2_glucose")),
  M = c("M_1", "M_2"),
  Y = c("Y_1", "Y_2"),
  cens = c("C_1", "C_2")
)

tmp <-
  test |>
  select(L_value_1_glucose,
         L_value_2_glucose,
         A_1, A_2,
         Z_value_1_glucose,
         Z_value_2_glucose,
         M_1, M_2,
         Y_1, Y_2,
         C_1, C_2)

tmp |>
  filter(C_1 == 0)


fit <- lcmmtp(test,
              vars,
              d_ap,
              d_as# ,
              # "glm",
              # folds=2
)

fit
