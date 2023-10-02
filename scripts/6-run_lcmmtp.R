#######################################################
#######################################################
### Run lcmmtp software on final dataset
### Kat Hoffman
#######################################################
#######################################################

library(tidyverse) 
library(lcmmtp) # remotes::install_github("nt-williams/lcmmtp")
library(mlr3superlearner) # remotes::install_github("nt-williams/mlr3superlearner")

all_wide <- read_rds("data/derived/all_wide.rds")

# set d prime and d star functions
d_ap <- function(data, trt) rep(1, length(data[[trt]]))
d_as <- function(data, trt) rep(0, length(data[[trt]]))

max_window <- 14
As <- paste0("A_", 1:max_window)
Ms <- paste0("M_", 1:max_window)
Ys <- paste0("Y_", 1:max_window)
Cs <- paste0("Observed_", 1:max_window)
Ls <- as.list(paste0("L_value_", 1:max_window, "_glucose"))
Zs <- as.list(paste0("Z_value_", 1:max_window, "_glucose"))

vars <- lcmmtp:::lcmmtp_variables$new(
  L = Ls,
  A = As,
  Z = Zs,
  M = Ms,
  Y = Ys,
  cens = Cs
)
# 
# # An `lcmmtp_variables` object mapping observed variables to the assumed variable structure.
# vars <- lcmmtp:::lcmmtp_variables$new(
#   L = list(c("L_value_1_glucose"), c("L_value_2_glucose")),
#   A = c("A_1", "A_2"),
#   Z = list(c("Z_value_1_glucose"), c("Z_value_2_glucose")),
#   M = c("M_1", "M_2"),
#   Y = c("Y_1", "Y_2"),
#   # Y = "Y_2",
#   cens = c("Observed_1", "Observed_2")
# )

lrnrs <- c("glm", "glmnet")
lrnrs <- c("glm")
fit <- lcmmtp(all_wide,
              vars,
              d_ap,
              d_as,
              .lcmmtp_control(learners_trt = lrnrs,
                              learners_mediator = lrnrs,
                              learners_QL = lrnrs,
                              learners_QZ = lrnrs,
                              learners_QM = lrnrs)
)

fit
