#######################################################
#######################################################
### Run lcmmtp software on final dataset
### Kat Hoffman
#######################################################
#######################################################

library(tidyverse) 
library(lcmmtp) # remotes::install_github("nt-williams/lcmmtp")
library(mlr3superlearner) # remotes::install_github("nt-williams/mlr3superlearner")
library(tictoc)
library(doFuture)

# plan(multisession, workers=50)

# cluster data paths
data_path <- "data/derived/"
results_path <- "results/"

all_wide <- read_rds(paste0(data_path, "all_wide.rds"))

### COMMENT/UNCOMMENT THE FACTORS -- mlr3 currently can't handle?
# set d prime and d star functions
d1 <- function(data, trt) {
  # extract time point
  tau <- readr::parse_number(trt)
  # get the col name of previous trt
  trt_prev <- paste0("A_", tau - 1)
  
  if(trt == "A_1") {
    # if first time point and intubated, set to 1
    # Could switch this to leave as intubated (if pos. violations for these immediately intubated people)
    
    # data[data[[trt]] == "2", trt] <- factor("1", levels=0:2)
    
    data[data[[trt]] == "2", trt] <- 1
    
  } else {
    # if intubated at time T but not T-1, set to 1
    # data[which(data[[trt]] == "2" & data[[trt_prev]] != "2"), trt] <- factor("1", levels=0:2)
    data[which(data[[trt]] == 2 & data[[trt_prev]] != 2), trt] <- 1
  }
  
  return(data[[trt]])
  # return(factor(data[[trt]], levels = 0:2)) # return the refactored treatment level
}

d2 <- function(data, trt) {
  return(data[[trt]]) # return the refactored treatment level
  # return(factor(data[[trt]], levels = 0:2))
} 

### COMMENT THIS FROM 2 to 8
# max_window <- 8
max_window <- 2
As <- paste0("A_", 1:max_window)
Ms <- paste0("M_", 1:max_window)
Ys <- paste0("Y_", 1:max_window)
Cs <- paste0("Observed_", 1:max_window)

L_names <- all_wide |> # get all letters for time varying
  dplyr::select(starts_with("L_")) |>
  names()

# turn time varying col names into a list (see ?lmtp_sdr tv argument input)
Ls <- map(1:max_window, function(x) { # list for time varying covars
  L_names[str_detect(L_names, paste0("_", x, "_"))]
})

Z_names <- all_wide |> # get all letters for time varying
  dplyr::select(starts_with("Z_")) |>
  names()

# turn time varying col names into a list (see ?lmtp_sdr tv argument input)
Zs <- map(1:max_window, function(x) { # list for time varying covars
  Z_names[str_detect(Z_names, paste0("_", x, "_"))]
})

vars <- lcmmtp:::lcmmtp_variables$new(
  L = Ls,
  A = As,
  Z = Zs,
  M = Ms,
  Y = Ys,
  cens = Cs
)

# Comment/uncomment to debug
# An `lcmmtp_variables` object mapping observed variables to the assumed variable structure.
vars <- lcmmtp:::lcmmtp_variables$new(
  L = list(c("L_value_1_glucose"), c("L_value_2_glucose")),
  A = c("A_1", "A_2"),
  Z = list(c("Z_value_1_glucose"), c("Z_value_2_glucose")),
  M = c("M_1", "M_2"),
  Y = c("Y_1", "Y_2"),
  # Y = "Y_2",
  cens = c("Observed_1", "Observed_2")
)

### COMMENT/UNCOMMENT TO DEBUG
lrnrs <- c("glm")
# lrnrs <- c("glm", "glmnet")
# lrnrs <- c("glm", "glmnet", "earth", "lightgbm")

tic()
fit_d1d2 <- lcmmtp(all_wide,
                   vars,
                   d1,
                   d2,
                   .lcmmtp_control(learners_trt = lrnrs,
                                   learners_mediator = lrnrs,
                                   learners_QL = lrnrs,
                                   learners_QZ = lrnrs,
                                   learners_QM = lrnrs)
)
toc()

tic()
fit_d1d1 <- lcmmtp(all_wide,
                   vars,
                   d1,
                   d1,
                   .lcmmtp_control(learners_trt = lrnrs,
                                   learners_mediator = lrnrs,
                                   learners_QL = lrnrs,
                                   learners_QZ = lrnrs,
                                   learners_QM = lrnrs)
)
toc()

tic()
fit_d2d2 <- lcmmtp(all_wide,
                   vars,
                   d2,
                   d2,
                   .lcmmtp_control(learners_trt = lrnrs,
                                   learners_mediator = lrnrs,
                                   learners_QL = lrnrs,
                                   learners_QZ = lrnrs,
                                   learners_QM = lrnrs)
)
toc()

write_rds(fit_d1d2, paste0(results_path, "fit_d1d2_8tp.rds"))
write_rds(fit_d1d1, paste0(results_path, "fit_d1d1_8tp.rds"))
write_rds(fit_d2d2, paste0(results_path, "fit_d2d2_8tp.rds"))
