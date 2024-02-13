#######################################################
#######################################################
### Run lcmmtp software on final dataset
#######################################################
#######################################################

## devtools::install_github("nt-williams/lcmmtp")
## devtools::install_github("nt-williams/mlr3superlearner@Screener")
## remotes::install_github('mlr-org/mlr3extralearners@*release')

library(tidyverse)
##library(lcmmtp)
library(devtools)
load_all(path = "lcmmtp")
library(mlr3superlearner) ## remotes::install_github("nt-williams/mlr3superlearner")
library(tictoc)
library(doFuture)
library(mlr3extralearners)

options(error = recover)
plan(multisession, workers=10)

# cluster data paths
data_path <- ""
results_path <- ""

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

    data[data[[trt]] == 2, trt] <- 1

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

max_window <- 4
As <- paste0("A_", 1:max_window)
Ms <- paste0("M_", 1:max_window)
Ys <- paste0("Y_", 1:max_window)
Cs <- paste0("Observed_", 1:max_window)

L_names <- all_wide %>% # get all letters for time varying
  dplyr::select(starts_with("L_")) %>%
  names()

# turn time varying col names into a list (see ?lmtp_sdr tv argument input)
Ls <- map(1:max_window, function(x) { # list for time varying covars
  L_names[str_detect(L_names, paste0("_", x, "_"))]
})

Z_names <- all_wide %>% # get all letters for time varying
  dplyr::select(starts_with("Z_")) %>%
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

### COMMENT/UNCOMMENT TO DEBUG
## lrnrs <- c("mean")
folds <- 3
lrnrs <- c("mean", "earth", "cv_glmnet", "lightgbm")

## tic()
## fit_d1d2 <- lcmmtp(all_wide,
##                    vars,
##                    d1,
##                    d2,
##                    .lcmmtp_control(folds = folds,
##                                    learners_trt = lrnrs,
##                                    learners_mediator = lrnrs,
##                                    learners_QL = lrnrs,
##                                    learners_QZ = lrnrs,
##                                    learners_QM = lrnrs)
## )
## toc()

tic()
fit_d2d1 <- lcmmtp(all_wide,
                   vars,
                   d2,
                   d1,
                   .lcmmtp_control(folds = folds,
                                   learners_trt = lrnrs,
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
                   .lcmmtp_control(folds = folds,
                                   learners_trt = lrnrs,
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
                   .lcmmtp_control(folds = folds,
                                   learners_trt = lrnrs,
                                   learners_mediator = lrnrs,
                                   learners_QL = lrnrs,
                                   learners_QZ = lrnrs,
                                   learners_QM = lrnrs)
)
toc()


## write_rds(fit_d1d2, paste0(results_path, paste0("fit_d1d2_", max_window, "tp.rds")))
write_rds(fit_d2d1, paste0(results_path, paste0("fit_d2d1_", max_window, "tp.rds")))
write_rds(fit_d1d1, paste0(results_path, paste0("fit_d1d1_", max_window, "tp.rds")))
write_rds(fit_d2d2, paste0(results_path, paste0("fit_d2d2_", max_window, "tp.rds")))

## Survival is also higher under a delay (in addition to AKI being lower, see LIDA paper)

## devtools::install_github("nt-williams/lcmmtp")
## devtools::install_github("nt-williams/mlr3superlearner@Screener")
## remotes::install_github('mlr-org/mlr3extralearners@*release')

library(tidyverse)
##library(lcmmtp)
library(devtools)
load_all(path = "lcmmtp")
library(mlr3superlearner) ## remotes::install_github("nt-williams/mlr3superlearner")
library(tictoc)
library(doFuture)
library(mlr3extralearners)

## options(error = recover)
## plan(multisession, workers=10)
max_window <- 4

## # cluster data paths
## data_path <- ""

results_path <- ""

## ## fit_d1d2 <- read_rds(paste0(results_path, paste0("fit_d1d2_", max_window, "tp.rds")))
fit_d2d1 <- read_rds(paste0(results_path, paste0("fit_d2d1_", max_window, "tp.rds")))
fit_d1d1 <- read_rds(paste0(results_path, paste0("fit_d1d1_", max_window, "tp.rds")))
fit_d2d2 <- read_rds(paste0(results_path, paste0("fit_d2d2_", max_window, "tp.rds")))

total <- -(fit_d2d2$theta - fit_d1d1$theta)
direct <- -(fit_d2d1$theta - fit_d1d1$theta)
mediated <- -(fit_d2d2$theta - fit_d2d1$theta)

se.total <- sqrt(var(fit_d2d2$S - fit_d1d1$S)/dim(all_wide)[1])
se.direct <- sqrt(var(fit_d2d1$S - fit_d1d1$S)/dim(all_wide)[1])
se.mediated <- sqrt(var(fit_d2d2$S - fit_d2d1$S)/dim(all_wide)[1])

pseudo_med <- -(fit_d2d2$S - fit_d2d1$S)
pseudo_dir <- -(fit_d2d1$S - fit_d1d1$S)
pseudo <- (pseudo_dir + pseudo_med)

vars <- all_wide %>% select(starts_with(c('L_1_', 'L_value_1_'))) %>% select(!contains('miss')) %>% select(!L_1_smoking_no) %>% select(!L_1_red_cap_source) %>% select(!L_1_ethnicity_not_hispanic_or_latino_or_spanish_origin)  %>% select(!L_1_race_other)

options(error=NULL)
library(glmnet)
coefs_med <- data.frame(name = names(vars), coef = NA)
for(j in 1:ncol(vars)) coefs_med[j, 2] <- coef(lm(pseudo_med ~ ., data = vars[,j]))[2]
res_med <- coefs_med %>% arrange(desc(abs(coef))) %>% head(20)
names_med <- c('ILD', 'Cancer', 'COPD', 'Lymph. count', 'CVA', 'CAD',
               'Current smoker', 'Cirrhosis', 'Creatinine', 'Home O2',
               'HIV', 'Former smoker', 'Male', 'Asthma', 'Hispanic', 'Immunosuppressed',
               'Asian', 'Diabetes', 'Black', 'White')
res_med[, 1] <- names_med
library(xtable)
print(xtable(res_med[1:10,]), include.rownames=FALSE)

coefs_dir <- data.frame(name = names(vars), coef = NA)
for(j in 1:ncol(vars)) coefs_dir[j, 2] <- coef(lm(pseudo_dir ~ ., data = vars[,j]))[2]
res_dir <- coefs_dir %>% arrange(desc(abs(coef))) %>% head(20)
names_dir <- c('ILD', 'Creatinine','Asthma', 'CVA','Cirrhosis',
               'Immunosuppressed', 'Bilirubin', 'Former smoker',
               'Home O2','White', 'Asian','COPD','Hypertension',
               'Cancer', 'Black', 'Current smoker', 'Diabetes',
               'Hypoxia', 'Arterial PCO2', 'CAD')
res_dir[, 1] <- names_dir
library(xtable)
print(xtable(res_dir[1:10,]), include.rownames=FALSE)

