#######################################################
########### Wrangle data for NEVERs
########### Kat Hoffman
#######################################################

# setup -------------------------------------------------------------------

library(tidyverse)
# library(tidylog)
library(janitor)

source(here::here("scripts/0-functions.R"))

cp_dt <- read_rds(here::here("data/derived/dts_cohort.rds"))
#labs <- read_rds(here::here("data", "raw", "2020-08-17", "wcm_labs_filtered.rds"))
all_labs_clean <- read_rds(here::here("data/derived/all_labs_clean.rds")) |>
  mutate(name = snakecase::to_snake_case(name))

### Filter for patients who were never intubated or AKI
nevers <-
  cp_dt |>
  filter((is.na(A_time) & is.na(M_time))) 

# Clean data set for determining time window with empis, l_start, l_end, z_start, and z_end
nevers_clean <-
  nevers |>
  nest(.by = empi) |>
  mutate(data = map(data, function(x){
    x |> 
      mutate(max_time = pmax(Y_time, Cens_time, na.rm = T),
             max_windows = difftime(max_time, t1_start),
             window = 1
      ) |>
      complete(window = 1:ceiling(max_windows)) |>
      fill(t1_start) |>
      mutate(l_start = t1_start + days(window - 1),
             l_end = l_start + hours(12) - seconds(1),
             z_start = l_start + hours(12),
             z_end = l_start + days(1) - seconds(1)) |>
      select(window, l_start, l_end, z_start, z_end) |>
      filter(window != 0) # filter out two windows that are 0 due to time ordering issues
  }) 
  ) |>
  unnest(cols = data)

# Create a grid of all covariates and time windows to map through
map_grid <- expand.grid("name" = unique(all_labs_clean$name), "window" = 1:28)
all_labs <- map2(map_grid$name, map_grid$window, ~divide_labs(nevers_clean, .x, .y))

# put together all L_t and Z_t in long format and create missingness indicators
Ls_and_Zs <- 
  data.table::rbindlist(all_labs) |> # bind all rows together
  mutate(covar = snakecase::to_snake_case(as.character(covar))) |> # clean up covariate names for wide col names later
  mutate(L_missing = ifelse(is.na(L_value), 1, 0), # create indicators for missing L and Z
         Z_missing = ifelse(is.na(Z_value), 1, 0)) |>
  group_by(empi, covar) |> 
  arrange(empi, window) |>
  fill(L_value, Z_value, .direction = "down") |> # Last Observation Carried Forward, by empi and covariate
  mutate(L_value = ifelse(is.na(L_value), -99999, L_value), # if no observations before, fill in with -99999 (could switch to median value)
         Z_value = ifelse(is.na(Z_value), -99999, Z_value))

# fill in M = 0 for all the mediators (this group never gets AKI)
# fill in A = 0 (change for oxygenation data)
M_and_A <- 
  nevers_clean |>
  mutate(M = 0,
         A = 0)

# Note I'll probably need to carry these Y's through (deterministic)
Ys <-
  nevers_clean |>
  select(empi, window, l_start, z_end) |>
  left_join(nevers |> select(empi, Y_time, Cens_time)) |>
  # if died in this period, outcome Y is 1
  mutate(Y = case_when(Y_time %within% interval(l_start, z_end) ~ 1,
                       # if they were discharged in this period, no outcome observed (NA)
                       Cens_time %within% interval(l_start, z_end) ~ NA_real_,
                       # otherwise, outcome is 0
                       TRUE ~ 0))

# Note I'll need to fill in the rest of the censoring variables as 0 (or leave as NA?) once in long format
Cs <- nevers_clean |>
  select(empi, window, l_start, z_end) |>
  left_join(nevers |> select(empi, Y_time, Cens_time)) |>
  # if they were discharged in this period, no outcome observed
  mutate(C = case_when(Cens_time %within% interval(l_start, z_end) ~ 0,
                       # otherwise, observed
                       TRUE ~ 1))

# Turn into wide form data ------------------------------------------------

max_window_data <- 28

M_and_A_wide <- 
  M_and_A |>
  filter(window < max_window_data) |>
  pivot_wider(id_cols=empi,
              names_from = window,
              values_from = c(A, M))

Ls_and_Zs_wide <-
  Ls_and_Zs |>
  filter(window < max_window_data) |>
  as_tibble() |>
  pivot_wider(id_cols=empi,
              names_from = c(window, covar),
              values_from = c(L_value, L_missing, Z_value, Z_missing))

Ys_wide <-
  Ys |>
  filter(window < max_window_data) |>
  pivot_wider(id_cols=empi,
              names_from = window,
              values_from = Y,
              names_prefix = "Y_")

Cs_wide <-
  Cs |>
  filter(window < max_window_data) |>
  pivot_wider(id_cols=empi,
              names_from = window,
              values_from = C,
              names_prefix = "C_")

# merge wide data set and save data ---------------------------------------

nevers_wide <-
  reduce(list(nevers_clean, M_and_A_wide, Ls_and_Zs_wide, Ys_wide, Cs_wide),
  ~left_join(.x, .y))

saveRDS(nevers_wide, here::here("data/derived/nevers_wide.rds"))


