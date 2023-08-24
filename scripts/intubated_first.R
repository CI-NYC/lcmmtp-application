#######################################################
########### Wrangle data for NEVERs
########### Kat Hoffman
#######################################################

# setup -------------------------------------------------------------------

library(tidyverse)
# library(tidylog)
library(janitor)

cp_dt <- read_rds(here::here("data/derived/dts_cohort.rds"))
#labs <- read_rds(here::here("data", "raw", "2020-08-17", "wcm_labs_filtered.rds"))
all_labs_clean <- read_rds(here::here("data/derived/all_labs_clean.rds")) |>
  mutate(name = snakecase::to_snake_case(name))

### Find patients for whom intubated before AKI (A_t < M_t)
intubated_first <- 
  cp_dt |>
  filter(A_time < M_time | (!is.na(A_time) & is.na(M_time))) |>
  filter(A_time >= t1_start)

### CREATE 24 HOUR PERIODS FROM ATIME TO T START
### THEN ATIME TO max time
  


create_periods  <- function(empi, t1_start, A_time, max_window_int, max_time, max_window_after_int){
  
  max_window_int <- max_window_int + 1
  
  periods_before <- 
    seq.POSIXt(from = t1_start, to = A_time, length.out = max_window_int) |>
    as_tibble() |>
    mutate(empi = empi)

  max_window_after_int <- max_window_after_int + 1
  
  periods_after <-
    seq.POSIXt(from = A_time, to = max_time, length.out = max_window_after_int) |>
    as_tibble() |>
    mutate(empi = empi)
  
  out <- full_join(periods_before, periods_after) |> select(empi, l_start = value)

  }

# go backwards to create intervals until t1_start
windows_before_intubation <-
  intubated_first |>
  head() |>
  mutate(max_time = pmax(Y_time, Cens_time, na.rm = T),
         max_windows = ceiling(as.numeric(difftime(max_time, t1_start, units = "days"))),
         max_window_int = ceiling(as.numeric(difftime(A_time, t1_start, units = "days"))),
         max_window_after_int = ceiling(as.numeric(difftime(max_time, A_time, units = "days"))),
         max_window_study_end = ceiling(as.numeric(difftime(max_time, t1_start, units = "days"))),
         window = 1) |> 
  select(empi, t1_start, A_time, max_window_int, max_time, max_window_after_int) |>
  pmap_dfr(create_periods)

windows_before_intubation |>
  group_by(empi) |>
  mutate(z_end = lead(l_start)) |>
  drop_na(z_end) |>
  mutate(z_end = case_when(row_number() == n() ~ z_end,
                           TRUE ~ z_end - seconds(1)))

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

# A function to define L_t and Z_t by covariate
# arguemnt covar is a string containing lab name, eg "ddimer"
# window_num indicates what t we're determining
# returns a list for every time point and L_t and Z_t values
never_labs <- function(covar, window_num){
  
  # extract labs within the window of interest
  rel_labs <-
    nevers_clean |>
    filter(window == window_num) |>
    left_join(all_labs_clean |> filter(name == covar)) |>
    filter(result_dt %within% interval(l_start, z_end))
  
  # if in L interval, record as such
  Ls <- rel_labs |>
    drop_na(result_value) |>
    filter(result_dt %within% interval(l_start, l_end)) |>
    arrange(result_dt) |>
    distinct(empi, .keep_all = T) |>
    mutate(missing = ifelse(is.na(result_value), 1, 0)) |>
    select(empi, window, L_value = result_value)
  
  # if in Z interval, document as such
  Zs <-  rel_labs |>
    drop_na(result_value) |>
    filter(result_dt %within% interval(z_start, z_end)) |>
    arrange(result_dt) |>
    distinct(empi, .keep_all = T) |>
    select(empi, window, Z_value = result_value)

  out <-
    nevers_clean |>
    select(empi, window) |>
    filter(window == window_num) |>
    mutate(covar = covar) |>
    left_join(Ls) |>
    left_join(Zs)
    
  return(out)

    }

# Create a grid of all covariates and time windows to map through
map_grid <- expand.grid("name" = unique(all_labs_clean$name), "window" = 1:28)
all_labs <- map2(map_grid$name, map_grid$window, ~never_labs(.x, .y))

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
              values_from = c(A,M))

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


