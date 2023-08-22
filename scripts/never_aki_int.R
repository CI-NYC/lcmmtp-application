#######################################################
########### Wrangle data for NEVERs
########### Kat Hoffman
#######################################################

# setup -------------------------------------------------------------------

library(tidyverse)
# library(tidylog)
library(janitor)

cp_dt <- read_rds(here::here("data/derived/dts_cohort.rds"))
labs <- read_rds(here::here("data", "raw", "2020-08-17", "wcm_labs_filtered.rds"))

### TO SWITCH OUT:
# Make this one long cleaned lab list

### All d-dimers, all pts
ddimers <- labs |>
  as_tibble() |>
  filter(result_name == "D-Dimer") |> # all d-dimers have this
  mutate(result_dt = as_datetime(result_time),
         result_value = readr::parse_number(ord_value)) |> 
  select(empi, result_dt, result_value)|>
  mutate(name = "ddimer")

### All potassium, all pts
potassium <- labs |>
  as_tibble() |>
  filter(result_name == "Potassium Level") |>
  mutate(result_dt = as_datetime(result_time),
         result_value = readr::parse_number(ord_value)) |> 
  select(empi, result_dt, result_value) |>
  mutate(name = "potassium")

all_labs_clean <- bind_rows(ddimers, potassium)

### Filter for patients who were never intubated or AKI
nevers <-
  cp_dt |>
  filter((is.na(A_time) & is.na(M_time)) 
  ) 

# 
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
      select(window, l_start, l_end, z_start, z_end)
  }) 
  ) |>
  unnest(cols = data)

nevers_clean |>
  filter(window == 1) |>
  left_join(ddimers) |>
  filter(result_dt %within% interval(l_start, z_end))

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

never_labs("ddimer", 1)

map_grid <- expand.grid("name" = c("ddimer","potassium"), "window" = 1:28)
map_grid
all_labs <- map2(map_grid$name, map_grid$window, ~never_labs(.x, .y))

Ls_and_Zs <- 
  data.table::rbindlist(all_labs) |>
  mutate(L_missing = ifelse(is.na(L_value), 1, 0),
         Z_missing = ifelse(is.na(Z_value), 1, 0)) |>
  mutate(L_value = ifelse(is.na(L_value), -99999, L_value),
         Z_value = ifelse(is.na(Z_value), -99999, Z_value))


# fill in M = 0 for all the mediators (this group never gets AKI)
# fill in A = 0 (change for oxygenation data)
M_and_A <- 
  nevers_clean |>
  mutate(M = 0,
         A = 0
         )

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
