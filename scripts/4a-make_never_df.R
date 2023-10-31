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

first_o2 <- read_rds("data/derived/first_o2.rds")

### Filter for patients who were never intubated or AKI
nevers <-
  cp_dt |>
  filter((is.na(A_time) & is.na(M_time))) 

create_periods  <- function(empi, t1_start, max_time, max_windows){
  
  max_windows <- max_windows + 1 # need to make a start/stop sequence that is 1 longer than the desired number of time windows 
  
  # divide up study start time to time of aki into equal intervals
  periods <- 
    seq.POSIXt(from = t1_start, to = max_time, length.out = max_windows) |>
    as_tibble() |>
    mutate(empi = empi)
  
  # return the unique sequence of times; these are our start times for each L window
  out <- periods |> select(empi, l_start = value)
  
}

nevers_clean <-
  nevers |>
  mutate(
    max_time = pmax(Y_time, Cens_time, na.rm = T),
    max_windows = ceiling(as.numeric(difftime(max_time, t1_start, units = "days"))), # get max number of 24 hour windows,
    window = 1) |>
  select(empi, t1_start, max_time, max_windows) |>
  pmap_dfr(create_periods) |>
  group_by(empi) |>
  mutate(z_end = lead(l_start)) |>
  drop_na(z_end) |>
  mutate(l_start = case_when(row_number() == 1 ~ l_start,
                             TRUE ~ l_start + seconds(1)) # we want a new L to start right *after* the person gets AKI, and current Z to end when they are intubated
  ) |>
  mutate(l_end = as.POSIXct((as.numeric(l_start) + as.numeric(z_end)) / 2, origin = '1970-01-01'), # get the midway point between L and Z start/end times, call this l_end
         z_start = l_end + seconds(1)) |> # Z for that t starts a second after L ends
  select(empi, l_start, l_end, z_start, z_end) |>
  mutate(window = row_number()) |>
  ungroup()

# Create a grid of all covariates and time windows to map through
map_grid <- expand.grid("name" = unique(all_labs_clean$name), "window" = 1:28)
all_labs <- map2(map_grid$name, map_grid$window, ~divide_labs(nevers_clean, .x, .y))

# put together all L_t and Z_t in long format and create missingness indicators
Ls_and_Zs <- 
  data.table::rbindlist(all_labs) |> # bind all rows together
  mutate(covar = snakecase::to_snake_case(as.character(covar))) |> # clean up covariate names for wide col names later
  group_by(empi, covar) |> 
  arrange(empi, window) |>
  fill(L_value, Z_value, .direction = "down") |> # Last Observation Carried Forward, by empi and covariate
  mutate(L_missing = ifelse(is.na(L_value), 1, 0), # create indicators for missing L and Z
         Z_missing = ifelse(is.na(Z_value), 1, 0)) |>
  mutate(L_value = ifelse(is.na(L_value), -99999, L_value), # if no observations before, fill in with -99999 (could switch to median value)
         Z_value = ifelse(is.na(Z_value), -99999, Z_value)) 


# fill in M = 0 for all the mediators (this group never gets AKI)
Ms <- 
  nevers_clean |>
  mutate(M = 0)

# fill in A = 0 
As <- 
  nevers_clean |>
  left_join(first_o2 |> ungroup()) |>
  mutate(A = case_when(first_o2 < z_end ~ 1, TRUE ~ 0)) |>
  select(-first_o2, -any_o2)

# Note I'll probably need to carry these Y's through (deterministic)
Ys <-
  nevers_clean |>
  select(empi, window, l_start, z_end) |>
  left_join(nevers |> select(empi, Y_time, Cens_time)) |>
  # if died in this period, outcome Y is 1
  mutate(Y = case_when(Y_time %within% interval(l_start, z_end) ~ 1,
                       # Y_time > l_start ~ 1,
                       # if they were discharged in this period, no outcome observed (NA)
                       Cens_time %within% interval(l_start, z_end) ~ NA_real_,
                       # otherwise, outcome is 0
                       TRUE ~ 0))

# Note I'll need to fill in the rest of the censoring variables as 0 (or leave as NA?) once in long format
Cs <- nevers_clean |>
  select(empi, window, l_start, z_end) |>
  left_join(nevers |> select(empi, Y_time, Cens_time)) |>
  # if they were discharged in this period, no outcome observed
  mutate(Observed = case_when(Cens_time %within% interval(l_start, z_end) ~ 0,
                       # otherwise, observed
                       TRUE ~ 1))


# Save long format data ------------------------------------------------

fp <- "data/derived/nevers"
saveRDS(Ls_and_Zs, here::here(fp, "Ls_and_Zs.rds"))
saveRDS(Ms, here::here(fp, "Ms.rds"))
saveRDS(As, here::here(fp, "As.rds"))
saveRDS(Cs, here::here(fp, "Cs.rds"))
saveRDS(Ys, here::here(fp, "Ys.rds"))
saveRDS(nevers_clean, here::here(fp, "nevers.rds"))

