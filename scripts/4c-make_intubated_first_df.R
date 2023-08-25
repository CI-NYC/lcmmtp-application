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
int_first <- 
  cp_dt |>
  filter(A_time < M_time | (!is.na(A_time) & is.na(M_time)))

# Function to create periods based on intubation (A_time)
# Finds the time of intubation and creates intervals in approximately 12 hour intervals from study start date
# Then makes intervals to max time in approximately 12 hour intervals
# Number of  intervals is determined by the number of intervals there would be if the patient had had their time divided by true 12 hour intervals (like the never group)
create_periods  <- function(empi, t1_start, A_time, windows_12h_to_int_odd, max_time, windows_12h_to_end_odd){

  windows_12h_to_int_odd <- windows_12h_to_int_odd + 1 # need to make a start/stop sequence that is 1 longer than the desired number of time windows 
  
  # divide up study start time to time of intubation into equal intervals
  periods_before <- 
    seq.POSIXt(from = t1_start, to = A_time, length.out = windows_12h_to_int_odd) |>
    as_tibble() |>
    mutate(empi = empi)
  
  windows_12h_to_end_odd <- windows_12h_to_end_odd + 1 # need to make a start/stop sequence that is 1 longer than the desired number of time windows 
  
  # divide up time of intubation to end of study follow up into equal intervals
  periods_after <-
    seq.POSIXt(from = A_time, to = max_time, length.out = windows_12h_to_end_odd) |>
    as_tibble() |>
    mutate(empi = empi)

  # return the entire data frame of l and z start/end times for that patient
  out <- full_join(periods_before, periods_after) |>
    mutate(l_or_z = case_when(row_number() %% 2 == 0 ~ "z_start",
                              row_number() == n() ~ "z_end",
                              TRUE ~ "l_start"
    )) |>
    groupdata2::group(n=2, method="greedy", col_name = "window") |>
    ungroup() |>
    mutate(window = case_when(row_number() == n() ~ max(as.numeric(window)) - 1,
                              TRUE ~ as.numeric(window))) |>
    pivot_wider(id_cols = c(empi, window),
                names_from = "l_or_z",
                values_from = "value",
                values_fill = NA) |>
    unnest(c(l_start, z_start, z_end)) |>
    mutate(z_end = case_when(!is.na(z_end) ~ z_end, # last z_end stays as max_time
                             TRUE ~ lead(l_start) - seconds(1)), # otherwise one second before the next l_start
           l_end = z_start - seconds(1) # l ends one second before z starts
           ) |>
    select(empi, window, l_start, l_end, z_start, z_end)
    
  return(out) # data frame with l/z start and end for one patient
  
}

# create Z/L intervals using above create_periods function (with pmap) - takes 53 seconds
int_windows <-
  int_first |>
  mutate(max_windows = ceiling(as.numeric(difftime(max_time, t1_start, units = "days"))),
         # divide into nearest odd 12-hours periods between study start and intubation (and then intubation to max_time), so that A_t falls at the end of L_t and the start of Z_t
         windows_12h_to_int = round(as.numeric(difftime(A_time, t1_start, units = "hours"))/12),
         windows_12h_to_int_odd = case_when(windows_12h_to_int == 0 ~ 1, # if closest interval is 0, then make this 1 window
                                            windows_12h_to_int %% 2 == 0 ~ windows_12h_to_int - 1,
                                            TRUE ~ windows_12h_to_int),
         windows_12h_to_end = round(as.numeric(difftime(max_time, A_time, units = "hours"))/12),
         windows_12h_to_end_odd = case_when(windows_12h_to_end == 0 ~ 1, # if closest interval is 0, then make this 1 window
                                            windows_12h_to_end %% 2 == 0 ~ windows_12h_to_end - 1,
                                            TRUE ~ windows_12h_to_end),
  ) |>
  select(empi, t1_start, A_time, windows_12h_to_int_odd, max_time, windows_12h_to_end_odd) |>
  pmap_dfr(create_periods)

# temporary data frame containing empi, window, l_start, l_end, z_start, z_end where l's and z's are dividied into equal intervals
# need to fix in next step for the n=XX intervals in which patients get AKI during that time period (so that Z ends at time of AKI and a new L starts)
int_windows_tmp <- int_windows |>
  group_by(empi) |>
  mutate(z_end = lead(l_start)) |>
  drop_na(z_end) |>
  mutate(l_start = case_when(row_number() == 1 ~ l_start,
                             TRUE ~ l_start + seconds(1)) # we want a new L to start right *after* the person is intubated, and Z to end when they are intubated
  ) |>
  mutate(l_end = as.POSIXct((as.numeric(l_start) + as.numeric(z_end)) / 2, origin = '1970-01-01'), # get the midway point between L and Z start/end times, call this l_end
         z_start = l_end + seconds(1)) |> # Z for that t starts a second after L ends
  select(empi, l_start, l_end, z_start, z_end) |>
  mutate(window = row_number()) |>
  ungroup() |>
  mutate(index = row_number()) # to modify intervals in which intubation occurs in next lines

# filter to windows in which intubation occurs (at some point after AKI)
# modify z start and l end to reflect this A_time
aki_windows_mod <-
  int_windows_tmp |>
  left_join(int_first |> select(empi, A_time)) |>
  drop_na(A_time) |>
  filter(A_time %within% interval(l_start, z_end)) |>
  mutate(z_start = A_time,
         l_end = z_start - seconds(1)) |>
  select(-A_time)

# merge modified data rows into final aki windows dataset 
aki_windows_clean <-
  aki_windows_tmp |>
  filter(!(index %in% aki_windows_mod$index)) |> # filter rows that need to be modified out of the original data
  full_join(aki_windows_mod) |> # add modified rows (from above) into final clean data
  arrange(index)

aki_windows_clean

# Create a grid of all covariates and time windows to map through
map_grid <- expand.grid("name" = unique(all_labs_clean$name), "window" = 1:28)
all_labs <- map2(map_grid$name, map_grid$window, ~divide_labs(aki_windows_clean, .x, .y))

# put together all L_t and Z_t in long format and create missingness indicators
Ls_and_Zs <- 
  data.table::rbindlist(all_labs) |> # bind all rows together
  mutate(L_missing = ifelse(is.na(L_value), 1, 0), # create indicators for missing L and Z
         Z_missing = ifelse(is.na(Z_value), 1, 0)) |>
  group_by(empi, covar) |> 
  arrange(empi, window) |>
  fill(L_value, Z_value, .direction = "down") |> # Last Observation Carried Forward, by empi and covariate
  mutate(L_value = ifelse(is.na(L_value), -99999, L_value), # if no observations before, fill in with -99999 (could switch to median value)
         Z_value = ifelse(is.na(Z_value), -99999, Z_value))


# fill in M = 1 for all the windows where intubation occurs on or after
# need to update this for different levels of supp o2 later
Ms <- 
  aki_windows_tmp |>
  left_join(aki_first |> select(empi, M_time)) |>
  mutate(M_this_window = case_when(M_time == z_end ~ window)) |> # mark which window AKI occurs in (it's the end of the Z interval)
  fill(M_this_window, .direction = "downup") |>
  mutate(M = case_when(window >= M_this_window ~ 1, TRUE ~ 0)) |>
  select(empi, window, M_this_window, M)

# mark which window intubation occurred --> switch to intubation = 2 later
As <-
  aki_windows_clean |>
  left_join(aki_first |> select(empi, A_time)) |>
  left_join(Ms) |>
  mutate(A = case_when(is.na(A_time) ~ 0, # never intubated
                       l_start < A_time ~ 0,
                       l_start >= A_time ~ 1)) 

# Note I'll probably need to carry these Y's through (deterministic)
Ys <-
  aki_windows_clean |>
  select(empi, window, l_start, z_end) |>
  left_join(aki_first |> select(empi, Y_time, Cens_time)) |>
  # if died in this period, outcome Y is 1
  mutate(Y = case_when(Y_time == z_end ~ 1,
                       # if they were discharged in this period, no outcome observed (NA)
                       Cens_time == z_end ~ NA_real_,
                       # otherwise, outcome is 0
                       TRUE ~ 0))

# Note I'll need to fill in the rest of the censoring variables as 0 (or leave as NA?) once in long format
Cs <- aki_windows_clean |>
  select(empi, window, l_start, z_end) |>
  left_join(aki_first |> select(empi, Y_time, Cens_time)) |>
  # if they were discharged in this period, no outcome observed
  mutate(C = case_when(Cens_time == z_end ~ 0,
                       # otherwise, observed
                       TRUE ~ 1))

# Turn into wide form data ------------------------------------------------

max_window_data <- 28

M_wide <- 
  Ms |>
  filter(window < max_window_data) |>
  pivot_wider(id_cols=empi,
              names_from = window,
              values_from = M)

A_wide <- 
  As |>
  filter(window < max_window_data) |>
  pivot_wider(id_cols=empi,
              names_from = window,
              values_from = A)

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

aki_first_wide <-
  reduce(list(aki_first, Ms, As, Ls_and_Zs_wide, Ys_wide, Cs_wide),
         ~left_join(.x, .y))

saveRDS(aki_first_wide, here::here("data/derived/aki_first_wide.rds"))


