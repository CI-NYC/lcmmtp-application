#######################################################
########### Wrangle data for NEVERs
########### Kat Hoffman
#######################################################

# setup -------------------------------------------------------------------

library(tidyverse)
# library(tidylog)
library(janitor)

source("scripts/0-functions.R")

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
  pmap_dfr(create_periods) |>
  arrange(empi, window) |>
  mutate(index = row_number()) |>  # index for modifying in next step
  left_join(int_first |> select(empi, M_time))  # add in aki time for next step

int_windows |> distinct(empi)

#### THREE TYPES OF PATIENTS
# 1. get intubated and never have AKI [not in check_empis]
# 2. get intubated then get AKI before the death date, but not their final time point [in check_empis] -- need to change z_end and l_start for two rows
# 3. get intubated then get AKI before the death date, in their final time point [in check_empis] -- need to add a row of time periods
# 3. get intubated then get AKI right before death [in check_empis] -- need to 

# 1. get intubated and never have AKI
int_windows_section1 <-
  int_windows |>
  filter(is.na(M_time))

# these empis AKI after intubation -- need to fix some
int_windows_section23 <-
  int_windows |>
  filter(!is.na(M_time)) |>
  group_by(empi) |>
  arrange(empi, l_start) |>
  filter(M_time %within% interval(l_start, z_end) | M_time %within% interval(lag(l_start), lag(z_end))) |>
  group_by(empi) |>
  arrange(l_start) |>
  add_count()

int_windows_section23 |> filter(n == 3)

nrow(int_windows)
nrow(int_windows_section23)
nrow(int_windows_section1)

distinct(int_windows, empi) |> nrow()
distinct(int_windows_section23, empi) |> nrow()
distinct(int_windows_section1, empi) |> nrow()

# Keep the row that they got AKI *and* the row after, because we need to change z_end and l_start
int_windows_section2 <-
  int_windows_section23 |>
  filter(n == 2) # remove patients that get AKI in final time point for now

# change 1st row (new z_end time)
int_windows_section2_mod1 <- 
  int_windows_section2 |>
  filter(row_number() == 1) |>
  mutate(z_end = M_time ) |> #,  # new end time, No new midpoints, otherwise it'll mess up intubation timing
                # l_end = as.POSIXct((as.numeric(l_start) + as.numeric(z_end)) / 2, origin = '1970-01-01'), # get the midway point between L and Z start/end times, call this l_end
                #  z_start = l_end + seconds(1))  |>
  select(empi, window, l_start, l_end, z_start, z_end, index, M_time)

int_windows_section2_mod1 |>
  filter(empi == "1000017735")

# change 2nd row (new l_start time)
int_windows_section2_mod2 <- 
  int_windows_section2 |>
  filter(row_number() == 2) |>
  mutate(l_start = M_time + seconds(1 )) |> #,  # new start time, No new midpoints, otherwise it'll mess up intubation timing
       #  l_end = as.POSIXct((as.numeric(l_start) + as.numeric(z_end)) / 2, origin = '1970-01-01'), # get the midway point between L and Z start/end times, call this l_end
       #  z_start = l_end + seconds(1))  |>
  select(empi, window, l_start, l_end, z_start, z_end, index, M_time)


# Keep the patients that got AKI in their final time point
int_windows_section3 <-
  int_windows_section23 |>
  filter(n == 1) # remove patients that get AKI in final time point for now

dup <- int_windows_section3


# modify existing final row so that M is the final time point in this window t
section3_mod_row <-
  int_windows_section3 |>
  mutate(z_end = M_time ) |> #,  # new end time, No new midpoints, otherwise it'll mess up intubation timing
         #l_end = as.POSIXct((as.numeric(l_start) + as.numeric(z_end)) / 2, origin = '1970-01-01'), # get the midway point between L and Z start/end times, call this l_end
         #z_start = l_end + seconds(1))  |>
  select(empi, window, l_start, l_end, z_start, z_end, index, M_time)

# create a new final row for these patients (so that there is an additional time window to measure labs after M_t)
section3_new_row <-
  dup |>
  mutate(l_start = M_time + seconds(1)) |> #,# new start time, No new midpoints, otherwise it'll mess up intubation timing
        # l_end = as.POSIXct((as.numeric(l_start) + as.numeric(z_end)) / 2, origin = '1970-01-01'), # get the midway point between L and Z start/end times, call this l_end
        # z_start = l_end + seconds(1))  |>
  select(empi, window, l_start, l_end, z_start, z_end, index, M_time)

# merge into final modified data set --------------------------------------

# merge old and new data for section 2
section2 <- int_windows |>
  filter(empi %in% int_windows_section2$empi) |> 
  filter(!(index %in% int_windows_section2$index)) |>
  bind_rows(int_windows_section2_mod1) |>
  bind_rows(int_windows_section2_mod2) |>
  arrange(index, l_start)

# merge old and new data for section 3
section3 <- int_windows |>
  filter(empi %in% int_windows_section3$empi) |> 
  filter(!(index %in% int_windows_section3$index)) |>
  bind_rows(section3_mod_row) |>
  bind_rows(section3_new_row) |>
  arrange(index, l_start)

# merge into final modified data set --------------------------------------

int_windows_clean <-
  int_windows_section1 |>
  full_join(section2) |>
  full_join(section3) |>
  group_by(empi) |>
  arrange(empi, l_start) |>
  mutate(window = row_number())

# # check if there is someone in ids but not section 2 or 3

# nrow(int_windows)
# nrow(int_windows_clean)
# 
# int_windows_section1 |> left_join(int_first |> select(empi, M_time)) |> filter(is.na(M_time)) |> nrow()
# int_windows |>
#   left_join(int_first |> select(empi, M_time))  |>
#   filter(!is.na(M_time)) |> nrow()


# Create a grid of all covariates and time windows to map through
map_grid <- expand.grid("name" = unique(all_labs_clean$name), "window" = 1:28)
all_labs <- map2(map_grid$name, map_grid$window, ~divide_labs(int_windows_clean, .x, .y))

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
  int_windows_clean |>
  left_join(int_first |> select(empi, M_time)) |>
  mutate(M_this_window = case_when(M_time == z_end ~ window)) |> # mark which window AKI occurs in (it's the end of the Z interval)
  fill(M_this_window, .direction = "downup") |>
  mutate(M = case_when(window >= M_this_window ~ 1, TRUE ~ 0)) |>
  select(empi, window, M_this_window, M)

# mark which window intubation occurred --> switch to intubation = 2 later
As <-
  int_windows_clean |>
  left_join(int_first |> select(empi, A_time)) |>
  mutate(A = case_when(is.na(A_time) ~ 0, # never intubated
                       l_start < A_time ~ 0,
                       l_start >= A_time ~ 1))

# Note I'll probably need to carry these Y's through (deterministic)
Ys <-
  int_windows_clean |>
  select(empi, window, l_start, z_end) |>
  left_join(int_first |> select(empi, Y_time, Cens_time)) |>
  # if died in this period, outcome Y is 1
  mutate(Y = case_when(Y_time == z_end ~ 1,
                       # if they were discharged in this period, no outcome observed (NA)
                       Cens_time == z_end ~ NA_real_,
                       # otherwise, outcome is 0
                       TRUE ~ 0))

# Note I'll need to fill in the rest of the censoring variables as 0 (or leave as NA?) once in long format
Cs <- int_windows_clean |>
  select(empi, window, l_start, z_end) |>
  left_join(int_first |> select(empi, Y_time, Cens_time)) |>
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
              values_from = M,
              names_prefix = "M_")

A_wide <- 
  As |>
  filter(window < max_window_data) |>
  pivot_wider(id_cols=empi,
              names_from = window,
              values_from = A,
              names_prefix = "A_")

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

int_first_wide <-
  reduce(list(int_first, M_wide, A_wide, Ls_and_Zs_wide, Ys_wide, Cs_wide),
         ~left_join(.x, .y))

saveRDS(int_first_wide, here::here("data/derived/int_first_wide.rds"))


