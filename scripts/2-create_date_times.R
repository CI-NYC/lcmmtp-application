#######################################################
#######################################################
### Create relevant date-times for computing necessary variables L_t, A_t, Z_t, M_t, C_t, Y_t
### Kat Hoffman
#######################################################
#######################################################

library(tidyverse)
library(tidylog)

cohort <- read_rds(here::here("data","derived","hospitalized_cohort_with_queens.rds"))

# select cohort date variables (cohort parsed date-times, cp_dt)
cp_dt <-
  cohort |>
  mutate(Cens_time = case_when(is.na(death_dt) ~ end_dt)) |>
  select(empi,
         t1_start = ed_adm_dt,
         A_time = intubation1_dt, # this should be changed to any respiratory support time
         M_time = aki_time,
         Y_time = death_dt,
         Cens_time
  ) |>
  mutate(max_time = pmax(Y_time, Cens_time, na.rm = T),
         M_time = case_when(M_time > Y_time ~ Y_time - hours(1), # fix temporality of a few AKIs that are after death
                            M_time > Cens_time ~ NA_POSIXct_, # some patients had AKI detected thru labs collected after discharge, mark these AKI instances as NA because we don't care about it (and will mess up our coding logic later)
                            TRUE ~ M_time),
         A_time = case_when(A_time > Cens_time ~ NA_POSIXct_, # some patients were intubated after discharge (readmitted), mark these as NA because we don't care about it (and will mess up our coding logic later)
                            (A_time > Y_time) & (A_time %within% interval(Y_time, Y_time + days(1))) ~ Y_time - hours(1), # assume patients with intubation date-time marked shortly after death date died during crash intubation and had data entry error
                            (A_time > Y_time) & (as.Date(A_time) == as.Date("2020-04-24 16:57:00")) ~ as.POSIXct("2020-03-24 16:57:00"), # month entry error
                            TRUE ~ A_time)
         ) |>
  filter(A_time >= t1_start | is.na(A_time), # make sure all times are sensical (n=23, 0.6% observations removed for data entry error)
         M_time >= t1_start | is.na(M_time),
         Y_time >= t1_start | is.na(Y_time),
         Y_time >= A_time | is.na(Y_time) | is.na(A_time), # 1 person with non-sensical time
         Cens_time > t1_start | is.na(Cens_time))

saveRDS(cp_dt, here::here("data/derived/dts_cohort.rds"))
