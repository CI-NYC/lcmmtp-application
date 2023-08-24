#######################################################
#######################################################
### Create relevant date-times for computing necessary variables L_t, A_t, Z_t, M_t, C_t, Y_t
### Kat Hoffman
#######################################################
#######################################################

library(tidyverse)

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
  mutate(t1_end = t1_start + hours(24),
         M_time = case_when(M_time > Y_time ~ Y_time - hours(1), # fix temporality of a few AKIs that are after death
                            M_time > Cens_time ~ NA_POSIXct_, # some patients had AKI detected thru labs collected after discharge, mark these AKI instances as NA because we don't care about it (and will mess up our coding logic later)
                            TRUE ~ M_time)) |>
  filter(A_time >= t1_start | is.na(A_time), # make sure all times are sensical (n=22, 0.6% observations removed for data entry error)
         M_time >= t1_start | is.na(M_time),
         Y_time >= t1_start | is.na(Y_time),
         Cens_time >= t1_start | is.na(Cens_time))

saveRDS(cp_dt, here::here("data/derived/dts_cohort.rds"))
