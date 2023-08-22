#### Crude code to test out one time point for temporality logic

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
                            TRUE ~ M_time)) |>
  filter(!(A_time < t1_start) |
         !(M_time < t1_start) |
         !(Y_time < t1_start) |
         !(Cens_time < t1_start)) # remove 3 rows with nonsensical times

saveRDS(cp_dt, here::here("data/derived/dts_cohort.rds"))
