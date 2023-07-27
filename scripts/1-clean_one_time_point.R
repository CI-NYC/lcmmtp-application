#### Crude code to test out one time point for temporality logic

library(tidyverse)

cohort <- read_rds(here::here("data","derived","hospitalized_cohort_with_queens.rds"))
labs <- read_rds(here::here("data", "raw", "2020-08-17", "wcm_labs_filtered.rds"))

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

### All d-dimers, all pts
ddimers <- labs |>
  as_tibble() |>
  filter(result_name == "D-Dimer") |> # all d-dimers have this
  mutate(result_dt = as_datetime(result_time),
         result_value = readr::parse_number(ord_value)) |> 
  select(empi, result_dt, result_value)

### All potassium, all pts
potassium <- labs |>
  as_tibble() |>
  filter(result_name == "Potassium Level") |>
  mutate(result_dt = as_datetime(result_time),
         result_value = readr::parse_number(ord_value)) |> 
  select(empi, result_dt, result_value)


### Find patients who were never intubated or AKI
nevers <-
  cp_dt |>
  filter((is.na(A_time) & is.na(M_time)) 
           ) |>
  mutate(z1_start = t1_start + hours(12)) # define L1 as first 12 hours, Z1 as second 12 hours

nevers_l1_ddimer <-
  nevers |>
  left_join(ddimers) |>
  filter(result_dt %within% interval(t1_start, z1_start)) |>
  arrange(desc(result_dt)) |>
  select(empi, L1_ddimer = result_value) |>
  distinct(empi, .keep_all = T)

nevers_z1_ddimer <-
  nevers |>
  left_join(ddimers) |>
  filter(result_dt %within% interval(z1_start, t1_end)) |>
  arrange(desc(result_dt)) |>
  select(empi, "Z1_ddimer" = result_value) |>
  distinct(empi, .keep_all = T) 

nevers_l1_potassium <-
  nevers |>
  left_join(potassium) |>
  filter(result_dt %within% interval(t1_start, z1_start)) |>
  arrange(desc(result_dt)) |>
  select(empi, L1_potassium = result_value) |>
  distinct(empi, .keep_all = T) 

nevers_z1_potassium <-
  nevers |>
  left_join(potassium) |>
  filter(result_dt %within% interval(z1_start, t1_end)) |>
  arrange(desc(result_dt)) |>
  select(empi, Z1_potassium = result_value) |>
  distinct(empi, .keep_all = T) 

nevers_full_int <- reduce(list(nevers_l1_potassium,
            nevers_l1_ddimer,
            nevers_z1_potassium,
            nevers_z1_ddimer), ~full_join(.x,.y)) |>
  full_join(nevers) |>
  select(empi, t1_start, t1_end, A_time, M_time, Y_time, Cens_time,
         contains("l1"), contains("z1")
         )

nevers_final <- nevers_full_int |>
  mutate(Y1 = case_when(Y_time < t1_start ~ 1, # (for when it's not the first time point) need to check if already died
                        Y_time %within% interval(t1_start, t1_end) ~ 1, # anyone who dies in this interval gets a 1
                        Cens_time %within% interval(t1_start, t1_end) ~ NA_real_, # anyone who is censored (discharged) in this interval has an unknown outcome 
                        TRUE ~ 0
                        ),
         C1 = case_when(Cens_time < t1_start ~ NA_real_, # for when it's not the first time point, need to check if already censored
                        Cens_time %within% interval(t1_start, t1_end) ~ 0,
                        TRUE ~ 1
                        ),
         A1 = 0, # this will be changed depending on oxygenation
         M1 = 0) |>
  select(empi, starts_with("L1"), A1, M1, starts_with("Z1"), C1, Y1, -z1_start)
  # mutate(C1 = ifelse(Cens_time < t1_end, 0, 1))

nevers_final
# 
# # intubated and aki within same 24 hour interval
# im_t1 <- cp_dt |>
#   filter(A_time %within% interval(t1_start, t1_end),
#          M_time %within% interval(t1_start, t1_end)) 
# 
# int_im_t1 <-
#   im_t1 |>
#   filter(A_time < M_time) |>
#   mutate(A1 = 1, M1 = 1)
# 
# int_im_t1_l1_ddimer <-
#   int_im_t1 |>
#   left_join(ddimers) |>
#   filter(result_dt %within% interval(t1_start, A_time)) |>
#   arrange(desc(result_dt)) |>
#   select(empi, L1_ddimer = result_value) |>
#   distinct(empi, .keep_all = T)
# 
# int_im_t1_z1_ddimer <-
#   int_im_t1 |>
#   left_join(ddimers) |>
#   filter(result_dt %within% interval(A_time, M_time)) |>
#   arrange(desc(result_dt)) |>
#   select(empi, Z1_ddimer = result_value) |>
#   distinct(empi, .keep_all = T)
# 
# int_im_t1_l1_potassium <-
#   int_im_t1 |>
#   left_join(potassium) |>
#   filter(result_dt %within% interval(t1_start, A_time)) |>
#   arrange(desc(result_dt)) |>
#   select(empi, L1_potassium = result_value) |>
#   distinct(empi, .keep_all = T)
# 
# int_im_t1_z1_potassium <-
#   int_im_t1 |>
#   left_join(potassium) |>
#   filter(result_dt %within% interval(A_time, M_time)) |>
#   arrange(desc(result_dt)) |>
#   select(empi, Z1_potassium = result_value) |>
#   distinct(empi, .keep_all = T)
# 
# # join all intubated first, intubated and aki in same time period, measures together
# int_im_t1_full_int <- reduce(list(int_im_t1_l1_potassium,
#                                int_im_t1_l1_ddimer,
#                                int_im_t1_z1_potassium,
#                                int_im_t1_z1_ddimer), ~full_join(.x,.y)) |>
#   full_join(int_im_t1) |>
#   select(empi, t1_start, t1_end, A_time, M_time, Y_time, Cens_time,
#          contains("l1"), contains("z1")
#   )
# 
# int_im_t1


### Find patients for whom intubated before AKI (A_t < M_t)
intubated_first <- 
  cp_dt |>
  filter(A_time < M_time | (!is.na(A_time) & is.na(M_time))) |>
  filter(A_time >= t1_start)

### Find patients for whom intubated before AKI (A_t < M_t)
aki_first <- 
  cp_dt |>
  filter(M_time < A_time | (!is.na(M_time) & is.na(A_time))) |>
  filter(M_time >= t1_start)

int_times <- intubated_first |>
  mutate(hours_to_intubation = (t1_start  %--%  A_time) %/% hours(1),
         hours_int_to_aki = (A_time  %--%  M_time) %/% hours(1)) 

max_hours_to_aki <- 12

# patients who are intubated at time "1" (less than 12 hours of admission)
int_immediate <-
  int_times |>
  filter(hours_to_intubation < 12) |>
  mutate(A1 = 2, # these patients were intubated at time 1
         M1 = case_when(hours_int_to_aki < max_hours_to_aki ~ 1,
                        TRUE ~ 0), # these patients were not yet AKI
         t1_end = case_when(hours_int_to_aki < max_hours_to_aki ~ M_time, # make the end of time for Z variables to occur within 12 hours of intubation, or as soon as AKI occurs
                            TRUE ~ A_time + hours(max_hours_to_aki)),
         z1_start = A_time)

int_immediate_l1_ddimer <-
  int_immediate |>
  left_join(ddimers) |>
  filter(result_dt %within% interval(t1_start, z1_start)) |>
  arrange(desc(result_dt)) |>
  select(empi, "L1_ddimer" = result_value) |>
  distinct(empi, .keep_all = T)

int_immediate_z1_ddimer <-
  int_immediate |>
  left_join(ddimers) |>
  filter(result_dt %within% interval(z1_start, t1_end)) |>
  arrange(desc(result_dt)) |>
  select(empi, "Z1_ddimer" = result_value) |>
  distinct(empi, .keep_all = T) 

  
int_immediate_l1_potassium <-
  int_immediate |>
  left_join(potassium) |>
  filter(result_dt %within% interval(t1_start, z1_start)) |>
  arrange(desc(result_dt)) |>
  select(empi, "L1_potassium" = result_value) |>
  distinct(empi, .keep_all = T)

int_immediate_z1_potassium <-
  int_immediate |>
  left_join(potassium) |>
  filter(result_dt %within% interval(z1_start, t1_end)) |>
  arrange(desc(result_dt)) |>
  select(empi, "Z1_potassium" = result_value) |>
  distinct(empi, .keep_all = T) 


int_immediate_full_int <- reduce(list(int_immediate_l1_potassium,
                                      int_immediate_l1_ddimer,
                                      int_immediate_z1_potassium,
                                      int_immediate_z1_ddimer), ~full_join(.x,.y)) |>
  full_join(int_immediate) |>
  select(empi, t1_start, t1_end, A_time, M_time, Y_time, Cens_time,
         contains("l1"), contains("z1"), contains("a1"), contains("m1"), -z1_start)

int_immediate_final <- int_immediate_full_int |>
  mutate(Y1 = case_when(Y_time < t1_start ~ 1, # (for when it's not the first time point) need to check if already died
                        Y_time %within% interval(t1_start, t1_end) ~ 1, # anyone who dies in this interval gets a 1 -- seems impossible if t1_end = M_time tho
                        Cens_time %within% interval(t1_start, t1_end) ~ NA_real_, # anyone who is censored (discharged) in this interval has an unknown outcome 
                        TRUE ~ 0
  ),
  C1 = case_when(Cens_time < t1_start ~ NA_real_, # for when it's not the first time point, need to check if already censored
                 Cens_time %within% interval(t1_start, t1_end) ~ 0,
                 TRUE ~ 1
  )) |>
  select(empi, starts_with("L1"), A1, M1, starts_with("Z1"), C1, Y1)
# mutate(C1 = ifelse(Cens_time < t1_end, 0, 1))

int_immediate_final


# patients who are intubated not at time 1, so less than 12 hours of admission
int_not_immediate <-
  int_times |>
  filter(hours_to_intubation >= 12) |>
  mutate(A1 = 0, # these patients were not intubated at time 1
         M1 = 0, # these patients were not yet AKI (intubated first group)
         z1_start = t1_start + hours(12))

int_not_immediate_l1_ddimer <-
  int_not_immediate |>
  left_join(ddimers) |>
  filter(result_dt %within% interval(t1_start, A_time)) |>
  arrange(desc(result_dt)) |>
  select(empi, "L1_ddimer" = result_value) |>
  distinct(empi, .keep_all = T)

int_not_immediate_z1_ddimer <-
  int_not_immediate |>
  left_join(ddimers) |>
  filter(result_dt %within% interval(A_time, t1_end)) |>
  arrange(desc(result_dt)) |>
  select(empi, "Z1_ddimer" = result_value) |>
  distinct(empi, .keep_all = T) 


int_not_immediate_l1_potassium <-
  int_not_immediate |>
  left_join(potassium) |>
  filter(result_dt %within% interval(t1_start, A_time)) |>
  arrange(desc(result_dt)) |>
  select(empi, "L1_potassium" = result_value) |>
  distinct(empi, .keep_all = T)

int_not_immediate_z1_potassium <-
  int_not_immediate |>
  left_join(potassium) |>
  filter(result_dt %within% interval(A_time, t1_end)) |>
  arrange(desc(result_dt)) |>
  select(empi, "Z1_potassium" = result_value) |>
  distinct(empi, .keep_all = T) 


int_not_immediate_full_int <- reduce(list(int_not_immediate_l1_potassium,
                                      int_not_immediate_l1_ddimer,
                                      int_not_immediate_z1_potassium,
                                      int_not_immediate_z1_ddimer), ~full_join(.x,.y)) |>
  full_join(int_not_immediate) |>
  select(empi, t1_start, t1_end, A_time, M_time, Y_time, Cens_time,
         contains("l1"), contains("z1"), contains("a1"), contains("m1"), -z1_start
  )

int_not_immediate_final <- int_not_immediate_full_int |>
  mutate(Y1 = case_when(Y_time < t1_start ~ 1, # (for when it's not the first time point) need to check if already died
                        Y_time %within% interval(t1_start, t1_end) ~ 1, # anyone who dies in this interval gets a 1 -- seems impossible if t1_end = M_time tho
                        Cens_time %within% interval(t1_start, t1_end) ~ NA_real_, # anyone who is censored (discharged) in this interval has an unknown outcome 
                        TRUE ~ 0
  ),
  C1 = case_when(Cens_time < t1_start ~ NA_real_, # for when it's not the first time point, need to check if already censored
                 Cens_time %within% interval(t1_start, t1_end) ~ 0,
                 TRUE ~ 1
  )) |>
  select(empi, starts_with("L1"), A1, M1, starts_with("Z1"), C1, Y1)
# mutate(C1 = ifelse(Cens_time < t1_end, 0, 1))

int_not_immediate_final

##### wrangle data for those immediately meeting aki criteria in the first 24 hours
## aki needs to be the end of the time period because the patient can't get intubated after, but marked as "before"

aki_immediate <- 
  aki_first |>
  mutate(hours_to_aki = (t1_start  %--%  M_time) %/% hours(1)) |> 
  filter(hours_to_aki < 24) |>
  mutate(A1 = 0, # these patients were not intubated before AKI
         M1 = 1, # these patients do have AKI at time 1
         t1_end = M_time, # end the time frame at the time of AKI
         z1_start = t1_start + hours_to_aki / 2) # z1 is halfway between the time from the last time interval
  

aki_immediate_l1_ddimer <-
  aki_immediate |>
  left_join(ddimers) |>
  filter(result_dt %within% interval(t1_start, z1_start)) |>
  arrange(desc(result_dt)) |>
  select(empi, "L1_ddimer" = result_value) |>
  distinct(empi, .keep_all = T)

aki_immediate_z1_ddimer <-
  aki_immediate |>
  left_join(ddimers) |>
  filter(result_dt %within% interval(z1_start, t1_end)) |>
  arrange(desc(result_dt)) |>
  select(empi, "Z1_ddimer" = result_value) |>
  distinct(empi, .keep_all = T) 


aki_immediate_l1_potassium <-
  aki_immediate |>
  left_join(potassium) |>
  filter(result_dt %within% interval(t1_start, z1_start)) |>
  arrange(desc(result_dt)) |>
  select(empi, "L1_potassium" = result_value) |>
  distinct(empi, .keep_all = T)

aki_immediate_z1_potassium <-
  aki_immediate |>
  left_join(potassium) |>
  filter(result_dt %within% interval(z1_start, t1_end)) |>
  arrange(desc(result_dt)) |>
  select(empi, "Z1_potassium" = result_value) |>
  distinct(empi, .keep_all = T) 


aki_immediate_full_int <- reduce(list(aki_immediate_l1_potassium,
                                      aki_immediate_l1_ddimer,
                                      aki_immediate_z1_potassium,
                                      aki_immediate_z1_ddimer), ~full_join(.x,.y)) |>
  full_join(aki_immediate) |>
  select(empi, t1_start, t1_end, A_time, M_time, Y_time, Cens_time,
         contains("l1"), contains("z1"), contains("a1"), contains("m1")
  )

aki_immediate_final <- aki_immediate_full_int |>
  mutate(Y1 = case_when(Y_time < t1_start ~ 1, # (for when it's not the first time point) need to check if already died
                        Y_time %within% interval(t1_start, t1_end) ~ 1, # anyone who dies in this interval gets a 1 -- seems impossible if t1_end = M_time tho
                        Cens_time %within% interval(t1_start, t1_end) ~ NA_real_, # anyone who is censored (discharged) in this interval has an unknown outcome 
                        TRUE ~ 0
  ),
  C1 = case_when(Cens_time < t1_start ~ NA_real_, # for when it's not the first time point, need to check if already censored
                 Cens_time %within% interval(t1_start, t1_end) ~ 0,
                 TRUE ~ 1
  )) |>
  select(empi, starts_with("L1"), A1, M1, starts_with("Z1"), C1, Y1, -z1_start)
# mutate(C1 = ifelse(Cens_time < t1_end, 0, 1))

aki_immediate_final |> count(Y1)






##### PATIENT WHO GET AKI first, but not in first 24 hours
aki_not_immediate <- 
  aki_first |>
  mutate(hours_to_aki = (t1_start  %--%  M_time) %/% hours(1)) |> 
  filter(hours_to_aki >= 24) |>
  mutate(A1 = 0, # these patients were not intubated before AKI
         M1 = 0, # these patients do have AKI at time 1
         z1_start = t1_start + hours(12)) # z1 is halfway between the time from the last time interval


aki_not_immediate_l1_ddimer <-
  aki_not_immediate |>
  left_join(ddimers) |>
  filter(result_dt %within% interval(t1_start, z1_start)) |>
  arrange(desc(result_dt)) |>
  select(empi, "L1_ddimer" = result_value) |>
  distinct(empi, .keep_all = T)

aki_not_immediate_z1_ddimer <-
  aki_not_immediate |>
  left_join(ddimers) |>
  filter(result_dt %within% interval(z1_start, t1_end)) |>
  arrange(desc(result_dt)) |>
  select(empi, "Z1_ddimer" = result_value) |>
  distinct(empi, .keep_all = T) 


aki_not_immediate_l1_potassium <-
  aki_not_immediate |>
  left_join(potassium) |>
  filter(result_dt %within% interval(t1_start, z1_start)) |>
  arrange(desc(result_dt)) |>
  select(empi, "L1_potassium" = result_value) |>
  distinct(empi, .keep_all = T)

aki_not_immediate_z1_potassium <-
  aki_not_immediate |>
  left_join(potassium) |>
  filter(result_dt %within% interval(z1_start, t1_end)) |>
  arrange(desc(result_dt)) |>
  select(empi, "Z1_potassium" = result_value) |>
  distinct(empi, .keep_all = T) 


aki_not_immediate_full_int <- reduce(list(aki_not_immediate_l1_potassium,
                                      aki_not_immediate_l1_ddimer,
                                      aki_not_immediate_z1_potassium,
                                      aki_not_immediate_z1_ddimer), ~full_join(.x,.y)) |>
  full_join(aki_not_immediate) |>
  select(empi, t1_start, t1_end, A_time, M_time, Y_time, Cens_time,
         contains("l1"), contains("z1"), contains("a1"), contains("m1")
  )

aki_not_immediate_final <- aki_not_immediate_full_int |>
  mutate(Y1 = case_when(Y_time < t1_start ~ 1, # (for when it's not the first time point) need to check if already died
                        Y_time %within% interval(t1_start, t1_end) ~ 1, # anyone who dies in this interval gets a 1 -- seems impossible if t1_end = M_time tho
                        Cens_time %within% interval(t1_start, t1_end) ~ NA_real_, # anyone who is censored (discharged) in this interval has an unknown outcome 
                        TRUE ~ 0
  ),
  C1 = case_when(Cens_time < t1_start ~ NA_real_, # for when it's not the first time point, need to check if already censored
                 Cens_time %within% interval(t1_start, t1_end) ~ 0,
                 TRUE ~ 1
  )) |>
  select(empi, starts_with("L1"), A1, M1, starts_with("Z1"), C1, Y1, -z1_start)
# mutate(C1 = ifelse(Cens_time < t1_end, 0, 1))


t1_covars <-
  reduce(list(aki_not_immediate_final,
  aki_immediate_final,
  int_not_immediate_final,
  int_immediate_final,
  nevers_final), ~full_join(.x,.y)) |>
  mutate(miss_L1_ddimer = ifelse(is.na(L1_ddimer), 1, 0),
         miss_L1_potassium = ifelse(is.na(L1_potassium), 1, 0),
         miss_Z1_ddimer = ifelse(is.na(Z1_ddimer), 1, 0),
         miss_Z1_potassium = ifelse(is.na(Z1_potassium), 1, 0),
         across(starts_with("L1"), ~ifelse(is.na(.x), mean(.x, na.rm=T), .x)),
         across(starts_with("Z1"), ~ifelse(is.na(.x), mean(.x, na.rm=T), .x)),
         ##### FIX THIS LATER
         A1 = case_when(A1 == 2 ~ 1, TRUE ~ A1)
         )

saveRDS(t1_covars, "data/derived/draft_t1_covars.rds")
