# clean baseline demographics

library(tidyverse)

cohort <- read_rds(here::here(fp, "hospitalized_cohort_with_queens.rds")) # contains baseline data
all_wide_pre_baseline <- read_rds(here::here(fp, "all_wide_pre_baseline.rds"))

### prep baseline data to merge
baseline <-
  cohort |>
  select(empi, red_cap_source,
         age, sex, race, ethnicity, bmi, smoking, 
         cad, home_o2_yn, dm, htn, cva, cirrhosis, ckd_or_esrd, asthma, copd, active_cancer, 
         immunosuppressed, ild, hiv, hypoxia_ed, hypoxia_ed_method) |>
  mutate(bmi_miss = case_when(is.na(bmi) ~ 1, TRUE ~ 0),
  home_o2_miss = case_when(is.na(home_o2_yn) ~ 1, TRUE ~ 0), 
  hypoxia_ed = case_when(hypoxia_ed == "No" ~ 0, hypoxia_ed == "Yes" ~ 1),
  race_miss = case_when(is.na(race) ~ 1, TRUE ~ 0),
  ethnicity_miss = case_when(is.na(ethnicity) ~ 1, TRUE ~ 0),
  race = case_when(is.na(race) ~ "Missing", TRUE ~ race),
  ethnicity = case_when(is.na(ethnicity) ~ "Missing", TRUE ~ ethnicity),
  smoking = case_when(is.na(smoking) ~ "No", TRUE ~ smoking),
  male = case_when(sex=="Male"~1, TRUE ~0),
  red_cap_source = case_when(red_cap_source == "QUEENS" ~ 0,
                             red_cap_source == "EAST" ~ 1),
  across(cad:hiv, ~case_when(.x == "No" ~ 0, TRUE ~ 1))) |>
  fastDummies::dummy_columns(select_columns = c("race", "ethnicity","hypoxia_ed_method","smoking")) |>
  janitor::clean_names() |>
  select(-race, -ethnicity, -hypoxia_ed_method, -smoking, -sex) |>
  rename_with(~str_c("baseline_", .), .cols = -empi)

all_wide <- 
  all_wide_pre_baseline |>
  left_join(baseline) |>
  select(empi, starts_with("baseline"), everything())


write_rds(all_wide, ("data/derived/all_wide.rds"))

