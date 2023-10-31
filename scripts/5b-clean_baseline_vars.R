#######################################################
#######################################################
### Clean baseline demographics, add in, save final data set
### Kat Hoffman
#######################################################
#######################################################

library(tidyverse)

cohort <- read_rds(here::here(fp, "hospitalized_cohort_with_queens.rds")) # contains baseline data
all_wide_pre_baseline <- read_rds(here::here(fp, "all_wide_pre_baseline.rds"))

### prep baseline data to merge
baseline <-
  cohort |>
  select(empi, red_cap_source,
         age, sex, race, ethnicity, bmi, smoking, 
         cad, home_o2_yn, dm, htn, cva, cirrhosis, ckd_or_esrd, asthma, copd, active_cancer, 
         immunosuppressed, ild, hiv, hypoxia_ed) |>
  mutate(bmi_miss = case_when(is.na(bmi) ~ 1, TRUE ~ 0),
  home_o2_miss = case_when(is.na(home_o2_yn) ~ 1, TRUE ~ 0), 
  race_miss = case_when(is.na(race) ~ 1, TRUE ~ 0),
  ethnicity_miss = case_when(is.na(ethnicity) ~ 1, TRUE ~ 0),
  race = case_when(is.na(race) ~ "Missing", TRUE ~ race),
  ethnicity = case_when(is.na(ethnicity) ~ "Missing", TRUE ~ ethnicity),
  smoking = case_when(is.na(smoking) ~ "No", TRUE ~ smoking),
  male = case_when(sex=="Male"~1, TRUE ~0),
  red_cap_source = case_when(red_cap_source == "QUEENS" ~ 0,
                             red_cap_source == "EAST" ~ 1),
  across(cad:hiv, ~case_when(.x == "No" ~ 0, TRUE ~ 1))) |>
  fastDummies::dummy_columns(select_columns = c("race", "ethnicity","smoking")) |>
  janitor::clean_names() |>
  select(-race, -ethnicity,  -smoking, -sex) |>
  rename_with(~str_c("L_1_", .), .cols = -empi)

# merge with the time-varying dataset and save
all_wide <- 
  all_wide_pre_baseline |>
  left_join(baseline) |> 
  mutate(id = row_number()) |>
  select(id, starts_with("L_1_"), everything(), -empi)

write_rds(all_wide, ("data/derived/all_wide.rds"))
