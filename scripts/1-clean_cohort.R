#######################################################
#######################################################
### Create cohort, apply eligibility criteria
### Kat Hoffman
#######################################################
#######################################################

library(lubridate)
library(tictoc)
library(dtplyr)
library(tidyverse)
library(tidylog)
library(zoo)
library(gtsummary)

# set ggplot theme
theme_set(theme_classic(base_family = "Roboto"))

file_date <- "2020-08-17"
file_path <- paste0("data/raw/",file_date)

# read in all redcap data, remove patients missing empi
redcap_raw <-
  read_rds(here::here(file_path, "wcm_redcap_data.rds")) %>%
  drop_na(empi)

aki <- read_rds(here::here("data/derived/kh_aki.rds"))
tz(aki$aki_time) <- "America/New_York"

# Apply inclusion criteria:
#   
#   1. Adult COVID-19 positive patients
#   2. Admitted to LMH, WCM, or Queens hospital (not just ED) and if admitted more than once, we keep their first visit.
#   3. Data available in the Data Lake from beginning of ED admission (i.e. not a transfer into the system)
#   4. No CKD or ESRD

redcap <-
  redcap_raw %>%
  # some charts not finished by abstractors yet
  filter(covid_complete != "Unverified")  %>%
  # only LMH and WCM
  # filter(red_cap_source == "EAST") %>%
  # No transfer patients
  filter(is.na(transfer_in_date)) %>%
  # Cleaning up redcap variables
  mutate(mrn = str_squish(mrn),
         age = as.numeric(age),
         bmi = as.numeric(bmi),
         active_cancer = case_when(str_detect(active_cancer, "No") ~ "No", TRUE ~ "Yes"),
         race = case_when(race == "Not Specified" ~ NA_character_, TRUE ~ race),
         ethnicity = case_when(ethnicity %in% c("Unknown or not specified", "Other") ~ NA_character_, TRUE ~ ethnicity),
         pulm_disease = case_when(pulm == "No" ~ "No", TRUE ~ "Yes"),
         renal_disease = case_when(renal == "No" ~ "No", TRUE ~ "Yes"),
         ckd_or_esrd = case_when(str_detect(renal,"CKD (Creat &gt;2 at baseline)") ~ "Yes",
                                 str_detect(renal,"ESRD") ~ "Yes",
                                 TRUE ~ "No"),
         copd = case_when(str_detect(pulm,"COPD") == T ~ "Yes", T ~ "No"),
         ild = case_when(str_detect(pulm,"Interstitial Lung Disease") == T ~ "Yes", T ~ "No"),
         asthma = case_when(str_detect(pulm,"Asthma") == T ~ "Yes", T ~ "No"),
         renal_disease = case_when(renal == "No" ~ "No", TRUE ~ "Yes"),
         hf = case_when(hf == "No" ~ "No", TRUE ~ "Yes"),
         hypoxia_ed_method = case_when(is.na(hypoxia_ed_method) ~ "None", TRUE ~ hypoxia_ed_method),
         hepatitis = case_when(hepatitis == "No" ~ "No", TRUE ~ "Yes"),
         first_cxr_results = case_when(
           str_detect(first_cxr_results, "Bilateral Infiltrates") ~ "Bilateral Infiltrates",
           str_detect(first_cxr_results, "Unilateral Infiltrates") ~ "Unilateral Infiltrates",
           str_detect(first_cxr_results, "Pleural Effusion") ~ "Pleural Effusion",
           first_cxr_results == "Clear" ~ "Clear",
           TRUE ~ "Other"),
         vte = case_when(vte == "No" ~ "No", TRUE ~ "Yes"),
         hypoxia_ed = case_when(str_detect(hypoxia_ed, "Yes") ~ "Yes",
                                TRUE ~ hypoxia_ed),
         hypoxia_ed_method = case_when(str_detect(hypoxia_ed_method, "Venti mask") ~ "Venti mask",
                                       TRUE ~ hypoxia_ed_method),
         home_o2_yn = case_when(str_detect(home_o2, "No") ~ "No",
                                is.na(home_o2) ~ NA_character_,
                                TRUE ~ "Yes"),
         # date time vars
         ed_adm_dt = paste0(ed_date, ed_time),
         admit_dt = paste0(admit_date, admit_time),
         intubation1_dt = paste(intubation1_date, intubation1_time),
         discharge_dt = paste0(discharge_date,
                               discharge_time),
         transfer_dt = paste0(transfer_out_date,
                              transfer_out_time),
         death_dt = paste0(death_date,
                           death_time),
         last_fu_dt = paste0(last_reviewed, "23:59")) %>%
  mutate_at(vars(ends_with("dt")), ~as.POSIXct(.,  format="%Y-%m-%d %H:%M", tz="America/New_York")) %>%
  mutate(hours_until_intubation = time_length(difftime(intubation1_dt, ed_adm_dt), unit="hour"),
         hours_until_intubation = time_length(difftime(intubation1_dt, ed_adm_dt), unit="hour")) %>% 
  arrange(empi, ed_adm_dt) %>%
  group_by(empi) %>%
  mutate(visit = 1, # for summing; meaningless variable
         visit_total = sum(visit),
         visit_n = cumsum(visit), # number of visits a patient has had (ever, admitted or not)
         # did patient discharge/transfer to outside hospital/die? if not, chart review is last follow up (fu) date
         end_dt = pmin(discharge_dt, death_dt, na.rm=T),
         end_dt = case_when(is.na(end_dt) ~ last_fu_dt, TRUE ~ end_dt)) %>%
  # inclusion: all patients admitted to the hospital
  filter(age >= 18,
    admit == "Yes",
    #pregnancy != "Yes",
  ) %>%
  mutate(adm_visit_n = cumsum(visit)) %>%  # number of admitted visits a patient has had
  arrange(adm_visit_n) %>%
  # keep only their first admitted visit in the GIM database
  filter(row_number() == 1) %>%
  mutate(day_death = ceiling(time_length(difftime(death_dt, ed_adm_dt), unit="day"))) %>%
  ungroup()

# apply CKD exclusion criteria --------------------------------------------

# remove pts with ckd or esrd already
cohort <-
  redcap %>%
  filter(ckd_or_esrd == "No") |>
  left_join(aki)

# save data
saveRDS(cohort, here::here("data","derived","hospitalized_cohort_with_queens.rds"))