######################################################
######################################################
## Wrangle labs and vitals into long clean data sets
## Kat Hoffman
######################################################
######################################################

# load packages -----------------------------------------------------------

library(lubridate)
library(tictoc)
library(dtplyr)
library(tidyverse)
library(tidylog)
library(zoo)
library(gtsummary)
library(gt)
library(labelled)
library(snakecase)
library(janitor)

# load data -----------------------------------------------------------

cp_dt <- read_rds(here::here("data/derived/dts_cohort.rds"))
labs <- read_rds(here::here("data", "raw", "2020-08-17", "wcm_labs_filtered.rds"))

# Define "worst" lab as low or high depending on lab

low_lab_params <- c(
  "Platelet", # decided to use low platelets as lowest (if two measures in 12 hours... unlikely) because we have all the cytostorm labs already
  "Absolute Neutrophil Count", 
  "Neutrophil Automated",
  "Glucose Whole Blood  Meter POC",
  "Glucose Whole Blood Meter POC",
  "pO2 (Arterial) - EPOC",
  "pO2 Arterial",
  "pCO2 (Arterial) - EPOC",
  "pCO2 Arterial",
  "P CO2 Arterial",
  "P O2 Arterial"
)

high_lab_params <- c(
  "Lymphocyte Automated",
  "Absolute Lymphocyte Count",
  "Fibrinogen",
  "D-Dimer",
  "C-Reactive Protein",
  "C-Reactive Protein High Sensitivity",
  "Ferritin",
  "Prothrombin Time",
  "Activated Partial Thromboplastin Time",
  "Troponin-I",
  "Glucose Whole Blood  Meter POC",
  "Glucose Whole Blood Meter POC",
  "Bilirubin Total",
  "BUN/Creatinine Ratio",
  "Creatine Kinase",
  "Creatinine"
)

labs_low_long <-
  labs %>%
  filter(result_name %in% low_lab_params) %>%
  select(empi, result_time, result_name, ord_value) %>%
  mutate(ord_value = parse_number(ord_value),
         result_name = case_when(str_detect(result_name, "Glucose") ~ "Glucose",
                                 str_detect(result_name, "Reactive") ~ "C-Reactive Protein",
                                 str_detect(result_name, "Reactive") ~ "C-Reactive Protein",
                                 str_detect(result_name, "Lactate") ~ "Lactate",
                                 str_detect(result_name, "pH") ~ "Arterial Ph",
                                 str_detect(result_name, "CO2") ~ "Arterial PCO2",
                                 str_detect(result_name, "O2") ~ "Arterial PaO2",
                                 str_detect(result_name, "Bilirubin") ~ "Bilirubin",
                                 TRUE ~ result_name)) |>
  select(empi, result_dt = result_time, result_value = ord_value, name = result_name)

labs_high_long <-
  labs %>%
  filter(result_name %in% high_lab_params) %>%
  select(empi, result_time, result_name, ord_value) %>%
  mutate(ord_value = parse_number(ord_value),
         result_name = case_when(str_detect(result_name, "Glucose") ~ "Glucose",
                                 str_detect(result_name, "Reactive") ~ "C-Reactive Protein",
                                 str_detect(result_name, "Reactive") ~ "C-Reactive Protein",
                                 str_detect(result_name, "Lactate") ~ "Lactate",
                                 str_detect(result_name, "pH") ~ "Arterial Ph",
                                 str_detect(result_name, "CO2") ~ "Arterial PCO2",
                                 str_detect(result_name, "O2") ~ "Arterial PO2",
                                 str_detect(result_name, "Bilirubin") ~ "Bilirubin",
                                 TRUE ~ result_name)) |>
  select(empi, result_dt = result_time, result_value = ord_value, name = result_name)

all_labs_clean <- bind_rows(labs_low_long, labs_high_long)

write_rds(all_labs_clean, here::here("data/derived/all_labs_clean.rds"))

