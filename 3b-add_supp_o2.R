######################################################
######################################################
## Add supplemental oxygen into 
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
vitals <- read_rds(here::here("data", "raw", "2020-08-17", "wcm_vitals_filtered_rel.rds"))

o2_clean <- 
  cp_dt %>%
  select(empi, t1_start, max_time) %>%
  left_join(vitals) %>%
  filter(observation_name %in%
           c("resp_non vent device",
             "resp_insp_gas_label_oxygen",
             "vs_resp_devicepacu"),
         recordeddtm %within% interval(t1_start, max_time)) %>%
  mutate(any_o2 = case_when(
    observation_name == "resp_insp_gas_label_oxygen" &
      parse_number(valuetext) > 0 ~ "Yes",
    observation_name == "resp_non vent device" &
      ((valuetext != "0") |  !str_detect(tolower(valuetext), "room"))  ~ "Yes",
    observation_name == "vs_resp_devicepacu" &
      !str_detect(tolower(valuetext), "room")  ~ "Yes",
    TRUE ~ "No"
  )) %>%
  group_by(empi) %>%
  select(empi, recordeddtm, any_o2) %>%
  filter(any_o2 == "Yes") |>
  mutate(any_o2 = 1)


first_o2 <- 
  o2_clean  |>
  arrange(recordeddtm) |>
  filter(row_number() == 1) |>
  rename(first_o2 = recordeddtm)

write_rds(first_o2, "data/derived/first_o2.rds")
