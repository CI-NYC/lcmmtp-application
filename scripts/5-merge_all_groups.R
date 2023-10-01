### Merge all three groups (never, intubated first, AKI first) data into one wide DF
### Kat Hoffman

library(tidyverse)

nevers_wide <- read_rds(here::here("data/derived/nevers_wide.rds")) |> mutate(group = "nevers")
aki_first_wide <- read_rds(here::here("data/derived/aki_first_wide.rds")) |> mutate(group = "aki_first")
int_first_wide <- read_rds(here::here("data/derived/int_first_wide.rds")) |> mutate(group = "int_first")

all_wide <-
  reduce(list(nevers_wide, aki_first_wide, int_first_wide),
       ~bind_rows(.x, .y))

write_rds(all_wide, ("data/derived/all_wide.rds"))

