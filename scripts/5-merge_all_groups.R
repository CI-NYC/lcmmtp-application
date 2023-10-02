#######################################################
#######################################################
### Merge vars for all three groups (never, intubated first, AKI first) data into one wide DF
### Kat Hoffman
#######################################################
#######################################################

library(tidyverse)

max_window <- 14 # max windows t for data set

na_list <- as.list(rep(0, max_window))
na_list
names(na_list) <- paste0("Observed_", 1:max_window)

fp <- "data/derived" # folder where derived data lives

int_first <- read_rds(here::here(fp, "int_first", "int_first.rds"))
aki_first <- read_rds(here::here(fp, "aki_first", "aki_first.rds"))
nevers <- read_rds(here::here(fp, "nevers", "nevers.rds"))
  
all_wide <-
  map_dfr(
  list("nevers", "aki_first", "int_first"),
  function(group) {
  
  Ls_and_Zs <- read_rds(here::here(fp, group, "Ls_and_Zs.rds"))
  As <- read_rds(here::here(fp, group, "As.rds"))
  Ms <- read_rds(here::here(fp, group, "Ms.rds"))
  Cs <- read_rds(here::here(fp, group, "Cs.rds"))
  Ys <- read_rds(here::here(fp, group, "Ys.rds"))
  
  # first make Ls and Zs wide, more complex
  Ls_and_Zs_wide <-
    Ls_and_Zs |>
    filter(window <= max_window) |>
    as_tibble() |>
    pivot_wider(id_cols=empi,
                names_from = c(window, covar),
                values_from = c(L_value, L_missing, Z_value, Z_missing))
  
  # compile the rest of the long data
  dat_long <- reduce(list(As, Ms, Cs, Ys), ~left_join(.x, .y)) |>
    filter(window <= max_window)
  
  # make the rest of the data wide
  dat_wide <- pivot_wider(dat_long, 
                          id_cols = empi,
                          names_from = window,
                          values_from = c(A, M, Observed, Y)) |>
    left_join(Ls_and_Zs_wide) |> # add Ls and Zs in
    mutate(group = group) |> # add group id for later debugging
    replace_na(na_list)
  
  return(dat_wide)
  }
)

write_rds(all_wide, ("data/derived/all_wide.rds"))

