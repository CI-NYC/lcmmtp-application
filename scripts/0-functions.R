# A function to define L_t and Z_t by covariate
# argument covar is a string containing lab name, eg "ddimer"
# window_num indicates what t we're determining
# returns a list for every time point and L_t and Z_t values

divide_labs <- function(df, covar, window_num){
  
  # extract labs within the window of interest
  rel_labs <-
    df |>
    filter(window == window_num) |>
    left_join(all_labs_clean |> filter(name == covar)) |>
    filter(result_dt %within% interval(l_start, z_end))
  
  # if in L interval, record as such
  Ls <- rel_labs |>
    drop_na(result_value) |>
    filter(result_dt %within% interval(l_start, l_end)) |>
    arrange(result_dt) |>
    distinct(empi, .keep_all = T) |>
    mutate(missing = ifelse(is.na(result_value), 1, 0)) |>
    select(empi, window, L_value = result_value)
  
  # if in Z interval, document as such
  Zs <-  rel_labs |>
    drop_na(result_value) |>
    filter(result_dt %within% interval(z_start, z_end)) |>
    arrange(result_dt) |>
    distinct(empi, .keep_all = T) |>
    select(empi, window, Z_value = result_value)
  
  out <-
    df |>
    select(empi, window) |>
    filter(window == window_num) |>
    mutate(covar = covar) |>
    left_join(Ls) |>
    left_join(Zs)
  
  return(out)
  
}