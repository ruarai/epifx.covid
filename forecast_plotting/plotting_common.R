

close_hdf_files <- function() {
  require(rhdf5)
  
  h5closeAll()
}

read_hdf_date <- function(date_time) {
  require(lubridate)
  
  date_time %>%
    fast_strptime("%Y-%m-%d %H:%M:%S", lt = FALSE) %>%
    as_date()
}

get_plotting_data <- function(forecast_exps_dir) {
  hdf_files <- tibble(file = list.files(forecast_exps_dir, full.names = TRUE)) %>%
    mutate(state = str_extract(file, "(?<=/).{2,3}(?=_FS-)"))
  
  if(nrow(hdf_files) != 8) print("Number of HDF files is not equal to 8.")
  
  require(rhdf5)
  
  dt_notifications <- hdf_files %>%
    pmap_dfr(function(file, state) {
      a <- h5read(file = file,
                  name = "/data/obs/Notifications") %>%
        select(date, value, pr_detect) %>% 
        mutate(date = read_hdf_date(date),
               state = state)
    })
  
  
  min_date <- dt_notifications %>% 
    pull(date) %>% min() %>% as_date()
  
  fs_date <- dt_notifications %>% 
    filter(pr_detect > 0.5) %>%
    pull(date) %>% max() %>% as_date()
  fs_horizon <- fs_date + 28
  
  states <- dt_notifications %>%
    pull(state) %>% unique()
  
  forecast_line <- geom_vline(xintercept = fs_date, 
                              linetype = 'dashed',
                              alpha = 0.8,
                              color = '#016c59')
  
  list(min_date = min_date,
       fs_date = fs_date,
       fs_horizon = fs_horizon,
       states = states,
       forecast_line = forecast_line,
       dt_notifications = dt_notifications,
       hdf_files = hdf_files,
       forecast_exps_dir = forecast_exps_dir)
}



