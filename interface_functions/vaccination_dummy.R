

produce_interim_vaccination_inputs <- function(vacc_timeseries_file) {
  
  VE_timeseries <- read_csv(vacc_timeseries_file)
  
  
  Ei_Et_timeseries <- VE_timeseries %>%
    mutate(mean_Ei = 1 - sqrt(effect),
           mean_Et = mean_Ei) %>%
    select(state, date, mean_Ei, mean_Et)
  
  source("vaccination/code/functions_age_classes.R")
  
  dummy_dose_data <- expand_grid(
    state = unique(Ei_Et_timeseries$state),
    date = ymd("2021-01-01")
  ) %>%
    
    left_join(state_populations) %>%
    rename(rate = population)
  
  resultant_data <- expand_grid(
    state = unique(Ei_Et_timeseries$state),
    date = seq(ymd("2020-01-01"), max(Ei_Et_timeseries$date), by = 'days')
  ) %>%
    
    left_join(dummy_dose_data) %>%
    mutate(rate = if_else(is.na(rate), 0, rate)) %>%
    
    left_join(Ei_Et_timeseries) %>%
    
    group_by(state) %>%
    fill(c(mean_Ei, mean_Et), .direction =  "updown") %>%
    ungroup() %>%
    
    select(date, state, rate, mean_Ei, mean_Et)
    
  
  resultant_data %>%
    write_csv("vaccination/vaccination_output_data.csv")
  
  
  
  for(i_state in unique(resultant_data$state)){
    file_out <- paste0("data/vacc-data-",tolower(i_state), ".ssv")
    
    print(paste0("Writing ", file_out))
    resultant_data %>%
      filter(state == i_state) %>%
      select(-state) %>%
      write_delim(file_out)
    
    
  }
    
}