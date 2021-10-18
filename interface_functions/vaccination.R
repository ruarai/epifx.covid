


produce_vaccination_inputs <- function(vacc_dose_data_file,
                                       vacc_timeseries_file,
                                       run_name) {
  source("vaccination/code/process_AIR_dose_data.R")
  source("vaccination/code/process_first_dose_timeseries.R")
  
  diag_plot_dir <- paste0("results/", run_name, "/diagnostics")
  
  Ei_Et_timeseries <- produce_Ei_Et_timeseries(vacc_dose_data_file,
                                               vacc_timeseries_file,
                                               run_name)
  
  
  first_dose_timeseries <- produce_first_dose_timeseries(vacc_dose_data_file,
                                                         run_name)
  
  
  resultant_Ei_Et_timeseries <- expand_grid(
    state = unique(Ei_Et_timeseries$state),
    date = seq(min(Ei_Et_timeseries$date), max(Ei_Et_timeseries$date), by = 'days')
  ) %>%
    left_join(Ei_Et_timeseries %>% select(-overall)) %>%
    mutate(mean_Ei = zoo::na.approx(mean_Ei, na.rm = FALSE),
           mean_Et = zoo::na.approx(mean_Et, na.rm = FALSE))
  
  
  ggplot() +
    geom_line(aes(x = date, y = mean_Ei, color = 'Ei'),
              resultant_Ei_Et_timeseries) +
    geom_line(aes(x = date, y = mean_Et, color = 'Et'),
              resultant_Ei_Et_timeseries) +
    
    facet_wrap(~state) +
    theme_minimal() +
    theme(legend.position = 'bottom') +
    
    coord_cartesian(ylim = c(0,0.8))
  
  ggsave(paste0(diag_plot_dir, "/Ei_Et_timeseries_final.png"),
         width = 8, height = 6)
  
  
  
  output_data <- first_dose_timeseries %>%
    left_join(resultant_Ei_Et_timeseries) %>%
    group_by(state) %>%
    fill(mean_Ei, mean_Et, .direction = 'downup') %>%
    select(date, state, rate = count_new, mean_Ei, mean_Et) %>%
    ungroup()
  
  output_data %>%
    write_csv("vaccination/vaccination_output_data.csv")
  
  
  
  for(i_state in unique(output_data$state)){
    file_out <- paste0("data/vacc-data-",tolower(i_state), ".ssv")
    
    print(paste0("Writing ", file_out))
    output_data %>%
      filter(state == i_state) %>%
      select(-state) %>%
      write_delim(file_out)
    
    
  }
  
}

