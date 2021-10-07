
process_local_cases <- function(local_cases_file, 
                                pr_detect_threshold = 0.5,
                                output_path = "data") {
  
  local_cases_cols <- cols(
    date_onset = col_date(format = ""),
    detection_probability = col_double(),
    state = col_character(),
    count = col_double(),
    acquired_in_state = col_double()
  )
  
  local_cases <- read_csv(local_cases_file,
                          col_types = local_cases_cols)
  
  
  local_cases %>%
    filter(detection_probability > pr_detect_threshold) %>%
    filter(date_onset >= max(date_onset ) - 3) %>%
    arrange(state, desc(date_onset)) %>%
    print(n = 100)
  
  
  states <- unique(local_cases$state)
  
  for(i_state in states) {
    out_file <- paste0(output_path, "/daily-covid-cases-50-", tolower(i_state), ".ssv")
    
    out_data <- local_cases %>%
      filter(state == i_state,
             detection_probability > pr_detect_threshold) %>%
      select(state, date = date_onset, cases = count, pr_detect = detection_probability)
    
    print(paste0("Writing ", out_file ))
    
    out_data %>%
      write_delim(out_file)
    
  }
}


process_reff_trajs <- function(reff_proj_file,
                               vaccine_effect_file,
                               output_path = "data",
                               truncation_date = NA) {
  
  reff_col_types <- cols(
    date = col_date(format = ""),
    state = col_character(),
    date_onset = col_date(format = ""),
    .default = col_double()
  )
  
  ve_col_types <- cols(
    state = col_character(),
    date = col_date(format = ""),
    effect = col_double()
  )
  
  reff_proj_data <- vroom::vroom(reff_proj_file,
                                 col_types = reff_col_types)
  
  vaccine_effect_data <- read_csv(vaccine_effect_file,
                                  col_types = ve_col_types)
  
  devacc_reff <- reff_proj_data %>%
    left_join(vaccine_effect_data) %>%
    group_by(state) %>%
    fill(effect, .direction = 'downup') %>%
    mutate(across(starts_with("sim"), ~ . / effect)) %>%
    ungroup() %>%
    select(-effect) %>%
    select(1:1003)
  
  if(!is.na(truncation_date)){
    print(paste0("Cutting off at ", truncation_date))
    
    devacc_reff <- devacc_reff %>%
      filter(date <= truncation_date)
  }
  
  if(any(is.na(devacc_reff))) {
    stop("Reff trajectories data contains NA values.")
  }
  
  all_state_row_counts_equal <- devacc_reff %>%
    group_by(state) %>%
    summarise(nrow = nrow(.)) %>%
    summarise(equal_row_count = var(nrow) == 0) %>%
    pull(equal_row_count)
  
  if(!all_state_row_counts_equal) {
    stop("Reff trajectories of different length between states.")
  }
  
  states <- unique(devacc_reff$state)
  
  for(i_state in states) {
    out_file <- paste0(output_path, "/reff-proj-", tolower(i_state), ".ssv")
    
    
    out_data <- devacc_reff %>%
      filter(state == i_state) %>%
      select(-state, - date_onset)
    
    print(paste0("Writing ", out_file ))
    
    out_data %>%
      write_delim(out_file)
  }
}

process_external_exposures <- function(ext_exps_file,
                                       output_path = "data") {
  ext_exps_col_types <- cols(
    jurisdiction = col_character(),
    date = col_character(),
    value = col_double()
  )
  
  
  ext_exps_all <- read_csv(ext_exps_file,
                           col_types = ext_exps_col_types)
  
  jurisdictions <- c("wa", "vic", "nt", "nsw", "tas", "qld", "act", "sa")
  
  
  for(jurisdiction_i in jurisdictions) {
    i_data <- ext_exps_all %>%
      filter(jurisdiction == jurisdiction_i)
    
    i_table <- tibble(date = "2020-07-01", value = 0)
    
    if(nrow(i_data) > 0){
      for(j_row in 1:nrow(i_data)) { 
        date_case <- dmy(i_data[j_row, c("date")])
        date_formatted <- format(date_case, "%Y-%m-%d")
        date_plus_one <- format(date_case + days(1), "%Y-%m-%d")
        
        value_case <- i_data[j_row, ]$value
        
        j_data <- tribble(~date         , ~value,
                          date_formatted, value_case,
                          date_plus_one , 0          )
        
        
        
        i_table <- bind_rows(i_table,
                             j_data)
      }
    }
    print(paste0(jurisdiction_i,":"))
    print(i_table)
    
    i_file <- paste0(output_path,"/daily-external-exposures-", jurisdiction_i, ".ssv")
    
    write.table(i_table,
                i_file, 
                quote = FALSE, 
                row.names = FALSE)
  }
}


get_forecast_dates <- function(local_cases_file) {
  local_cases <- read_csv(local_cases_file)
  
  date_minimum_onset <- local_cases %>%
    pull(date_onset) %>%
    min()
  
  date_last_onset_50 <- local_cases %>%
    filter(detection_probability > 0.5) %>%
    pull(date_onset) %>% max()
  
  date_last_infection_50 <- date_last_onset_50 - 5
  
  date_forecast_horizon = date_last_onset_50 + 28
  
  tibble(
    date_minimum_onset = date_minimum_onset,
    date_last_onset_50 = date_last_onset_50,
    date_last_infection_50 = date_last_infection_50,
    date_forecast_horizon = date_forecast_horizon,
  )
}


