

produce_interim_vaccination_inputs <- function(vacc_timeseries_file) {
  
  # VE_timeseries <- read_csv(vacc_timeseries_file)
  VE_timeseries <- interpolated_vaccine_effect(vacc_timeseries_file)

  Ei_Et_timeseries <- VE_timeseries %>%
    mutate(mean_Ei = 1 - sqrt(effect),
           mean_Et = mean_Ei) %>%
    select(state, date, mean_Ei, mean_Et)
  
  source("vaccination/code/functions_age_classes.R")

  # NOTE: we need to vaccinate each jurisdiction at different times.
  # The forecasts start at different dates in each jurisdiction.
  dummy_dose_data <- tibble(
    state = c('ACT', 'NSW', 'NT', 'QLD', 'SA', 'TAS', 'VIC', 'WA'),
    date = ymd('2021-08-06', '2021-06-05', '2021-10-28', '2021-06-10',
               '2021-11-19', '2021-12-04', '2021-06-27', '2021-12-14')
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


delta_omicron_transition_dates <- function() {
    tribble(
        ~state, ~delta_until, ~omicron_from,
        'ACT', as.Date('2021-12-01'), as.Date('2022-01-01'),
        'NSW', as.Date('2021-12-01'), as.Date('2022-01-01'),
        'NT',  as.Date('2021-12-01'), as.Date('2022-01-01'),
        'QLD', as.Date('2021-12-01'), as.Date('2022-01-01'),
        'SA',  as.Date('2021-12-01'), as.Date('2022-01-01'),
        'TAS', as.Date('2021-12-01'), as.Date('2022-01-01'),
        'VIC', as.Date('2021-12-01'), as.Date('2022-01-01'),
        'WA',  as.Date('2022-12-31'), as.Date('2023-01-01'),
        )
}


interpolated_vaccine_effect <- function(vacc_timeseries_file) {
    VE_timeseries <- read_csv(vacc_timeseries_file) %>%
        mutate(effect_delta = if_else(is.na(effect_delta), 1, effect_delta))

    if (any(is.na(VE_timeseries))) {
        stop('NA values in VE_timeseries')
    }

    ratio_interpolation <- delta_omicron_transition_dates()

    VE_effect_timeseries <- inner_join(
        VE_timeseries, ratio_interpolation, by = 'state') %>%
        mutate(
            effect_change = effect_omicron - effect_delta,
            interp_numer = as.numeric(date - delta_until),
            interp_denom = as.numeric(omicron_from - delta_until),
            interp_scale = interp_numer / interp_denom,
            effect = case_when(
                date <= delta_until ~ effect_delta,
                date >= omicron_from ~ effect_omicron,
                TRUE ~ effect_delta + effect_change * interp_scale)) %>%
        select(state, date, effect, effect_delta, effect_omicron)

    write_csv(VE_effect_timeseries, vacc_timeseries_file)

    interpolation_plot <- ggplot(VE_effect_timeseries,
                                 aes(date, effect, colour = state)) +
        geom_line() +
        geom_line(mapping = aes(date, effect_delta),
                  linetype = 'dashed') +
        geom_line(mapping = aes(date, effect_omicron),
                  linetype = 'dotted') +
        xlab(NULL) +
        ylab('Vaccine effect') +
        scale_colour_hue(NULL) +
        theme_minimal()

    ggsave('VE_effect_interpolation.png', width = 8, height = 6, bg = 'white')

    VE_effect_timeseries %>% select(state, date, effect)
}
