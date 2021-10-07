

produce_Ei_Et_timeseries <- function(vacc_dose_data_file,
                                     vacc_timeseries_file,
                                     run_name) {
  
  require(tidyverse)
  require(lubridate)
  
  diag_plot_dir <- paste0("results/", run_name, "/diagnostics")
  
  dir.create(diag_plot_dir, recursive = TRUE)
  
  
  source("vaccination/code/functions_age_classes.R")
  source("vaccination/code/functions_matrices.R")
  
  dose_data <- read_csv(vacc_dose_data_file) %>%
    mutate(date = date + 6)
  
  # ATAGI Delta estimates as of 2021-08-19
  assign_vacc_Ei <- function(vaccine_type, dose) {
    case_when(vaccine_type == "pf" & dose == 1 ~ 0.3,
              vaccine_type == "pf" & dose == 2 ~ 0.79,
              vaccine_type == "az" & dose == 1 ~ 0.18,
              vaccine_type == "az" & dose == 2 ~ 0.6)
  }
  assign_vacc_Et <- function(vaccine_type, dose) {
    case_when(vaccine_type == "pf" & dose == 1 ~ 0.46,
              vaccine_type == "pf" & dose == 2 ~ 0.65,
              vaccine_type == "az" & dose == 1 ~ 0.48,
              vaccine_type == "az" & dose == 2 ~ 0.65)
  }
  
  assign_vacc_overall<- function(vaccine_type, dose) {
    Ei <- assign_vacc_Ei(vaccine_type, dose)
    Et <- assign_vacc_Et(vaccine_type, dose)
    
    1 - (1 - Ei) * (1 - Et)
  }
  
  dose_data_with_frac <- dose_data %>%
    group_by(state, date, age_class) %>%
    mutate(effective_doses_any_vacc = sum(effective_doses)) %>%
    left_join(age_distribution_by_state) %>%
    mutate(effective_fraction_any_vacc = effective_doses_any_vacc / population,
           effective_proportion = effective_doses / effective_doses_any_vacc)
  
  efficacy_data <- dose_data_with_frac %>%
    group_by(state, date, age_class) %>%
    summarise(effective_mean_Ei = sum(effective_proportion * assign_vacc_Ei(vaccine, dose_number), na.rm = TRUE),
              effective_mean_Et = sum(effective_proportion * assign_vacc_Et(vaccine, dose_number), na.rm = TRUE),
              
              effective_overall = sum(effective_proportion * assign_vacc_overall(vaccine, dose_number), na.rm = TRUE),
              
              
              effective_coverage_any_vaccine = first(effective_fraction_any_vacc),
              age_class_frac = first(fraction)) %>%
    drop_na(effective_coverage_any_vaccine)
  
  calc_vaccination_effect <- function(
    age_coverage,
    mean_Et,
    mean_Ei,
    mean_overall,
    next_generation_matrix,
    fraction,
    age_class
  ) {
    
    print(age_class)
    
    # For comparison, the Reff model estimate method:
    NG_effect <- 1 - age_coverage * mean_overall
    
    
    NG_next_gen <-  sweep(
      next_generation_matrix,
      2,
      NG_effect,
      FUN = "*"
    )
    overall_NG <- get_R(NG_next_gen) / get_R(next_generation_matrix)
    
    
    
    ## Our method:
    
    age_transmission_reduction <- 1 - age_coverage * mean_Et
    age_infection_reduction <- 1 - age_coverage * mean_Ei
    
    
    inf_next_gen_matrix <- sweep(
      next_generation_matrix,
      2,
      age_infection_reduction,
      FUN = "*"
    )
    inf_overall <- get_R(inf_next_gen_matrix) / get_R(next_generation_matrix)
    
    trans_next_gen_matrix <- NG_next_gen
    trans_overall <- get_R(trans_next_gen_matrix) / get_R(inf_next_gen_matrix)
    
    
    mean_Et_adj <- (1 - trans_overall) / sum(age_coverage * fraction)
    mean_Ei_adj <- (1 - inf_overall) / sum(age_coverage * fraction)
    
    overall <- inf_overall *  trans_overall
    
    naive_Et = sum(fraction * mean_Et)
    naive_Ei = sum(fraction * mean_Ei)
    
    tibble_row(
      overall_NG = overall_NG,
      overall = overall,
      mean_Et_adj = mean_Et_adj,
      mean_Ei_adj = mean_Ei_adj,
      naive_Et = naive_Et,
      naive_Ei = naive_Ei
    )
  }
  
  vaccination_effect <- efficacy_data %>%
    
    # Important: ensure ordering by age class
    mutate(age_first_val = as.numeric(str_extract(age_class, "^[0-9]{1,2}(?=[-,+])"))) %>%
    arrange(age_first_val) %>% select(-age_first_val) %>%
    
    group_by(
      state, date
    ) %>%
    summarise(
      transmission_effect = calc_vaccination_effect(
        age_coverage = effective_coverage_any_vaccine,
        mean_Ei = effective_mean_Ei,
        mean_Et = effective_mean_Et,
        mean_overall = effective_overall,
        next_generation_matrix = baseline_matrix(),
        fraction = age_class_frac,
        age_class
      ),
      .groups = "drop"
    ) %>%
    unpack(transmission_effect)
  
  
  
  vacc_timeseries <- read_csv(vacc_timeseries_file)
  
  
  ggplot() +
    
    geom_line(aes(x = date, y = overall_NG, color = 'our calculation of TP model'),
              vaccination_effect) +
    
    geom_line(aes(x = date, y = effect, color = 'provided timeseries'),
              linetype = 'dotted',
              vacc_timeseries) +
    
    geom_line(aes(x = date, y = overall, color = 'our method'),
              linetype = 'dashed',
              vaccination_effect) +
    
    facet_wrap(~state) +
    theme_minimal() +
    theme(legend.position = 'bottom')
  
  ggsave(paste0(diag_plot_dir, "/vaccine_effect_comparison.png"),
         device = ragg::agg_png, bg = 'white',
         width = 8, height = 6)
  
  ggplot() +
    geom_line(aes(x = date, y = mean_Ei_adj, color = 'Ei'),
              vaccination_effect) +
    geom_line(aes(x = date, y = naive_Ei, color = 'Ei', linetype = 'naive'),
              vaccination_effect) +
    
    geom_line(aes(x = date, y = mean_Et_adj, color = 'Et', linetype = 'adjusted'),
              vaccination_effect) +
    geom_line(aes(x = date, y = naive_Et, color = 'Et', linetype = 'naive'),
              vaccination_effect) +
    
    facet_wrap(~state) +
    theme_minimal() +
    theme(legend.position = 'bottom')
  
  ggsave(paste0(diag_plot_dir, "/Ei_Et_naive_adjusted_comparison.png"),
         device = ragg::agg_png, bg = 'white',
         width = 8, height = 6)
  
  output_data <- vaccination_effect %>%
    select(date, state, mean_Ei = mean_Ei_adj, mean_Et = mean_Et_adj, overall = overall)
  
  output_data
}

