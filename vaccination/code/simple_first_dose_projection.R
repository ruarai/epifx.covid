project_first_doses <- function(AIR_first_dose_counts_interp, over12_pop) {
  
  
  simple_logit_models <- AIR_first_dose_counts_interp %>%
    left_join(over12_pop) %>%
    mutate(prop = count / population,
           t = as.numeric(date - ymd("2020-01-01"))) %>%
    
    filter(date >= max(date) - 7) %>%
    
    do(model = glm(prop ~ t, family = 'quasibinomial', data = .))
  
  
  simple_logit_lookup <- simple_logit_models$model
  names(simple_logit_lookup) <- simple_logit_models$state
  
  
  pred_data <- map_dfr(unique(AIR_first_dose_counts$state),
                       function(i_state) {
                         prediction_data <- expand_grid(
                           date = seq(max(AIR_first_dose_counts$date) + 1, max(AIR_first_dose_counts$date) + 28, by = 'day')
                         ) %>%
                           mutate(t = as.numeric(date - ymd("2020-01-01"))) %>%
                           mutate(pred_prop = predict(simple_logit_lookup[i_state][[1]],
                                                      newdata = .,
                                                      type = 'response'),
                                  state = i_state)
                       }) %>%
    select(state, date, pred_prop) %>%
    left_join(over12_pop) %>%
    mutate(count = pred_prop * population)
  
  
  
  combined_pred_data <- pred_data %>%
    bind_rows(AIR_first_dose_counts_interp) %>%
    
    ungroup() %>%
    arrange(state, date) %>%
    
    group_by(state) %>%
    mutate(count_new = count - lag(count, default = 0))
  
  
  
  
  
  ggplot() +
    geom_line(aes(x = date, y = count),
              alpha = 0.5, linetype='dotted',
              AIR_first_dose_counts_interp) +
    
    geom_line(aes(x = date, y = count),
              pred_data)  +
    
    geom_hline(aes(yintercept = population, color = 'Over 12'),
               over12_pop)  +
    
    geom_hline(aes(yintercept = population, color = 'Total population'),
               state_populations) +
    
    geom_vline(xintercept = max(AIR_first_dose_counts$date),
               alpha = 0.5, linetype = 'dashed') +
    
    scale_y_continuous(breaks = scales::breaks_extended(),
                       labels = scales::label_number()) +
    
    scale_color_brewer(type = 'qual', palette = 7) +
    
    xlab("Date (of effective protection)") + ylab("Number of people with at least 1 dose") +
    
    facet_wrap(~state,
               scale = "free_y") +
    
    coord_cartesian(xlim = c(today() - 60, today() + 40)) +
    
    theme_minimal() +
    theme(legend.position = 'bottom')
  
  
  ggsave(paste0(diag_plot_dir,"/dose_timeseries_pop_proj_logit.png"),
         bg = 'white',
         width = 8, height = 8)
  
  combined_pred_data
}


