



plot_trace_incidence <- function(dt_sim_obs,
                                 date_min, date_max) {
  dt_sim_obs_plot <- dt_sim_obs %>%
    filter(date >= date_min, date <= date_max) %>%
    filter(index %% 20 == 0)
  
  
  ggplot() +
    geom_line(aes(x = date, y = value, group = index),
              alpha = 0.2,
              size = 0.1,
              dt_sim_obs_plot) +
    
    theme_minimal() +
    
    facet_wrap(~state,
               ncol = 2,
               scales = "free_y") +
    scale_y_continuous(breaks = scales::breaks_extended(),
                       labels = scales::label_comma(),
                       position = 'right') +
    
    scale_x_date(breaks = scales::breaks_pretty(12),
                 labels = scales::label_date_short()) +
    
    xlab("") + ylab("")
}