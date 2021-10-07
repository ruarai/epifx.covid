

make_obs_and_CI_plot <- function(dt_notifications,
                                 dt_forecasts,
                                 date_min, date_max,
                                 use_columns = TRUE,
                                 probs_plot = c(50, 95)) {
  
  dt_forecasts_plot <- dt_forecasts %>%
    filter(date >= date_min, date <= date_max)
  
  dt_notifications_plot <- dt_notifications %>%
    filter(date >= date_min, date <= date_max) %>%
    mutate(adj_value = value / pr_detect)
  
  notifications_columns <- list(
    geom_linerange(aes(x = date, ymin = 0, ymax = value),
                   dt_notifications_plot,
                   col = 'grey50',
                   size = 2)
  )
  
  notifications_line <- list(
    geom_line(aes(x = date, y = value),
                   dt_notifications_plot,
                   col = 'grey50')
  )
  
  
  ggplot() +
    
    ifelse(use_columns, notifications_columns, notifications_line) +
    
    geom_point(aes(x = date, y = adj_value),
                   dt_notifications_plot %>% filter(pr_detect < 0.95, adj_value > 0),
                   col = 'grey50',
               pch = 1,
               size = 1) +
    
    geom_ribbon(aes(x = date, ymin = ymin, ymax = ymax, fill = prob, color = prob, group = prob),
                dt_forecasts_plot %>% filter(prob %in% probs_plot),
                alpha = 0.5) +
    
    geom_line(aes(x = date, y = ymin),
              dt_forecasts_plot %>% filter(prob == 0), 
              alpha = 0.5,
              color = '#3182bd') +
    
    scale_fill_brewer(palette = "Blues") +
    
    scale_color_brewer(palette = "Blues") +
    
    expand_limits(y = 10) +
    scale_y_continuous(breaks = scales::breaks_extended(),
                       labels = scales::label_comma(),
                       position = 'right') +
    
    scale_x_date(breaks = scales::breaks_pretty(12),
                 labels = scales::label_date_short()) +
    
    xlab("") + ylab("") +
    
    theme_minimal() +
    theme(legend.position = 'none') +
    
    facet_wrap(~state,
               scales = "free_y",
               ncol = 2)
  
  
}
