

make_obs_and_CI_plot <- function(dt_notifications,
                                 dt_forecasts,
                                 date_min, date_max,
                                 use_columns = TRUE,
                                 probs_plot = c(50, 95)) {
  
  dt_forecasts_plot <- dt_forecasts %>%
    filter(date >= date_min, date <= date_max)

  dt_start <- aggregate(date ~ state, data = dt_forecasts_plot, FUN = min) %>%
    rename(start_date = date)
  
  dt_notifications_plot <- dt_notifications %>%
    filter(date >= date_min, date <= date_max) %>%
    left_join(dt_start) %>%
    group_by(state) %>%
    filter(date >= start_date) %>%
    ungroup() %>%
    mutate(start_date = NULL) %>%
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

  draw_nothing <- list(
    geom_blank(aes(x = date, y = value),
               data = data.frame(date = dt_start$start_date, value = 0, state = dt_start$state))
  )
  
  
  ggplot() +
    
    ifelse(use_columns, notifications_columns, draw_nothing) +

    geom_ribbon(aes(x = date, ymin = ymin, ymax = ymax, fill = prob, color = prob, group = prob),
                dt_forecasts_plot %>% filter(prob %in% probs_plot),
                alpha = 0.5) +
    
    geom_line(aes(x = date, y = ymin),
              dt_forecasts_plot %>% filter(prob == 0), 
              alpha = 0.5,
              color = '#3182bd') +
    
    ifelse(use_columns, draw_nothing, notifications_line) +

    geom_point(aes(x = date, y = adj_value),
                   dt_notifications_plot %>% filter(pr_detect < 0.95, adj_value > 0),
                   col = 'grey50',
               pch = 1,
               size = 1) +

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
               scales = "free",
               ncol = 2)
  
  
}
