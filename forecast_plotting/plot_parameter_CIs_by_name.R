
make_parameter_CI_plot_by_name <- function(dt_parameter_CIs,
                                   param_name,
                                   date_min, date_max) {
  
  
  dt_parameter_plot <- dt_parameter_CIs %>%
    filter(date >= date_min, date <= date_max,
           name == param_name,
           prob %in% c(95, 50, 0))
  
  
  ggplot() +
    geom_ribbon(aes(x = date, ymin = ymin, ymax = ymax, fill = prob, group = prob, color = prob),
                dt_parameter_plot %>% filter(prob != 0),
                alpha = 0.5) +
    
    geom_line(aes(x = date, y = ymin),
              dt_parameter_plot %>% filter(prob == 0), 
              alpha = 0.5,
              color = '#3182bd') +
    
    facet_grid(rows = vars(state),
               switch = "both") +
    
    scale_fill_brewer(palette = "Blues") +
    
    scale_color_brewer(palette = "Blues") +
    
    scale_y_continuous(breaks = scales::breaks_extended(),
                       labels = scales::label_comma(),
                       position = 'right') +
    
    scale_x_date(breaks = scales::breaks_pretty(12),
                 labels = scales::label_date_short()) +
    
    xlab("") + ylab("") +
    
    theme_minimal() +
    theme(legend.position = 'none')
  
  
}