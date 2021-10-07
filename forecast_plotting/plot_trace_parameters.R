

plot_trace_parameters <- function(dt_statevec,
                                  date_min, date_max) {
  dt_statevec_plot <- dt_statevec %>%
    filter(date >= date_min, date <= date_max) %>%
    filter(index %% 200 == 0)
  
  
  ggplot() +
    geom_line(aes(x = date, y = value),
              dt_statevec_plot) +
    
    facet_wrap(~state,
               ncol = 2,
               scales = "free_y")
  
  
}