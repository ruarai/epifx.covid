
plot_vaccine_effect <- function(plotting_data, plot_subdir, ve_ts_file) {
  
  require(tidyverse)
  require(lubridate)
  require(rhdf5)
  require(Cairo)
  require(cowplot)
  require(furrr)
  
  vaccine_effect_timeseries <- read_csv(ve_ts_file)
  
  vaccine_effect_timeseries <- expand_grid(date = seq(min(plotting_data$dt_notifications$date),
                                                      max(plotting_data$dt_notifications$date) + 70,
                                                      by = 'days'),
                                           state = plotting_data$states) %>%
    arrange(desc(date)) %>%
    left_join(vaccine_effect_timeseries) %>%
    group_by(state) %>%
    fill(effect, .direction = 'updown')
  
  plan(multisession, workers = 8)
  
  dt_statevec <- plotting_data$hdf_files %>%
    future_pmap_dfr(function(file, state) {
      h5read(file = file,
             name = "/data/seeiir") %>%
        mutate(E_U = E1_U + E2_U,
               E_V = E1_V + E2_V,
               I_U = I1_U + I2_U,
               I_V = I1_V + I2_V,
               n = S_U + S_V + E_U + E_V + I_U + I_V + R_U + R_V,
               prop_vacc = (S_V + E_V + I_V + R_V) / n,
               r_vacc = (1 - prop_vacc * mean_Ei) * (1 - prop_vacc * mean_Et)) %>%
        mutate(date = read_hdf_date(date),
               R_with_vacc = R_eff * r_vacc,
               state = !!state) %>%
        left_join(vaccine_effect_timeseries) %>%
        mutate(R_original = R_eff * effect)  %>%
        select(state, date, weight, ix, R_original, R_eff, mean_Et, mean_Ei, prop_vacc, r_vacc, R_with_vacc) %>%
        pivot_longer(cols = -c(state, date,ix,weight)) %>%
        
        group_by(state, name, date) %>%
        summarise(lower_90 = quantile(value,0.05,na.rm=T),
                  upper_90 = quantile(value,0.95,na.rm=T),
                  lower_50 = quantile(value,0.25,na.rm=T),
                  upper_50 = quantile(value,0.75,na.rm=T),
                  median = quantile(value,0.5,na.rm=T))
    })
  
  
  
  
  horizon_line <- geom_vline(xintercept = plotting_data$fs_date + 28, 
                             linetype = 'dashed',
                             alpha = 0.8,
                             color = 'grey40')
  
  CairoPDF(paste0(plot_subdir, "/vaccination.pdf"), 
           width = 12, height = 14)
  
  for(state_i in plotting_data$states) {
    
    plot_data_statevec_CI <- dt_statevec %>%
      filter(state == state_i)
    
    color_pal <- c("R_original" = "#66c2a5",
                   "R_eff" = "#fc8d62",
                   "R_with_vacc" = "#8da0cb")
    
    
    p1 <- ggplot() +
      geom_hline(yintercept = 1, linetype = 'longdash') +
      geom_line(aes(x = date, y = median, color = 'R_eff'),
                plot_data_statevec_CI %>% filter(name == "R_eff")) +
      geom_ribbon(aes(x = date, ymin = lower_50, ymax = upper_50, fill = 'R_eff'),
                  plot_data_statevec_CI %>% filter(name == "R_eff"),
                  alpha = 0.3) +
      
      
      geom_line(aes(x = date, y = median, color = 'R_original'),
                plot_data_statevec_CI %>% filter(name == "R_original")) +
      geom_ribbon(aes(x = date, ymin = lower_50, ymax = upper_50, fill = 'R_original'),
                  plot_data_statevec_CI %>% filter(name == "R_original"),
                  alpha = 0.5) +
      
      ylab("R") + xlab("") +
      
      scale_x_date(breaks = scales::breaks_pretty(12),
                   labels = scales::label_date_short()) +
      coord_cartesian(xlim = c(dmy("01-Feb-2021"), plotting_data$fs_horizon)) +
  
      scale_y_continuous(breaks = scales::breaks_extended(),
                         labels = scales::label_comma(),
                         position = 'right') +
      
      scale_color_manual("Trajectories",
                         values = color_pal,
                         labels = c("R_original" = "R from C1/2",
                                    "R_eff" = "R without vaccination effect")) +
      
      scale_fill_manual(values = color_pal) +
      theme_minimal() +
      theme(legend.position = 'bottom') +
      guides(fill = 'none') +
      plotting_data$forecast_line +
      horizon_line +
      ggtitle(paste0("Population-level vaccination measures ", state_i))
    
    
    p2 <- ggplot() +
      geom_hline(yintercept = 1, linetype = 'longdash') +
      geom_line(aes(x = date, y = median, color = 'R_eff'),
                plot_data_statevec_CI %>% filter(name == "R_eff")) +
      geom_ribbon(aes(x = date, ymin = lower_50, ymax = upper_50, fill = 'R_eff'),
                  plot_data_statevec_CI %>% filter(name == "R_eff"),
                  alpha = 0.3) +
      
      
      geom_line(aes(x = date, y = median, color = 'R_with_vacc'),
                plot_data_statevec_CI %>% filter(name == "R_with_vacc")) +
      geom_ribbon(aes(x = date, ymin = lower_50, ymax = upper_50, fill = 'R_with_vacc'),
                  plot_data_statevec_CI %>% filter(name == "R_with_vacc"),
                  alpha = 0.3) +
      
      ylab("R") + xlab("") +
      
      scale_x_date(breaks = scales::breaks_pretty(12),
                   labels = scales::label_date_short()) +
      coord_cartesian(xlim = c(dmy("01-Feb-2021"), plotting_data$fs_horizon)) +
      scale_y_continuous(breaks = scales::breaks_extended(),
                         labels = scales::label_comma(),
                         position = 'right') +
      
      scale_color_manual("Trajectories",
                         values = color_pal,
                         labels = c("R_with_vacc" = "R with EpiFX vaccination effect",
                                    "R_eff" = "R without vaccination effect")) +
      
      scale_fill_manual(values = color_pal) +
      
      theme_minimal() +
      theme(legend.position = 'bottom') +
      guides(fill = 'none') +
      plotting_data$forecast_line +
      horizon_line
    
    p3 <- ggplot() +
      geom_hline(yintercept = 1, linetype = 'longdash') +
  
      
      
      geom_line(aes(x = date, y = median, color = 'R_with_vacc'),
                plot_data_statevec_CI %>% filter(name == "R_with_vacc")) +
      geom_ribbon(aes(x = date, ymin = lower_50, ymax = upper_50, fill = 'R_with_vacc'),
                  plot_data_statevec_CI %>% filter(name == "R_with_vacc"),
                  alpha = 0.3) +
  
      geom_line(aes(x = date, y = median, color = 'R_original'),
                plot_data_statevec_CI %>% filter(name == "R_original")) +
      geom_ribbon(aes(x = date, ymin = lower_50, ymax = upper_50, fill = 'R_original'),
                  plot_data_statevec_CI %>% filter(name == "R_original"),
                  alpha = 0.5) +
      
      ylab("R") + xlab("") +
      
      scale_x_date(breaks = scales::breaks_pretty(12),
                   labels = scales::label_date_short()) +
      coord_cartesian(xlim = c(dmy("01-Feb-2021"), plotting_data$fs_horizon)) +
      scale_y_continuous(breaks = scales::breaks_extended(),
                         labels = scales::label_comma(),
                         position = 'right') +
      
      scale_color_manual("Trajectories",
                         values = color_pal,
                         labels = c("R_original" = "R from C1/2",
                                    "R_with_vacc" = "R with EpiFX vaccination effect")) +
      
      scale_fill_manual(values = color_pal) +
      
      theme_minimal() +
      theme(legend.position = 'bottom') +
      guides(fill = 'none') +
      plotting_data$forecast_line +
      horizon_line
    
    
    params_of_interest <- c("mean_Ei",
                            "mean_Et",
                            "prop_vacc",
                            "r_vacc")
    
    
    plot_data_params <- plot_data_statevec_CI %>%
      filter(name %in% params_of_interest) %>%
      mutate(name = factor(name, levels = params_of_interest))
    
    p4 <- ggplot() +
      geom_line(aes(x = date, y = median),
                plot_data_params) +
      
      facet_grid(rows = vars(name)) +
      
      ylab("Value") + xlab("") +
      
      scale_x_date(breaks = scales::breaks_pretty(12),
                   labels = scales::label_date_short()) +
      coord_cartesian(xlim = c(dmy("01-Feb-2021"), plotting_data$fs_horizon)) +
  
      scale_y_continuous(breaks = scales::breaks_extended(3),
                         labels = scales::label_comma(),
                         position = 'right') +
      
      theme_minimal() +
      theme(legend.position = 'bottom') +
      guides(fill = 'none') +
      plotting_data$forecast_line +
      horizon_line
    
    
    cowplot::plot_grid(p1,p2,p3,p4,
                       ncol = 1,
                       align = 'v', axis = 'lr') %>%
      print()
  }
  
  
  
  dev.off()

}
