


factor_table <- expand_grid(
  scenario_name = scenarios$name,
  run_name = run_name
)


local_cases <- read_csv("data/in/local_cases_input.csv") %>%
  rename(date = date_onset) %>%
  filter(detection_probability > 0.5)

ensemble_spec <- "cDdddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd"

ensemble_samples <- pmap_dfr(
  factor_table,
  function(scenario_name, run_name){
    file_name <- str_c("ensemble_samples/",
                       "moss_forecast_samples_vacc_",
                       scenario_name, 
                       "_", run_name, ".csv")
    
    vroom::vroom(file_name, col_types = ensemble_spec) %>%
      mutate(scenario_name = !!scenario_name,
             run_name = !!run_name)
    })


## Plotting

require(matrixStats)
require(Cairo)

fs_date <- forecasting_dates$date_last_onset_50

filter_state <- . %>%
  filter(state %in% c("VIC", "NSW", "ACT"),
         date >= fs_date - 14, date <= fs_date + 28)

data_cols <- ensemble_samples %>%
  select(-starts_with("sim"))

ensemble_quants <- ensemble_samples %>% select(starts_with("sim")) %>%
  as.matrix() %>%
  rowQuantiles(probs = c(0.05,0.25, 0.5, 0.75,0.95)) %>%
  bind_cols(data_cols, as_tibble(.)) %>%
  select(-1) %>%
  `colnames<-`(c(names(data_cols),"lower_90", "lower_50", "median", "upper_50", "upper_90")) %>%
  mutate(plot_group = interaction(scenario_name, run_name)) %>%
  filter_state %>%
  filter(date >= fs_date - 14, date <= fs_date + 28) %>%
  
  group_by(state) # %>%
  # mutate(upper_90_lim = max(upper_90) * 0.5,
  #        upper_90 = if_else(upper_90 > upper_90_lim, upper_90_lim, upper_90))



ggplot() +


  geom_ribbon(aes(x = date, ymin = lower_90, ymax = upper_90,
                  group = plot_group,
                  fill = scenario_name),
              alpha = 0.2,
              ensemble_quants) +

  
  geom_ribbon(aes(x = date, ymin = lower_50, ymax = upper_50,
                  group = plot_group,
                  fill = scenario_name),
              alpha = 0.5,
              ensemble_quants) +
  
  geom_line(aes(x = date, y = median, group = plot_group,
                color = scenario_name),
            ensemble_quants) +
  
  facet_grid(cols = vars(scenario_name),
             rows = vars(state),
             scales = "free_y") + #,
             #labeller = as_labeller(scenario_labels)) +
  scale_y_continuous(breaks = scales::breaks_extended(5),
                     labels = scales::label_comma(),
                     position = 'right') +
  
  scale_x_date(breaks = scales::breaks_pretty(5),
               labels = scales::label_date_short()) +
  
  ylab("Notifications") + xlab("Date") +
  
  scale_color_brewer("Forecast",
                     palette = 2, type = 'qual') +
  scale_fill_brewer(palette = 2, type = 'qual') +
  
  theme_minimal() +
  theme(legend.position = 'none') +
  guides(fill = 'none')

ggsave(paste0("results/", run_name, "/scenario_comparison.png"),
       width = 8, height = 6, bg = 'white')



ggplot() +


  geom_ribbon(aes(x = date, ymin = lower_50, ymax = upper_50,
                  group = plot_group,
                  fill = scenario_name),
              alpha = 0.5,
              ensemble_quants) +
  
  geom_line(aes(x = date, y = median, group = plot_group,
                color = scenario_name),
            ensemble_quants) +
  
  facet_grid(cols = vars(scenario_name),
             rows = vars(state),
             scales = "free_y") + #,
             #labeller = as_labeller(scenario_labels)) +
  scale_y_continuous(breaks = scales::breaks_extended(5),
                     labels = scales::label_comma(),
                     position = 'right') +
  
  scale_x_date(breaks = scales::breaks_pretty(5),
               labels = scales::label_date_short()) +
  
  ylab("Notifications") + xlab("Date") +
  
  scale_color_brewer("Forecast",
                     palette = 2, type = 'qual') +
  scale_fill_brewer(palette = 2, type = 'qual') +
  
  theme_minimal() +
  theme(legend.position = 'none') +
  guides(fill = 'none')

ggsave(paste0("results/", run_name, "/scenario_comparison_50.png"),
       width = 8, height = 6, bg = 'white')
