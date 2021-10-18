


factor_table <- expand_grid(
  scenario_name = "no_reversion",
  run_name = c("2021-09-10","2021-09-17", "2021-09-25", "2021-10-02")
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


date_range <- c(ensemble_samples$date %>% min() - 7,
                local_cases$date %>% max())

state_limits <- local_cases %>%
  filter(date >= max(date) - 28) %>%
  group_by(state) %>%
  summarise(max = max(count) * 1.5) %>%
  pull(max) %>%
  `names<-`(local_cases$state %>% unique())


## Plotting

require(matrixStats)
require(Cairo)

fs_date <- forecasting_dates$date_last_onset_50

filter_state <- . %>%
  filter(state %in% c("VIC", "NSW", "ACT"),
         date <= fs_date + 28)

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
  filter(date <= fs_date + 28) %>%
  
  group_by(state) %>%
  mutate(upper_90_lim = state_limits[state],
         upper_90 = if_else(upper_90 > upper_90_lim, upper_90_lim, upper_90),
         upper_50 = if_else(upper_50 > upper_90_lim, upper_90_lim, upper_50))



ggplot() +
  
  geom_linerange(aes(x = date, ymin = 0, ymax = count),
           local_cases %>% filter_state,
           size = 1,
           color = 'gray70') +
  
  geom_point(aes(x = date, y = count / detection_probability),
             local_cases %>% filter_state %>% filter(detection_probability < 0.95),
             pch = 1) +
  
  geom_ribbon(aes(x = date, ymin = lower_90, ymax = upper_90,
                  group = plot_group,
                  fill = "90"),
              alpha = 0.5,
              ensemble_quants) +
  
  
  geom_ribbon(aes(x = date, ymin = lower_50, ymax = upper_50,
                  group = plot_group,
                  fill = "50"),
              alpha = 0.5,
              ensemble_quants) +
  
  geom_line(aes(x = date, y = median, group = plot_group,
                color = "median"),
            ensemble_quants) +
  
  facet_grid(cols = vars(run_name),
             rows = vars(state),
             scales = "free_y") + #,
  #labeller = as_labeller(scenario_labels)) +
  scale_y_continuous(breaks = scales::breaks_extended(10),
                     labels = scales::label_comma(),
                     position = 'right') +
  
  scale_x_date(breaks = scales::breaks_pretty(5),
               labels = scales::label_date_short(),
               limits = date_range) +
  
  ylab("Notifications") + xlab("Date") +
  
  scale_color_brewer("Forecast",
                     palette = "Blues") +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  
  theme_minimal() +
  theme(legend.position = 'bottom') +
  guides(fill = 'none', color = 'none')

ggsave(paste0("results/", run_name, "/comparison_performance.png"),
       width = 12, height = 8, bg = 'white')
