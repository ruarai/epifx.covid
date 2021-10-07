


factor_table <- expand_grid(
  scenario_name = scenarios$name,
  run_name = c("2021-10-02")
)


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

fs_date <- last_onset_date

filter_state <- . %>%
  filter(state %in% c("VIC", "NSW", "ACT"))

data_cols <- ensemble_samples %>%
  select(-starts_with("sim"))

ensemble_quants <- ensemble_samples %>% select(starts_with("sim")) %>%
  as.matrix() %>%
  rowQuantiles(probs = c(0.05,0.25, 0.5, 0.75,0.95)) %>%
  bind_cols(data_cols, as_tibble(.)) %>%
  select(-1) %>%
  `colnames<-`(c(names(data_cols),"lower_90", "lower_50", "median", "upper_50", "upper_90")) %>%
  mutate(plot_group = interaction(scenario_name, run_name),
         upper_90 = ifelse(upper_90 > 2000, 2000, upper_90)) %>%
  filter_state %>%
  filter(date >= fs_date, date <= fs_date + 28)




CairoPDF(paste0("results/", run_name, "/comparison_reversion.pdf"), 
         width = 8, height = 6)

# scenario_labels <- c("no_reversion" = "No reversion to C1",
#                      "with_reversion" = "With reversion to C1")

ggplot() +
  
  
  geom_ribbon(aes(x = date, ymin = lower_90, ymax = upper_90,
                  group = scenario_name,
                  fill = scenario_name),
              alpha = 0.2,
              ensemble_quants) +
  
  
  geom_ribbon(aes(x = date, ymin = lower_50, ymax = upper_50,
                  group = scenario_name,
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
  scale_y_continuous(breaks = scales::breaks_extended(10),
                     labels = scales::label_comma(),
                     position = 'right') +
  
  scale_x_date(breaks = scales::breaks_pretty(5),
               labels = scales::label_date_short()) +
  coord_cartesian(xlim = c(fs_date, fs_date + 28)) +
  
  ylab("Notifications") + xlab("Date") +
  
  scale_color_brewer("Scenario",
                     palette = 1, type = 'qual') +
  scale_fill_brewer(palette = 1, type = 'qual') +
  
  theme_minimal() +
  theme(legend.position = 'bottom') +
  guides(fill = 'none')

dev.off()
