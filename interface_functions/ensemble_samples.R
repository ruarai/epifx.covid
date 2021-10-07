

create_ensemble_samples <- function(scenario_name,
                                    run_name,
                                    plotting_data){
  
  ensemble_data <- plotting_data$hdf_files %>%
    pmap(function(file, state) {
      h5read(file = file,
             name = "/data/sim_obs") %>%
        select(date, value) %>%
        mutate(date = read_hdf_date(date),
               state = !!state,
               ix = rep(1:2000, length.out = nrow(.))) %>%
        
        filter(date >= plotting_data$fs_date) %>%
        
        pivot_wider(names_from = ix,
                    names_prefix = "sim",
                    values_from = value,
                    values_fn = first)
    })
  
  ensemble_data <- ensemble_data %>%
    bind_rows() %>%
    relocate(state, date)
  
  
  ensemble_file_name <- str_c("ensemble_samples/",
                              "moss_forecast_samples_",
                              "vacc_",
                              scenario_name, "_",
                              run_name,
                              ".csv")
  
  vroom::vroom_write(ensemble_data, ensemble_file_name, delim = ",")
}

create_mixture_file <- function(scenarios,
                                run_name) {
  
  
  ensembles_for_mixing <- str_c("ensemble_samples/",
                                "moss_forecast_samples_vacc_",
                                scenarios$name, 
                                "_", run_name, ".csv")
  
  ensemble_spec <- "cDdddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd"
  
  ensemble_samples_A <- vroom::vroom(ensembles_for_mixing[1], col_types = ensemble_spec)
  ensemble_samples_B <- vroom::vroom(ensembles_for_mixing[2], col_types = ensemble_spec)
  
  
  stopifnot(ensemble_samples_A$state == ensemble_samples_B$state)
  stopifnot(ensemble_samples_A$date == ensemble_samples_B$date)
  
  data_cols <- ensemble_samples_A %>%
    select(date, state)
  
  
  rand_cols_A <- sample(1:2000, 1000)
  rand_cols_B <- sample(1:2000, 1000)
  
  ensemble_mix <-
    bind_cols(
      data_cols,
      ensemble_samples_A %>% select(num_range("sim", rand_cols_A)),
      ensemble_samples_B %>% select(num_range("sim", rand_cols_B))
    )
  colnames(ensemble_mix) <- c("date", "state", str_c("sim", 1:2000))
  
  
  ensemble_file_name <- str_c("ensemble_samples/",
                              "moss_forecast_samples_",
                              "vacc_",
                              "reversion_mix_",
                              run_name,
                              ".csv")
  
  vroom::vroom_write(ensemble_mix, ensemble_file_name, delim = ",")
  
  
  
  
  
  
  ## Plotting
  
  require(matrixStats)
  require(Cairo)
  
  filter_state <- . %>%
    filter(state %in% c("ACT", "VIC", "NSW"))
  
  ensemble_quants_A <- ensemble_samples_A %>% select(starts_with("sim")) %>%
    as.matrix() %>%
    rowQuantiles(probs = c(0.25, 0.5, 0.75)) %>%
    bind_cols(data_cols, as_tibble(.)) %>%
    select(-1) %>%
    `colnames<-`(c("date", "state", "lower_50", "median", "upper_50")) %>%
    filter_state
  
  ensemble_quants_B <- ensemble_samples_B %>% select(starts_with("sim")) %>%
    as.matrix() %>%
    rowQuantiles(probs = c(0.25, 0.5, 0.75)) %>%
    bind_cols(data_cols, as_tibble(.)) %>%
    select(-1) %>%
    `colnames<-`(c("date", "state", "lower_50", "median", "upper_50")) %>%
    filter_state
  
  mix_quants <- ensemble_mix %>% select(starts_with("sim")) %>%
    as.matrix() %>%
    rowQuantiles(probs = c(0.25, 0.5, 0.75)) %>%
    bind_cols(data_cols, as_tibble(.)) %>%
    select(-1) %>%
    `colnames<-`(c("date", "state", "lower_50", "median", "upper_50")) %>%
    filter_state
  
  
  
  CairoPNG(paste0("results/", run_name, "/mix_comparison.png"), 
           width = 12, height = 6, dpi = 200, units = "in")
  
  p <- ggplot() +
    
    
    geom_ribbon(aes(x = date, ymin = lower_50, ymax = upper_50, fill = 'A'),
                alpha = 0.3,
                ensemble_quants_A) +
    geom_ribbon(aes(x = date, ymin = lower_50, ymax = upper_50, fill = 'B'),
                alpha = 0.3,
                ensemble_quants_B) +
    geom_ribbon(aes(x = date, ymin = lower_50, ymax = upper_50, fill = 'C'),
                alpha = 0.3,
                mix_quants) +
    
    geom_line(aes(x = date, y = median, color = 'C'),
              mix_quants) +
    geom_line(aes(x = date, y = median, color = 'A'),
              ensemble_quants_A) +
    geom_line(aes(x = date, y = median, color = 'B'),
              ensemble_quants_B) +
    
    facet_wrap(~state, scales = "free_y", nrow = 1) +
    scale_color_brewer("",
                       type = 'qual', palette = 2,
                       labels = c("A" = scenarios$name[1],
                                  "B" = scenarios$name[2],
                                  "C" = "Mix for ensemble")) +
    scale_y_continuous(breaks = scales::breaks_extended(10),
                       labels = scales::label_comma(),
                       position = 'right') +
    
    scale_x_date(breaks = scales::breaks_pretty(5),
                 labels = scales::label_date_short()) +
    coord_cartesian(xlim = c(fs_date, fs_date + 28)) +
    
    ylab("Notifications") + xlab("Date") +
    
    scale_fill_brewer(type = 'qual', palette = 2) +
    
    theme_minimal() +
    theme(legend.position = 'bottom') +
    guides(fill = 'none')
  
  print(p)
  
  dev.off()
}




