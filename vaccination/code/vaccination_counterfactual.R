




vacc_data_vic <- read_table2("data/vacc-data-vic.ssv")
vacc_data_vic$rate <- 0
vacc_data_vic %>% write_delim("data/vic_cf_vacc/vacc-data-vic.ssv")

# Plotting functions
source("forecast_plotting/plot_parameter_CIs.R")
source("forecast_plotting/plot_parameter_CIs_by_name.R")
source("forecast_plotting/plot_incidence.R")

source("forecast_plotting/plotting_common.R")
source("forecast_plotting/produce_CI_plots.R")
source("forecast_plotting/produce_state_summaries.R")
source("forecast_plotting/plot_vaccine_effect.R")

source("interface_functions/ensemble_samples.R")


job_list <- listenv()


plot_subdir <- paste0("results/", "vic_cf_vacc")

# Load in basic information about exps in forecast_exps_dir
plotting_data <- get_plotting_data("exps/vic_cf_vacc/")
plotting_data$fs_horizon <- ymd("2021-10-01")

dir.create(plot_subdir, recursive = TRUE, showWarnings = FALSE)

job_list[[length(job_list) + 1]] %<-% { produce_CI_plots(plotting_data, plot_subdir) }
#job_list[[length(job_list) + 1]] %<-% { produce_state_summaries(plotting_data, plot_subdir) }
#job_list[[length(job_list) + 1]] %<-% { plot_vaccine_effect(plotting_data, plot_subdir,
#                                                            "data/in/vaccine_effect_timeseries.csv") }


# Run every job in parallel now
job_list <- as.list(job_list)