
## Initialize the workspace
setwd("~/reff-covid-epifx-aus-2020-dev")
rm(list = ls())
library(tidyverse)
library(lubridate)
library(future.callr)
library(listenv)
library(furrr)


### Downloading latest data
## This includes the latest local_cases_input, r_eff_12_local_samples, external exposures, 
## as well as the vaccination time-series from the Curtain model

source("interface_functions/dropbox.R")
source("interface_functions/data.R")

dbx_files <- tribble(
  ~remote_file, ~local_file,
  "/covid_output/local_cases_input.csv",                 "data/in/local_cases_input.csv",
  "/covid_output/projection/r_eff_12_local_samples.csv", "data/in/r_eff_proj_samples.csv",
  "/covid_output/epifx_in/external_exposures_all.csv",   "data/in/external_exposures_all.csv",
)

# The vaccine files have different names for each week, so try and work out what they are:
dbx_files <- bind_rows(dbx_files, get_vaccine_files())

# And then download all of them
download_files(dbx_files)



# Create a 'run_name' that we use to label our results
# This is the date of the latest case data date with pr_detect > 0.5
forecasting_dates <- get_forecast_dates("data/in/local_cases_input.csv")

run_name <- forecasting_dates$date_last_onset_50 %>% format("%Y-%m-%d")

run_name <- str_c(run_name, "_proj_vacc")

### Processing vaccine data

source("interface_functions/vaccination.R")


produce_vaccination_inputs("data/in/effective_dose_data.csv",
                           "data/in/vaccine_effect_timeseries.csv",
                           run_name)


### Processing R_eff, case counts:


# Read in the local cases file and split into SSV files for each state
process_local_cases("data/in/local_cases_input.csv", 
                    pr_detect_threshold = 0.5)

# As above, for our external exposures listing
process_external_exposures("data/in/external_exposures_all.csv")


### Defining our different scenarios for ensembling:

args_with_reversion <- tribble(~key,          ~value,
                               "%exps_dir%",  "exps/with_reversion",
                               "%reffs_dir%", "with_reversion")
args_no_reversion <- tribble(~key,          ~value,
                             "%exps_dir%",  "exps/no_reversion",
                             "%reffs_dir%", "no_reversion")

scenarios <- tribble(
  ~name,            ~args,                     ~reff_truncation_date,
  "with_reversion", list(args_with_reversion), ymd(NA),
  "no_reversion",   list(args_no_reversion),   forecasting_dates$date_last_infection_50,
) %>%
  rowwise() %>%
  mutate(exps_dir = args[[1]] %>% filter(key == "%exps_dir%") %>% pull(value))



# Produce C1/2 trajectories, with truncation at reff_truncation_date, if that is not NULL
for(i_scen in 1:nrow(scenarios)) {
  scenario <- scenarios[i_scen,]
  dir.create(paste0("data/",scenario$name), showWarnings = FALSE)
  
  process_reff_trajs(
    reff_proj_file = "data/in/r_eff_proj_samples.csv",
    vaccine_effect_file = "data/in/vaccine_effect_timeseries.csv",
    output_path = paste0("data/",scenario$name),
    truncation_date = scenario$reff_truncation_date)
}

source("interface_functions/epifx.R")

run_forecast <- function(name, args, reff_truncation_date, exps_dir) {
  print(paste0("Running scenario forecast ", name))
  
  dir.create(exps_dir, recursive = TRUE, showWarnings = FALSE)
  file.remove(list.files(exps_dir, full.names = TRUE))
  
  scen_template_file <- paste0("epifx_params/forecast_", name, ".toml")
  
  # Create a custom TOML file for this scenario:
  create_forecast_file("epifx_params/forecast_template.toml",
                       args[[1]],
                       scen_template_file)
  
  system(paste0(venv_prefix, "python3 run_forecast.py ", "--ff ", scen_template_file),
                intern = TRUE)
}



### Running the forecast across our scenarios:
plan(callr)
future_pmap(scenarios, run_forecast)



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

for(i_scen in 1:nrow(scenarios)) {
  scenario <- scenarios[i_scen,]
  
  plot_subdir <- paste0("results/", run_name, "/", scenario$name)
  
  # Load in basic information about exps in forecast_exps_dir
  plotting_data <- get_plotting_data(scenario$exps_dir)
  
  dir.create(plot_subdir, recursive = TRUE, showWarnings = FALSE)
  
  job_list[[length(job_list) + 1]] %<-% { create_ensemble_samples(scenario$name, run_name, plotting_data) }
  job_list[[length(job_list) + 1]] %<-% { produce_CI_plots(plotting_data, plot_subdir) }
  job_list[[length(job_list) + 1]] %<-% { produce_state_summaries(plotting_data, plot_subdir) }
  job_list[[length(job_list) + 1]] %<-% { plot_vaccine_effect(plotting_data, plot_subdir,
                                                              "data/in/vaccine_effect_timeseries.csv") }
}

# Run every job in parallel now
job_list <- as.list(job_list)


### Upload ensemble samples

ensembles_for_upload <- str_c("ensemble_samples/",
                              "moss_forecast_samples_vacc_",
                              "no_reversion",
                              "_", run_name, ".csv")


for(i_file in ensembles_for_upload) {
  print(paste0("Uploading ", i_file))
  system(paste0("bash upload.sh ", i_file))
}



### Upload plots

source("interface_functions/dropbox.R")
upload_plot_files(paste0("results/", run_name))


### Backup our exps directory in case we ever care about that!

source("interface_functions/spartan_backup.R")

backup_exps_files(run_name)
