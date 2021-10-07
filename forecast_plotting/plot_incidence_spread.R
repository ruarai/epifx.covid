
library(tidyverse)
library(lubridate)
library(rhdf5)
library(Cairo)


#setwd("/data/gpfs/projects/punim1541/forecast_plotting/")

forecast_exps_dir <- "exps"
output_name <- "2021-09-03-improved-Ei-Et-reversion"

source("forecast_plotting/plot_incidence.R")