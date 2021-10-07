

state_populations <- read_csv("vaccination/data/demographics/population_by_state.csv") %>%
  filter(!is.na(state))

age_distribution_by_state <- read_csv("vaccination/data/demographics/age_distribution_by_state.csv")
adult_population_by_state <- read_csv("vaccination/data/demographics/adult_population_by_state.csv")

age_lookup <- tibble::tribble(
  ~age_5y, ~age, ~proportion_of_group,
  "0-4", "0-14", 5/15,
  "5-9", "0-14", 5/15,   
  "10-14", "0-14", 5/15,
  "15-19", "15-29", 5/15,
  "20-24", "15-29", 5/15,
  "25-29", "15-29", 5/15,
  "30-34", "30-39", 5/10,
  "35-39", "30-39", 5/10,
  "40-44", "40-49", 5/10,
  "45-49", "40-49", 5/10,
  "50-54", "50-59", 5/10,
  "55-59", "50-59", 5/10,
  "60-64", "60-69", 5/10,
  "65-69", "60-69", 5/10,
  "70-74", "70-79", 5/10,
  "75-79", "70-79", 5/10,
  "80+", "80+", 1/1
)



age_classes <- function(final_age_bin = 80, by = 5) {
  
  # compute age classes based on this spec
  ages_lower = seq(0, final_age_bin, by = by)
  n_ages <- length(ages_lower)
  ages_upper = ages_lower + by - 1
  ages_upper[n_ages] <- Inf
  
  age_classes <- c(
    paste(
      ages_lower[-n_ages],
      ages_upper[-n_ages],
      sep = "-"
    ),
    paste0(
      final_age_bin,
      "+"
    )
  )
  
  tibble::tibble(
    classes = age_classes,
    lower = ages_lower,
    upper = ages_upper
  )
  
}
