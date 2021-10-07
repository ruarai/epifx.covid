
find_m <- function(R_target, transition_matrix, stable_age = NULL) {
  # this function from STEPS
  #https://github.com/steps-dev/steps/blob/74c5359dd4470c4056cd799c53ef56d503ba69da/R/growth_transition_functions-class.R#L266
  #
  # compute the m that calibrates a next generation matrix to R0
  
  obj <- function (m, R_target, transition_matrix, stable_age = NULL) {
    new_transition_matrix <- m*transition_matrix
    R_current <- get_R(new_transition_matrix, stable_age = stable_age)
    (R_current - R_target) ^ 2
  } 
  
  out <- stats::optimise(f = obj,
                         interval = c(0, 1),
                         R_target,
                         transition_matrix,
                         stable_age)
  out$minimum
  
  return(out$minimum)
}
get_R <- function (transition_matrix, stable_age = NULL, tolerance = 0.001, max_iter = 1000) {
  #Re(eigen(x)$value[1])
  # function from STEPS
  # https://github.com/steps-dev/steps/blob/74c5359dd4470c4056cd799c53ef56d503ba69da/R/growth_transition_functions-class.R#L211
  # compute R from a transition (next generation) matrix
  
  if (is.null(stable_age)) {
    stable_age <- rep(1, ncol(transition_matrix))
  }
  old_stages <- stable_age
  converged <- FALSE
  iter <- 0
  old_Rs <- rep(.Machine$double.eps, ncol(transition_matrix))
  
  while (!converged & iter < max_iter) {
    new_stages <- transition_matrix %*% old_stages
    Rs <- new_stages / old_stages
    errors <- abs(1 - (Rs / old_Rs))
    converged <- all(errors < tolerance)
    old_Rs <- Rs
    old_stages <- new_stages
    iter <- iter + 1
  }
  
  if (!converged) {
    warning(
      "estimation of growth rate did not converge in ",
      max_iter,
      " iterations"
    )
  }
  
  # return the intrinsic growth rate
  Rs[1]
  
}

baseline_matrix <- function(R0 = 2.5, final_age_bin = 80) {
  # construct a next generation matrix for Australia from Prem matrix
  
  
  # Prem 2017 contact matrix
  contact_matrix_raw <- readxl::read_xlsx(
    path = "vaccination/data/MUestimates_all_locations_1.xlsx",
    sheet = "Australia",
    col_types = rep("numeric", 16)
  ) %>%
    as.matrix
  
  # expand out to add an 80+ category the same as the 75-80 category
  contact_matrix <- matrix(NA, 17, 17)
  contact_matrix[17, 17] <- contact_matrix_raw[16, 16]
  contact_matrix[17, 1:16] <- contact_matrix_raw[16, ]
  contact_matrix[1:16, 17] <- contact_matrix_raw[, 16]
  contact_matrix[1:16, 1:16] <- contact_matrix_raw
  
  # set names
  bin_names <- age_classes(80)$classes 
  dimnames(contact_matrix) <- list(
    bin_names,
    bin_names
  )
  
  # relative infectiousness data from Trauer et al 2021
  age_susceptability <- readr::read_csv(
    file = "vaccination/data/trauer_2021_supp_table5.csv",
    col_names = c(
      "age_group",
      "clinical_fraction",
      "relative_susceptability",
      "infection_fatality_rate",
      "proportion_symtomatic_hospitalised"
    ),
    col_types = cols(
      age_group = col_character(),
      clinical_fraction = col_double(),
      relative_susceptability = col_double(),
      infection_fatality_rate = col_double(),
      proportion_symtomatic_hospitalised = col_double()
    ),
    skip = 1
  ) %>%
    # duplicate the final row, and re-use for the 74-79 and 80+ classes
    add_row(
      .[16, ]
    ) %>%
    mutate(
      age_class = bin_names
    ) %>%
    dplyr::select(
      age_class,
      everything(),
      -age_group
    )
  
  # calculate relative infectiousness - assume asymptomatics are 50% less
  # infectious, and use age-stratified symptomaticity
  relative_infectiousness <- age_susceptability$clinical_fraction*1 + 0.5*(1 - age_susceptability$clinical_fraction)
  q <- relative_infectiousness
  q_scaled <- q/max(q)
  
  # apply the q scaling before computing m
  contact_matrix_scaled <- sweep(contact_matrix, 2, q_scaled, FUN = "*")
  
  # calculate m - number of onward infections per relative contact 
  m <- find_m(
    R_target = R0,
    transition_matrix = contact_matrix_scaled
  )
  
  contact_matrix_scaled * m
  
}