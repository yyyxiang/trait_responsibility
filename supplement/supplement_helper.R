load_fit_with_loo <- function(fit_file, loo_file) {
  fit <- readRDS(fit_file)
  fit <- add_criterion(fit, criterion = "loo", file = loo_file)
  fit
}

fig_labels <- c(
  actual_contribution = "Actual contribution",
  focal_mutation_counterfactual = "Focal trait\n mutation counterfactual",
  nonfocal_mutation_counterfactual = "Non-focal trait\n mutation counterfactual",
  focal_population_counterfactual = "Focal trait\n population counterfactual",
  nonfocal_population_counterfactual = "Non-focal trait\n population counterfactual",
  ensemble = "Ensemble",
  focal_perceived_population_counterfactual = "Focal trait\n perceived population counterfactual",
  nonfocal_perceived_population_counterfactual = "Non-focal trait\n perceived population counterfactual"
)
