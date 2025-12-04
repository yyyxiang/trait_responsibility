sigmoid <- function(x) 1 / (1 + exp(-x))

responsibility_focal <- function(df, cf_fun, 
                                 all_traits = c('openness', 'conscientiousness', 'extraversion', 'agreeableness', 'neuroticism'), 
                                 all_scores = c(0, 20, 40, 60, 80, 100)) {
  out <- df
  for (tr in all_traits) {
    observed <- df[[tr]]
    others <- setdiff(all_traits, tr)
    
    # sum_{j!=i} w_j * X_j
    non_focal_contribution <- Reduce(`+`, lapply(others, function(t) {
      df[[paste0('weight_', t)]] * df[[t]]
    }))
    
    out[[paste0('responsibility_', tr)]] <- sapply(seq_along(observed), function(i) {
      cf_tbl <- cf_fun(tr, observed[i])
      
      # 1 - SUM_k P(a=1 | X_i' = cf_k, X_/i fixed) * P(cf_k)
      1 - sum(
        sigmoid(df[i, 'intercept'] + non_focal_contribution[i] + 
                  df[i, paste0('weight_', tr)] * cf_tbl$counterfactual) * cf_tbl$p
        )
    })
  }
  out
}

responsibility_nonfocal <- function(df, cf_fun, 
                                    all_traits = c('openness', 'conscientiousness', 'extraversion', 'agreeableness', 'neuroticism'), 
                                    all_scores = c(0, 20, 40, 60, 80, 100)) {
  out <- df
  for (tr in all_traits) {
    observed <- df[[tr]]
    others <- setdiff(all_traits, tr)
    
    # w_i * X_i
    focal_contribution <- df[[paste0('weight_', tr)]] * df[[tr]]
    
    out[[paste0('responsibility_', tr)]] <- sapply(seq_along(observed), function(i) {
      
      cf_list <- lapply(others, function(t) {
        cf_fun(t, df[i, t]) %>% 
          select(counterfactual, p) %>% 
          setNames(c(paste0('cf_', t), paste0('p_', t)))
      })
      
      # enumerate combinations of non-focal counterfactuals
      grid <- Reduce(function(a, b) tidyr::crossing(a, b), cf_list)
      
      # joint probability across non-focal traits
      joint_p <- Reduce(`*`, lapply(setdiff(all_traits, tr), function(t) grid[[paste0('p_', t)]]))
      
      # sum_{j!=i} w_j * cf_j
      non_focal_contribution <- Reduce(`+`, lapply(others, function(t) {
        df[i, paste0('weight_', t)] * grid[[paste0('cf_', t)]]
      }))
      
      # SUM_k P(a=1 | X_i fixed, X_/i' = combo_k) * P(combo_k)
      sum(sigmoid(df[i, 'intercept'] + focal_contribution[i] + non_focal_contribution) * joint_p)
    })
  }
  out
}

# mutation-based: p \prop 1/|cf - observed|
make_cf_mutation <- function(all_scores) {
  force(all_scores)
  function(tr, observed) {
    tibble(counterfactual = all_scores[all_scores != observed]) %>% 
      mutate(p = 1 / abs(counterfactual - observed),
             p = p / sum(p))
  }
}

# population-based: look up from pop_df
make_cf_from_population <- function(pop_df) {
  force(pop_df)
  function(tr, observed) {
    pop_df %>% 
      filter(trait == tr, observed_score == observed) %>% 
      select(counterfactual, p)
  }
}

rank_loo <- function(loo_cmp) as.data.frame(loo_cmp) %>% rownames_to_column("model") %>% as_tibble()

boot_mean <- function(x, idx) mean(x[idx], na.rm = TRUE)

fit_with_loo <- function(formula, data, file_stub, seed = 1,
                         reloo = FALSE,
                         preset = c("fast","default","careful")) {
  preset <- match.arg(preset)
  
  # preset settings
  cfg <- switch(
    preset,
    "fast" = list(iter = 1000, warmup = 500, chains = 2,
                  adapt_delta = 0.8, max_treedepth = 10),
    "default" = list(iter = 2000, warmup = 1000, chains = 4,
                     adapt_delta = 0.9, max_treedepth = 10),
    "careful" = list(iter = 4000, warmup = 1000, chains = 4,
                     adapt_delta = 0.95, max_treedepth = 12)
  )
  
  fit <- brm(
    formula = formula,
    data = data,
    file = glue::glue("cache/{file_stub}"),
    seed = seed,
    cores = max(1, parallel::detectCores() - 1),
    iter = cfg$iter, warmup = cfg$warmup, chains = cfg$chains,
    control = list(adapt_delta = cfg$adapt_delta, max_treedepth = cfg$max_treedepth)
  )
  
  fit <- add_criterion(
    fit, criterion = "loo", reloo = reloo,
    file = glue::glue("cache/{file_stub}_loo")
  )
  fit
}

fig_labels <- c(
  production = "Production",
  focal_mutation_counterfactual = "Focal trait\n mutation counterfactual",
  nonfocal_mutation_counterfactual = "Non-focal trait\n mutation counterfactual",
  focal_population_counterfactual = "Focal trait\n population counterfactual",
  nonfocal_population_counterfactual = "Non-focal trait\n population counterfactual",
  ensemble = "Ensemble"
)

rank_loo <- function(loo_cmp) as.data.frame(loo_cmp) %>% rownames_to_column("model") %>% as_tibble()

# plot ELPD diff
loo_plot <- function(loo_cmp,
                     model_levels = NULL) {
  
  df <- as.data.frame(loo_cmp) %>%
    tibble::rownames_to_column("model") %>% 
    mutate(model = factor(model, levels = rev(model_levels), ordered = TRUE))
  
  ggplot(df, aes(x = model, y = elpd_diff)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_pointrange(
      aes(ymin = elpd_diff - 1.96*se_diff, ymax = elpd_diff + 1.96*se_diff),
      color = 'black', shape = 16, size = 0.6
    ) +
    coord_flip() +
    labs(x = NULL, y = expression(Delta~ELPD)) +
    scale_x_discrete(labels = fig_labels) +
    theme_classic(base_size = 12)
}

scatter_plot <- function(data_model) {
  ggplot(data_model, aes(x = predicted_responsibility, y = responsibility, color = focal_trait)) +
    geom_point(size = 0.5) +
    geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0, alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
    facet_wrap(~model, ncol = 3,
               labeller = as_labeller(fig_labels)) +
    theme_classic() +
    scale_x_continuous(limits = c(0,100), breaks = seq(0,100,25)) +
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100,25)) +
    coord_equal() +
    geom_text(
      data = corr_labels,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1, vjust = 0,
      fontface = "italic",
      size = 3
    ) +
    theme(legend.position = "bottom") + 
    labs(x = "Model predictions", y = "Data", color = "Focal trait")
}

plot_contribution <- function(data, model_fit,
                              xvar = c("wi_Xi","sum_wjXj"),
                              nbins_x = 5, nbins_cond = 5,
                              ncol = 4) {
  xvar <- match.arg(xvar)
  condvar <- setdiff(c("wi_Xi","sum_wjXj"), xvar)
  
  panels <- c("production",
              "focal_mutation_counterfactual","nonfocal_mutation_counterfactual",
              "Data",
              "focal_population_counterfactual","nonfocal_population_counterfactual",
              "ensemble")
  
  # breaks + rounded-midpoint labels
  br_x    <- seq(min(data[[xvar]],   na.rm=TRUE), max(data[[xvar]],   na.rm=TRUE), length.out = nbins_x    + 1)
  br_cond <- seq(min(data[[condvar]],na.rm=TRUE), max(data[[condvar]],na.rm=TRUE), length.out = nbins_cond + 1)
  labs_x    <- round((head(br_x,-1)    + tail(br_x,-1))/2, 1)
  labs_cond <- round((head(br_cond,-1) + tail(br_cond,-1))/2, 1)
  
  tmp_data <- data %>%
    mutate(bin_x    = cut(.data[[xvar]],    breaks = br_x,    include.lowest = TRUE, labels = labs_x),
           bin_cond = cut(.data[[condvar]], breaks = br_cond, include.lowest = TRUE, labels = labs_cond)) %>%
    group_by(bin_cond, bin_x) %>%
    summarise(value = mean(responsibility, na.rm = TRUE), .groups = "drop") %>%
    mutate(panel = "Data")
  
  tmp_models <- model_fit %>%
    filter(model %in% panels[panels != "Data"]) %>% # exclude data from panel, keep only models
    mutate(bin_x    = cut(.data[[xvar]],    breaks = br_x,    include.lowest = TRUE, labels = labs_x),
           bin_cond = cut(.data[[condvar]], breaks = br_cond, include.lowest = TRUE, labels = labs_cond)) %>%
    group_by(model, bin_cond, bin_x) %>%
    summarise(value = mean(predicted_responsibility, na.rm = TRUE), .groups = "drop") %>%
    rename(panel = model)
  
  tmp_all <- bind_rows(tmp_data, tmp_models) %>%
    mutate(panel = factor(panel, levels = panels, ordered = TRUE))
  
  # axis and legend labels
  x_lab <- if (xvar == "wi_Xi")
    expression(paste(w[i]*X[i], " (Focal contribution)"))
  else
    expression(paste(sum(w[j]*X[j]), " (Non-focal contribution)"))
  col_lab <- if (condvar == "wi_Xi")
    expression(paste(w[i]*X[i], " (Focal contribution)"))
  else
    expression(paste(sum(w[j]*X[j]), " (Non-focal contribution)"))
  
  ggplot(tmp_all, aes(x = bin_x, y = value, group = bin_cond, color = bin_cond)) +
    geom_line() +
    facet_wrap(~ panel, ncol = ncol, labeller = as_labeller(c(fig_labels, "Data" = "Data"))) +
    theme_classic() +
    labs(x = x_lab, y = "Responsibility", color = col_lab)
}
