library(boot)
library(tidyverse)
library(brms)
library(glue)
library(patchwork)
library(BayesFactor)
source("helper.R")

dat <- read.csv("./../data/exp2.csv", header = T, na.strings = c("", "NA", "null")) %>% 
  filter(total_failed_attention_checks == 0) %>% 
  select(-total_failed_attention_checks) 

all_traits <- c("openness", "conscientiousness", "extraversion", "agreeableness", "neuroticism")
all_scores <- c(0, 20, 40, 60, 80, 100)

weights_df <- read.csv("output/trait_weights.csv", header = T)
dat_weight <- dat %>% 
  select(-starts_with("responsibility_")) %>% 
  left_join(weights_df, by = "scenario")

########## Model 1: production ##########
production_predictions <- dat_weight
for (tr in all_traits) {
  production_predictions[[paste0("responsibility_", tr)]] <- sigmoid(production_predictions[[paste0("weight_", tr)]] * production_predictions[[tr]])
}

production_predictions <- production_predictions %>% 
  select(-starts_with("weight")) %>% 
  mutate(model = "production")

########## Models 2-3: mutation-based counterfactual ##########
cf_mutation <- make_cf_mutation(all_scores)

# focal
focal_mutation_predictions <- responsibility_focal(dat_weight, cf_mutation) %>%
  select(-starts_with("weight")) %>% 
  mutate(model = "focal_mutation_counterfactual")

# non-focal
nonfocal_mutation_predictions <- responsibility_nonfocal(dat_weight, cf_mutation) %>% 
  select(-starts_with("weight")) %>% 
  mutate(model = "nonfocal_mutation_counterfactual")

########## Models 4-5: population-based counterfactual ##########
# dataset: https://github.com/automoto/big-five-data

pop_dat <- read.csv("./../data/population_data.csv", header = T, na.strings = c("", "NA", "null")) %>% 
  dplyr::rename(agreeableness_score = agreeable_score)
bin <- 10
pop_df <- NULL
for (t in all_traits) {
  sub_df <- pop_dat %>% 
    select(paste0(t, "_score")) %>% 
    setNames("score") %>% 
    mutate(score = score * 100) %>% 
    mutate(score = case_when(score < bin ~ 0,
                             score >= (20 - bin) & score < (20 + bin) ~ 20,
                             score >= (40 - bin) & score < (40 + bin) ~ 40,
                             score >= (60 - bin) & score < (60 + bin) ~ 60,
                             score >= (80 - bin) & score < (80 + bin) ~ 80,
                             score >= 100-bin ~ 100)) %>% 
    filter(!is.na(score)) %>% 
    mutate(trait = t)
  
  pop_df <- rbind(pop_df, sub_df)
}
pop_df <- pop_df %>% 
  count(trait, score) %>% 
  group_by(trait) %>%
  mutate(p = n / sum(n)) %>% 
  ungroup() %>% 
  select(-n)

pop_df <- pop_df %>%
  dplyr::rename(counterfactual = score) %>%
  crossing(observed_score = all_scores) %>%
  complete(trait, observed_score = all_scores, counterfactual = all_scores,
           fill = list(p = 0)) %>% 
  filter(observed_score != counterfactual) %>%
  select(trait, observed_score, counterfactual, p) %>%
  group_by(trait, observed_score) %>%
  mutate(p = if (sum(p) > 0) p / sum(p) else 1 / n()) %>%  # normalize within trait × observed_score; uniform if bucket empty
  ungroup() %>% 
  arrange(trait, observed_score, counterfactual)

cf_pop <- make_cf_from_population(pop_df)

# focal
focal_population_predictions <- responsibility_focal(dat_weight, cf_pop) %>% 
  select(-starts_with("weight")) %>% 
  mutate(model = "focal_population_counterfactual")

# non-focal
nonfocal_population_predictions <- responsibility_nonfocal(dat_weight, cf_pop) %>% 
  select(-starts_with("weight")) %>% 
  mutate(model = "nonfocal_population_counterfactual")

########## combine all models ##########
model_predictions <- rbind(production_predictions,
                           focal_mutation_predictions,
                           nonfocal_mutation_predictions,
                           focal_population_predictions,
                           nonfocal_population_predictions)

# local save:
if (!file.exists("output/exp2_raw_model_predictions.csv")) {
  write.csv(model_predictions,"output/exp2_raw_model_predictions.csv", row.names = FALSE)
}

########## linear transformation ##########
model_predictions <- read.csv("output/exp2_raw_model_predictions.csv", header = T, na.strings = c("", "NA", "null"))

all_models <- c(
  "production",
  "focal_mutation_counterfactual","nonfocal_mutation_counterfactual",
  "focal_population_counterfactual","nonfocal_population_counterfactual"
)

dir.create("cache", showWarnings = FALSE)

predictions <- model_predictions %>% 
  mutate(model = factor(model, levels = all_models, ordered = T)) %>% 
  filter(!is.na(model)) %>% 
  pivot_longer(cols = starts_with("responsibility_"), names_to = "focal_trait", values_to = "predicted_responsibility") %>% 
  mutate(focal_trait = sub("^responsibility_", "", focal_trait)) 

predictions_wide <- predictions %>% 
  pivot_wider(id_cols = c("subject", "scenario", "person", "focal_trait"), names_from = "model", values_from = "predicted_responsibility")

dat_long <- dat %>% 
  pivot_longer(cols = starts_with("responsibility_"), names_to = "focal_trait", values_to = "responsibility") %>% 
  mutate(focal_trait = sub("^responsibility_", "", focal_trait)) %>% 
  left_join(predictions_wide, by = c("subject", "scenario", "person", "focal_trait"))

# test situational account
subject_mean_resp <- dat_long %>% 
  group_by(subject) %>% 
  summarise(mean_resp = mean(responsibility))
bf <- ttestBF(subject_mean_resp$mean_resp, mu = 0)
bf
# 2.3392e+55 ±0%
set.seed(123)
chains = posterior(bf, iterations = 10000)
summary(chains)
# delta 2.5% 3.056, median 3.595, 95% 4.144

cf_models <- setdiff(all_models, "production")
dat_long <- dat_long %>% mutate(subject = as.factor(subject))

# fit regression models
fits_single <- setNames(vector("list", length(all_models)), all_models)
for (m in all_models) {
  fml <- bf(glue("responsibility ~ 1 + {m}"))
  fits_single[[m]] <- fit_with_loo(
    fml, dat_long,
    file_stub = glue("resp_single_{m}")
  )
}

########## Model 6: ensemble ##########
dat_long <- dat_long %>%
  mutate(ensemble = (production + (focal_mutation_counterfactual + nonfocal_mutation_counterfactual)/2) / 2)

fml_ensemble <- bf(responsibility ~ 1 + ensemble)
fit_ensemble <- fit_with_loo(
  fml_ensemble, dat_long,
  file_stub = "resp_ensemble"
)
fits_single[["ensemble"]] <- fit_ensemble

loo_all <- loo::loo_compare(c(
  lapply(fits_single, \(f) f$criteria$loo)
))
print(loo_all)

if (!file.exists("output/exp2_loo_comparison.csv")) {
  write.csv(rank_loo(loo_all),"output/exp2_loo_comparison.csv", row.names = FALSE)
}

all_models <- c(all_models, "ensemble")
p0 <- loo_plot(loo_all,
               model_levels = all_models)

# Figure 2: ELPD differences
pdf('./../figures/fig2_ELPD_diff.pdf', onefile = T, width = 10, height = 4)
p0
dev.off()

# push model predictions through regression models
for (m in all_models) {
  coefs <- fixef(fits_single[[m]])
  intercept <- coefs["Intercept", "Estimate"]
  slope <- coefs[m, "Estimate"]

  dat_long[[paste0(m, "_fitted")]] <- intercept + slope * dat_long[[m]]
}

if (!file.exists("output/exp2_fitted_models.csv")) {
  write.csv(dat_long,"output/exp2_fitted_models.csv", row.names = FALSE)
}

dat_long <- read.csv("output/exp2_fitted_models.csv", header = T, na.strings = c("", "NA", "null"))

# comparing production and ensemble predictions
a <- dat_long$production_fitted
b <- dat_long$ensemble_fitted
bf <- ttestBF(a, b, paired = TRUE)
bf
# 0.005101658 ±1.1%
set.seed(123)
chains = posterior(bf, iterations = 10000)
summary(chains)
# delta 2.5% -0.009224, median -0.0008023, 95% 0.008006

########## scatter plot and correlations ##########
# collapse fitted predictions per scenario × trait × model
model_summary_fitted <- dat_long %>% 
  pivot_longer(cols = ends_with("_fitted"), names_to = "model", values_to = "predicted_responsibility") %>%
  mutate(model = sub("_fitted$", "", model)) %>%
  group_by(model, scenario, focal_trait) %>%
  summarise(predicted_responsibility = mean(predicted_responsibility, na.rm = TRUE), .groups = "drop")

# empirical data summary (bootstrap CI)
set.seed(123)
dat_summary <- dat_long %>%
  # 1) within-subject means
  group_by(subject, scenario, focal_trait) %>%
  summarise(responsibility = mean(responsibility, na.rm = TRUE), .groups = "drop") %>%
  # 2) bootstrap CI across subjects
  reframe({
    b <- boot(responsibility, statistic = boot_mean, R = 2000)
    ci <- boot.ci(b, type = "perc")$percent[4:5]
    tibble(responsibility = mean(responsibility, na.rm = TRUE),
           lowerCI = ci[1], upperCI = ci[2])
  }, .by = c(scenario, focal_trait))

dat_model <- left_join(model_summary_fitted, dat_summary, by = c("scenario", "focal_trait")) %>% 
  mutate(model = factor(model, levels = all_models, ordered = T))

corr_labels <- dat_model %>%
  group_by(model) %>%
  summarise(
    r = cor(responsibility, predicted_responsibility, use = "complete.obs"),
    p = cor.test(responsibility, predicted_responsibility)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("r = ", format(round(r, 2), nsmall = 2)),
    x = 95,
    y = 5
  )

# model                                 r               p
# production                            0.87840679      3.439515e-33
# focal_mutation_counterfactual         -0.02144818     8.322556e-01
# nonfocal_mutation_counterfactual      -0.13843340     1.695806e-01
# focal_population_counterfactual       0.72692182      1.107543e-17
# nonfocal_population_counterfactual    -0.60608417     2.360006e-11
# ensemble                              0.87945914      2.306389e-33

# Figure 3: correlation scatter plots
pdf('./../figures/fig3_scatter_plot.pdf', onefile = T, width = 10, height = 6)
scatter_plot(dat_model)
dev.off()

########## focal vs. non-focal contributions ##########
# weights in long form keyed by focal_trait name
weights_long_ref <- weights_df %>%
  pivot_longer(starts_with("weight_"),
               names_to = "focal_trait", values_to = "wi",
               names_prefix = "weight_")

# compute sum_all
dat_weight$sum_all <- rowSums(sapply(all_traits, function(tr) {
  dat_weight[[paste0("weight_", tr)]] * dat_weight[[tr]]
}))

# make a long table with Xi for each focal_trait, then attach wi and build contributions
contrib_long <- dat_weight %>%
  select(subject, scenario, person, all_of(all_traits), sum_all) %>%
  pivot_longer(all_of(all_traits), names_to = "focal_trait", values_to = "Xi") %>%
  left_join(weights_long_ref, by = c("scenario", "focal_trait")) %>%
  mutate(
    wi_Xi    = wi * Xi,
    sum_wjXj = sum_all - wi_Xi
  ) %>%
  select(subject, scenario, person, focal_trait, wi_Xi, sum_wjXj) 

# add contributions to dat_long
dat_long <- dat_long %>%
  left_join(contrib_long, by = c("subject", "scenario", "person", "focal_trait"))

fitted_long <- dat_long %>%
  pivot_longer(ends_with("_fitted"), names_to = "model", values_to = "predicted_responsibility") %>%
  mutate(
    model = sub("_fitted$", "", model),
    model = factor(model, levels = all_models, ordered = TRUE)
  )

p2a <- plot_contribution(dat_long, fitted_long, xvar = "wi_Xi",    nbins_x = 5, nbins_cond = 5)
p2b <- plot_contribution(dat_long, fitted_long, xvar = "sum_wjXj", nbins_x = 5, nbins_cond = 5)

# Figure 4: focal and non-focal trait contributions
pdf('./../figures/fig4_contribution.pdf', onefile = T, width = 10, height = 8)
(p2a / p2b) + plot_annotation(tag_levels = 'A')
dev.off()
