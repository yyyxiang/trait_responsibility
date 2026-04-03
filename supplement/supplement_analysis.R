library(boot)
library(tidyverse)
library(brms)
library(glue)
library(patchwork)
library(BayesFactor)
source("./../code/helper.R")
source("supplement_helper.R")

########## perceived vs. actual trait score distribution in population ##########
perceived_pop_dat <- read.csv("perceived_pop_pilot.csv", header = T, na.strings = c("", "NA", "null"))

all_traits <- c("openness", "conscientiousness", "extraversion", "agreeableness", "neuroticism")
all_scores <- c(0, 20, 40, 60, 80, 100)

perceived_pop_dat <- perceived_pop_dat %>% 
  pivot_longer(cols = starts_with("p_"), names_to = "score", values_to = "p") %>% 
  mutate(score = as.integer(sub("^p_", "", score))) %>% 
  mutate(type = "Perceived distribution")

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
  select(-n) %>% 
  mutate(trait = as_factor(trait)) %>% 
  complete(trait, score = all_scores, fill = list(p = 0)) %>%
  mutate(type = "Actual distribution") %>% 
  arrange(trait, score)

# Figure S1
pdf('figS1.pdf', onefile = T, width = 8, height = 5)
p1 <- perceived_pop_dat %>% 
  ggplot(aes(x = score, y = p, color = type)) +
  geom_point(data = pop_df) +
  geom_line(data = pop_df) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0, color = "black") +
  stat_summary(fun = 'mean', geom = 'point', size = 2) +
  stat_summary(fun = 'mean', geom = 'line') +
  facet_grid(~trait) +
  theme_classic() +
  ylim(0, 1)
p1 + labs(x = 'Trait score', y = 'Probability', color = NULL)
dev.off()

perceived_pop_dat_summary <- perceived_pop_dat %>%
  group_by(trait, score) %>%
  summarise(p = mean(p, na.rm = TRUE), .groups = "drop") %>% 
  arrange(trait, score)
cor.test(perceived_pop_dat_summary$p, pop_df$p)
# r = 0.581, p = .0008

########## perceived-population counterfactual models ##########
dir.create("cache", showWarnings = FALSE)
dir.create("perceived_population_output", showWarnings = FALSE)

dat <- read.csv("./../data/exp2.csv", header = T, na.strings = c("", "NA", "null")) %>% 
  filter(total_failed_attention_checks == 0) %>% 
  select(-total_failed_attention_checks) 

all_traits <- c("openness", "conscientiousness", "extraversion", "agreeableness", "neuroticism")
all_scores <- c(0, 20, 40, 60, 80, 100)
all_models <- c("actual_contribution", "focal_mutation_counterfactual", "nonfocal_mutation_counterfactual",
                     "focal_population_counterfactual", "nonfocal_population_counterfactual", "ensemble",
                     "focal_perceived_population_counterfactual", "nonfocal_perceived_population_counterfactual")

weights_df <- read.csv("./../code/output/exp1_trait_weights.csv", header = T)
dat_weight <- dat %>% 
  select(-starts_with("responsibility_")) %>% 
  left_join(weights_df, by = "scenario")

# perceived population counterfactual distribution
perceived_pop_df <- perceived_pop_dat_summary %>% 
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

cf_perceived_pop <- make_cf_from_population(perceived_pop_df)

# new raw model predictions
focal_perceived_population_predictions <- responsibility_focal(dat_weight, cf_perceived_pop) %>%
  select(-starts_with("weight")) %>%
  mutate(model = "focal_perceived_population_counterfactual")

nonfocal_perceived_population_predictions <- responsibility_nonfocal(dat_weight, cf_perceived_pop) %>%
  select(-starts_with("weight")) %>%
  mutate(model = "nonfocal_perceived_population_counterfactual")

# new raw model predictions
model_predictions_main <- read.csv("./../code/output/exp2_raw_model_predictions.csv", header = T, na.strings = c("", "NA", "null"))

model_predictions_supp <- rbind(model_predictions_main, 
                                focal_perceived_population_predictions, 
                                nonfocal_perceived_population_predictions)

if (!file.exists("perceived_population_output/exp2_raw_model_predictions.csv")) {
  write.csv(model_predictions_supp,"perceived_population_output/exp2_raw_model_predictions.csv", row.names = FALSE)
}

# prepare long data for regression
predictions <- model_predictions_supp %>%
  mutate(model = factor(model, levels = all_models, ordered = T)) %>%
  filter(!is.na(model)) %>%
  pivot_longer(cols = starts_with("responsibility_"), names_to = "focal_trait", values_to = "predicted_responsibility") %>%
  mutate(focal_trait = sub("^responsibility_", "", focal_trait))

predictions_wide <- predictions %>%
  pivot_wider(id_cols = c("subject", "scenario", "person", "focal_trait"), names_from = "model", values_from = "predicted_responsibility")

dat_long <- dat %>% 
  pivot_longer(cols = starts_with("responsibility_"), names_to = "focal_trait", values_to = "responsibility") %>% 
  mutate(focal_trait = sub("^responsibility_", "", focal_trait)) %>% 
  left_join(predictions_wide, by = c("subject", "scenario", "person", "focal_trait")) %>% 
  mutate(subject = as.factor(subject)) %>%
  mutate(ensemble = (actual_contribution + (focal_mutation_counterfactual + nonfocal_mutation_counterfactual)/2) / 2)

# load existing fits for first 6 models
fits_single <- setNames(vector("list", length(all_models)), all_models)

fits_single[["actual_contribution"]] <- load_fit_with_loo(
  "./../code/cache/resp_single_actual_contribution.rds",
  "./../code/cache/resp_single_actual_contribution_loo"
)
fits_single[["focal_mutation_counterfactual"]] <- load_fit_with_loo(
  "./../code/cache/resp_single_focal_mutation_counterfactual.rds",
  "./../code/cache/resp_single_focal_mutation_counterfactual_loo"
)
fits_single[["nonfocal_mutation_counterfactual"]] <- load_fit_with_loo(
  "./../code/cache/resp_single_nonfocal_mutation_counterfactual.rds",
  "./../code/cache/resp_single_nonfocal_mutation_counterfactual_loo"
)
fits_single[["focal_population_counterfactual"]] <- load_fit_with_loo(
  "./../code/cache/resp_single_focal_population_counterfactual.rds",
  "./../code/cache/resp_single_focal_population_counterfactual_loo"
)
fits_single[["nonfocal_population_counterfactual"]] <- load_fit_with_loo(
  "./../code/cache/resp_single_nonfocal_population_counterfactual.rds",
  "./../code/cache/resp_single_nonfocal_population_counterfactual_loo"
)
fits_single[["ensemble"]] <- load_fit_with_loo(
  "./../code/cache/resp_ensemble.rds",
  "./../code/cache/resp_ensemble_loo"
)

# fit 2 new models and compare
for (m in c("focal_perceived_population_counterfactual", "nonfocal_perceived_population_counterfactual")) {
  fml <- bf(glue("responsibility ~ 1 + {m}"))
  fits_single[[m]] <- fit_with_loo(
    fml, dat_long,
    file_stub = glue("resp_single_{m}")
  )
}

loo_all <- loo::loo_compare(c(
  lapply(fits_single, \(f) f$criteria$loo)
))
print(loo_all)

if (!file.exists("perceived_population_output/exp2_loo_comparison.csv")) {
  write.csv(rank_loo(loo_all),"perceived_population_output/exp2_loo_comparison.csv", row.names = FALSE)
}

# ELPD plot
p2a <- loo_plot(loo_all, model_levels = all_models)

# push model predictions through regression models
for (m in all_models) {
  coefs <- fixef(fits_single[[m]])
  intercept <- coefs["Intercept", "Estimate"]
  slope <- coefs[m, "Estimate"]
  
  dat_long[[paste0(m, "_fitted")]] <- intercept + slope * dat_long[[m]]
}

if (!file.exists("perceived_population_output/exp2_fitted_models.csv")) {
  write.csv(dat_long,"perceived_population_output/exp2_fitted_models.csv", row.names = FALSE)
}

dat_long <- read.csv("perceived_population_output/exp2_fitted_models.csv", header = T, na.strings = c("", "NA", "null"))

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
  filter(model %in% c("focal_perceived_population_counterfactual", "nonfocal_perceived_population_counterfactual")) %>%
  mutate(model = as_factor(model))

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

# scatter plot
p2b <- scatter_plot(dat_model)

pdf("figS2.pdf", onefile = TRUE, width = 16, height = 5)
(p2a | p2b) + plot_annotation(tag_levels = "A")
dev.off()
