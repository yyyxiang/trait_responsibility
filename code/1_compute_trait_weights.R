library(tidyverse)
library(patchwork)

dat <- read.csv("./../data/exp1.csv", header = T, stringsAsFactors = T) %>% 
  filter(total_failed_attention_checks == 0) %>%
  mutate(subject = as_factor(subject))

epsilon <- 0.0001
dat <- dat %>%
  mutate(prob = action_probability / 100,
         prob_adj = pmin(pmax(prob, epsilon), 1 - epsilon), # clipping to handle P=0 or P=1
         logit_prob = qlogis(prob_adj)) # logit: log P/(1-P)

# empty data frame to store results
weights_df <- data.frame(
  scenario = integer(),
  trait = character(),
  weight = numeric(),
  stringsAsFactors = T
)

for (sc in unique(dat$scenario)) {
  dat_sc <- subset(dat, scenario == sc)
  
  m <- lm(logit_prob ~ openness + conscientiousness + extraversion + agreeableness + neuroticism, data = dat_sc)
  coefs <- coef(m)
  
  temp <- data.frame(
    scenario = sc,
    trait = names(coefs),
    weight = as.numeric(coefs),
    stringsAsFactors = T
  )
  
  weights_df <- rbind(weights_df, temp)
}

weights_df <- weights_df %>% 
  mutate(trait = factor(trait, ordered = T,
                        levels = c("(Intercept)", "openness", "conscientiousness", "extraversion", "agreeableness", "neuroticism"))) %>% 
  arrange(scenario, trait) %>% 
  pivot_wider(names_from = "trait", values_from = "weight", names_prefix = "weight_") %>% 
  dplyr::rename("intercept" = `weight_(Intercept)`)

if (!file.exists("output/exp1_trait_weights.csv")) {
  write.csv(weights_df,"output/exp1_trait_weights.csv", row.names = FALSE)
}

# plot the trait weights
long <- weights_df %>%
  pivot_longer(
    cols = -scenario,
    names_to = "trait",
    values_to = "weight"
  ) %>% 
  mutate(
    trait = dplyr::recode(
      trait,
      "intercept" = "Intercept",
      "weight_openness"  = "Openness",
      "weight_conscientiousness" = "Conscientiousness",
      "weight_extraversion" = "Extraversion",
      "weight_agreeableness" = "Agreeableness",
      "weight_neuroticism" = "Neuroticism"
    )
  ) %>% 
  mutate(trait = factor(trait, ordered = T,
                        levels = c("Intercept", "Openness", "Conscientiousness", 
                                   "Extraversion", "Agreeableness", "Neuroticism")),
         ymin = ifelse(trait == "Intercept", -5,  -0.01),
         ymax = ifelse(trait == "Intercept", -2.5, 0.1)
  )
  
trait_colors <- c(
  "Intercept"   = "#323031",
  "Openness"    = "#ED8146",
  "Conscientiousness"    = "#FFC857",
  "Extraversion"    = "#DB3A34",
  "Agreeableness"    = "#A73F40",
  "Neuroticism"    = "#177E89"
)

p1a <- ggplot(long, aes(x = scenario, y = weight, color = trait)) +
  geom_blank(aes(y = ymin)) +   # enforce facet-specific lower limit
  geom_blank(aes(y = ymax)) +   # enforce facet-specific upper limit
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  facet_wrap(~trait, ncol = 3, scales = "free_y") +
  scale_color_manual(values = trait_colors) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none") +
  labs(x = "Scenario", y = "Weight")

# check if the weights predict data
df <- dat %>% 
  left_join(weights_df, by = "scenario")

all_traits <- c("openness","conscientiousness","extraversion","agreeableness","neuroticism")
wcols <- paste0("weight_", all_traits)

# elementwise multiply each trait by its weight, then sum across the row
df$logit <- rowSums(df[all_traits] * df[wcols])
df$logit <- df$logit + df$intercept

# sigmoid
df$predicted_probability <- plogis(df$logit)

p1b <- df %>%
  mutate(bin = cut(predicted_probability, breaks = 10)) %>%
  group_by(bin) %>%
  summarise(predicted = mean(predicted_probability),
            actual = mean(action_probability)) %>%
  ggplot(aes(predicted, actual)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic(base_size = 13) +
  xlim(0, 1.01) +
  labs(x = "Predicted Values", y = "Actual Values")

# Figure 1: estimates and validation of trait weights
pdf("./../figures/fig1_trait_weight.pdf", onefile = T, width = 12, height = 5)
(p1a | p1b) + plot_annotation(tag_levels = "A")
dev.off()
