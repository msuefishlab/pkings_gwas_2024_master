library(lme4)
library(lmerTest)
library(emmeans) # For post-hoc tests

# Step 1: Model only peak regions to test for BP vs TP differences
glm_peak <- lm(avg_pi ~ phenotype_pop1 + population + peak, 
               data = pi_res_group %>% filter(category == "peak"))

# View model summary
summary(glm_peak)

# Perform an ANOVA to check significance of BP vs TP effect
anova(glm_peak)

# Step 2: Post-hoc tests to identify which populations are different
emmeans(glm_peak, pairwise ~ population, adjust = "bonferroni")

# Step 3: Model including peak vs. random for significant populations
significant_pops <- c("APA", "BAM")  # Example; replace with actual significant populations

glm_full <- lm(avg_pi ~ phenotype_pop1 * category + population + peak, 
               data = pi_res_group %>% filter(population %in% significant_pops))

# View model summary
summary(glm_full)

# Perform an ANOVA to check significance of peak vs. random effect within BP populations
anova(glm_full)
