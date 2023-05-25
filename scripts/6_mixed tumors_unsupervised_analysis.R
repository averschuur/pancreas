## Analysis with Mixed tumors


# Load required packages and sources
library(minfi)
library(RFpurify)

library(tidyverse)

library(Rtsne)
library(umap)
library(keras)
library(caret)

source("./scripts/0_helpers.R")



### Prepare data ----------------------------------------------------------------


# load anno, RF model and methylation data
anno <- readRDS("./input/sample_annotation.rds")
rf_model <- readRDS(file = "./output/rf_model_default.rds")
betas <- readRDS(file = "./input/betas_pancreas_everything.rds")

# include mixed tumors
anno <- anno %>% 
  mutate(subtype = tumorType, 
         tumorType = str_replace(tumorType, "Mixed.*", "Mixed"), 
         mixed = as.factor(ifelse(tumorType == "Mixed", 1, 0))) %>% 
  filter(tumorType %in% c("PB", "ACC", "PanNEC", "PanNET", "NORM", "PDAC", "SPN", "Mixed")) %>% 
  filter(location %in% c("primary", "pancreas"))

# filter betas
betas <- betas[, anno$arrayId]
betas_model <- betas[names(rf_model$trainingData)[1:5000], ]

# rf predictions
rf_pred_class <- predict(rf_model, newdata = t(betas))
rf_pred_scores <- predict(rf_model, newdata = t(betas), type = "prob") %>% 
  as_tibble()
rf_pred_scores <- rename_with(rf_pred_scores, ~ paste0("scores_", .x))
rf_pred_scores_max <- apply(rf_pred_scores, 1, max)

# rf score entropy
rf_entropy <- apply(rf_pred_scores, 1, entropy)
beta_var <- apply(betas_model, 2, var)
beta_entropy <- apply(betas_model, 2, entropy)

# purity estimates
estimate <- RFpurify::predict_purity_betas(betas = betas, method = "ESTIMATE")
absolute <- RFpurify::predict_purity_betas(betas = betas, method = "ABSOLUTE")

# add to annotation
anno <- anno %>% 
  mutate(estimate = estimate, 
         absolute = absolute, 
         rf_class = rf_pred_class, 
         rf_winning_score = rf_pred_scores_max, 
         rf_entropy = rf_entropy,
         beta_entropy = beta_entropy, 
         beta_var = beta_var)

anno <- bind_cols(anno, rf_pred_scores)

# check correlations between variables
anno %>% 
  select(tumorType, mixed, absolute, estimate, 
         starts_with("rf"), starts_with("beta")) %>% 
  GGally::ggpairs()

ggsave(filename = "./plots/mixed_tumors_pairwise_plots.pdf")

# heatmap of sample-wise correlations
sample_cor <- cor(betas_model)
betas_cor <- cor(t(betas_model))

rownames(sample_cor) <- ifelse(anno$tumorType == "Mixed", "Mixed", ".")
#pdf(file = "plots/mixed_tumors_sample_correlation_heatmap.pdf")
pheatmap::pheatmap(sample_cor, show_colnames = FALSE)
#dev.off()

# run UMAP umap
set.seed(12341234)
umap <- umap(d = sample_cor)



# look into detail regarding mixed tumors --------------------------------------

# RF score entropy
anno %>% 
  ggplot(aes(mixed, rf_entropy, fill = mixed)) +
  geom_boxplot() +
  theme_classic(base_size = 18) +
  labs(x = NULL, y = "RF score entropy")

anno %>% 
  ggplot(aes(mixed, estimate, fill = mixed)) +
  geom_boxplot() +
  theme_bw(base_size = 18) +
  labs(x = NULL, y = "Tumor Purity")


anno %>% 
  select(tumorType, starts_with("scores")) %>% 
  pivot_longer(cols = starts_with("scores")) %>% 
  mutate(name = str_replace_all(name, "scores_", "")) %>% 
  mutate(name = str_replace_all(name, "Pan", "")) %>% 
  group_by(tumorType, name) %>% 
  summarise(mean_score = mean(value)) %>% 
  ungroup %>% 
  ggplot(aes(name, mean_score, fill = name)) +
  geom_col() +
  theme_bw() +
  labs(y = "Random Forest Score (mean)", fill = "Tumor Type") +
  facet_wrap(facets = vars(tumorType))


anno %>% 
  filter(tumorType %in% c("PB", "ACC", "Mixed")) %>% 
  ggplot(aes(scores_ACC, scores_PB)) +
  geom_point(aes(size = mixed, col = tumorType)) +
  theme_bw(base_size = 18) +
  labs(x = "Score ACC", y = "Score PB")


browseVignettes("EpiSCORE")

