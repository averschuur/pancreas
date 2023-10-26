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
         mixed = as.factor(ifelse(tumorType == "Mixed", 1, 0)),
         annotation = tumorType) %>% 
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


# heatmap of sample-wise correlations
sample_cor <- cor(betas_model)
betas_cor <- cor(t(betas_model))

rownames(sample_cor) <- ifelse(anno$tumorType == "Mixed", "Mixed", ".")
#pdf(file = "plots/mixed_tumors_sample_correlation_heatmap.pdf")
pheatmap::pheatmap(sample_cor, show_colnames = FALSE)
#dev.off()

# run UMAP umap
set.seed(45098)
umap_settings <- umap.defaults
umap_settings$n_neighbors = 15
umap_settings$min_dist = 0.2

#umap <- umap(d = sample_cor)
umap <- umap(d = t(betas_model), config = umap_settings, ret_model = TRUE)



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

# plot UMAP
anno <- anno %>% 
  mutate(umap_x = umap$layout[, 1], 
         umap_y = umap$layout[, 2])

anno %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType)) +
  geom_point(size = 4)


# outlier detection

anno_Mixed <- anno %>% 
  mutate(winning_class = c(rf_class), 
         winning_score = c(rf_winning_score), 
         tumorType = anno$tumorType,
         class_int = rep(1, length(rf_class)),
         class_char = ifelse(class_int == 0, "outlier", "pancreas"))

od_pred1 <- anno_Mixed %>%
  predict(object = od_model, newdata = ., type = "response")

anno_Mixed <- anno_Mixed %>% 
  mutate(od_prob1 = od_pred1, 
         od_class = ifelse(od_prob1 > 0.5, "pancreas", "outlier"))

#saveRDS(object = anno_Mixed, file = "./output/anno_incl_mixed tumors.rds")


### add TCGA data ----------------------------------------------------------
# processing was performed on server
anno_tcga <- readRDS("./output/sample_anno_tcga.rds")
betas_tcga <- readRDS(file = "./input/betas_tcga_modelprobes.rds")
betas_tcga <- betas_tcga[, anno_tcga$basename]

# rename TCGA columnns
anno_tcga <- anno_tcga %>% 
  select(basename, tissue) 

colnames(anno_tcga) <- c("arrayId","tumorType")

# rename TCGA ACC (adrenal carcinoma) to ADC (prevent clash with acinar carcinoma)
anno_tcga$tumorType[anno_tcga$tumorType == "ACC"] <- "ADC"

anno_tcga <- anno_tcga %>% 
  mutate(mixed = as.factor(ifelse(tumorType == "Mixed", 1, 0)),
         subtype = tumorType,
         tumorType = "TCGA",
         annotation = subtype)

# RF predictions
rf_scores_tcga <- predict(object = rf_model, t(betas_tcga), type = "prob")
rf_scores_tcga <- rename_with(rf_scores_tcga, ~ paste0("scores_", .x))

# winning score
rf_scores_tcga_max <- apply(rf_scores_tcga, 1, max)

# class
rf_class_tcga <- predict(object = rf_model, t(betas_tcga), type = "raw")

# rf score entropy
rf_entropy_tcga <- apply(rf_scores_tcga, 1, entropy)
beta_var_tcga <- apply(betas_tcga, 2, var)
beta_entropy_tcga <- apply(betas_tcga, 2, entropy)

# combine data
rf_data <- bind_rows(as_tibble(rf_pred_scores), 
                     as_tibble(rf_scores_tcga))

# add annotaion
rf_data <- rf_data %>% 
  mutate(winning_class = c(rf_pred_class, rf_class_tcga), 
         winning_score = c(rf_pred_scores_max, rf_scores_tcga_max),
         rf_entropy = c(rf_entropy, rf_entropy_tcga),
         beta_entropy = c(beta_entropy, beta_entropy_tcga), 
         beta_var = c(beta_var, beta_var_tcga),
         class_int = c(rep(1, length(rf_pred_class)), rep(0, length(rf_class_tcga))),
         class_char = ifelse(class_int == 0, "outlier", "pancreas"))

# combine data
anno_all <- bind_rows(as_tibble(anno), 
                      as_tibble(anno_tcga))

anno_all <- bind_cols(anno_all, rf_data)


# RF score entropy
anno_all %>% 
  ggplot(aes(tumorType, rf_entropy, fill = tumorType)) +
  geom_boxplot() +
  theme_classic(base_size = 18) +
  labs(x = NULL, y = "RF score entropy")



# combine betas
betas_all <- bind_cols(betas_model, betas_tcga)

# plot UMAP
umap_all <- umap(d = t(betas_all), config = umap_settings, ret_model = TRUE)


anno_all <- anno_all %>% 
  mutate(umap_x = umap_all$layout[, 1], 
         umap_y = umap_all$layout[, 2])

anno_all %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType)) +
  geom_point(size = 4)


