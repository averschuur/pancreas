# Christoph Geisenberger
# github: @cgeisenberger
# last edited 04/01/2023 by AV Verschuur


### Load required packages and sources ----------------------------------------------------------------

library(tidyverse)
library(ggplot2)

library(keras)

library(Rtsne)
library(umap)

source("./scripts/0_functions.R")
source("./scripts/0_branded_colors.R")

# load data for pancreas data set ----------------------------------------------

anno <- readRDS("./input/sample_annotation_umap_purity.rds")
betas <- readRDS(file = "./input/pancreas_betas_everything.rds")
betas <- betas[, anno$arrayId]


# load data for TCGA cases -----------------------------------------------------

# processing was performed on server
anno_tcga <- read_csv("./annotation/TCGA/annotation_tcga_subsampled.csv")
betas_tcga <- readRDS(file = "./input/betas_tcga_subsampled.rds")
betas_tcga <- betas_tcga[, anno_tcga$arrayId]


# remove some columns from TCGA annotation 
anno_tcga <- anno_tcga %>% 
  select(arrayId, sampleName, tumorType)

anno_tcga %>% 
  group_by(tumorType) %>% 
  summarise(n = n())

# rename adrenocortical carcinoma 
anno_tcga$tumorType[anno_tcga$tumorType == "ACC"] <- "ADC"



# load models and model probes -------------------------------------------------

top_var_probes <- readRDS(file = "./input/pancreas_top_variable_probes.rds")
nn_model <- load_model_hdf5(filepath = "./output/model_nn.hdf5")
rf_model <- readRDS(file = "./output/model_rf.rds")


# run data through neural network ----------------------------------------------

# determine NN scores
nn_scores_tcga <- predict(object = nn_model, t(betas_tcga[model_probes, ]))
nn_scores_pancreas <- predict(object = nn_model, t(betas[model_probes, ]))

# add labels to score matrixs
nn_labels <- c("ACC", "NORM", "PanNEC", "PanNET", "PB", "PDAC", "SPN")
colnames(nn_scores_tcga) <- nn_labels
colnames(nn_scores_pancreas) <- nn_labels

# calculate max score and label for TCGA cases
nn_scores_tcga_max <- apply(nn_scores_tcga, 1, max)
nn_classes_tcga <- apply(nn_scores_tcga, 1, which.max)
nn_classes_tcga <- nn_labels[nn_classes_tcga]

# calculate scores for Pancreas Samples
nn_scores_pancreas_max <- apply(nn_scores_pancreas, 1, max)

# collect data
nn_scores_all <- tibble(class_int = c(rep(0, nrow(nn_scores_tcga)), 
                                      rep(1, nrow(nn_scores_pancreas))), 
                        score = c(nn_scores_tcga_max, nn_scores_pancreas_max))

nn_scores_max <- nn_scores_all %>% 
  mutate(class_char = ifelse(class_int == 0, "tcga", "pancreas")) %>% 
  select(class_char, class_int, score)


# plot score distribution 
nn_scores_max %>% 
  ggplot(aes(class_char, score, fill = class_char)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  labs(x = NULL)

# plot ROC curve
nn_scores_roc <- roc(nn_scores_max$class_int, nn_scores_max$score)
plot.roc(nn_scores_roc)

# create tibble with NeuralNet results for each TCGA class
nn_label_prop <- anno_tcga %>% 
  select(tumorType, nn_class) %>% 
  group_by(tumorType, nn_class) %>% 
  summarise(n = n()/25) %>% 
  ungroup %>% 
  pivot_wider(id_cols = tumorType, names_from = nn_class, values_from = n)

# replace NAs with 0's, plot as heatmap
nn_label_prop_matrix <- as.matrix(nn_label_prop[, 2:7])
rownames(nn_label_prop_matrix) <- nn_label_prop$tumorType
nn_label_prop_matrix[is.na(nn_label_prop_matrix)] <- 0

pheatmap::pheatmap(nn_label_prop_matrix)



### run data through random forest ---------------------------------------------

# RF predictions
rf_scores_tcga <- predict(object = rf_model, t(betas_tcga[model_probes, ]), type = "prob")
rf_scores_pancreas <- predict(object = rf_model, t(betas[model_probes, ]), type = "prob")

# winning score
rf_scores_tcga_max <- apply(rf_scores_tcga, 1, max)
rf_scores_pancreas_max <- apply(rf_scores_pancreas, 1, max)

# class
rf_class_tcga <- predict(object = rf_model, t(betas_tcga[model_probes, ]), type = "response")
rf_class_pancreas <- predict(object = rf_model, t(betas[model_probes, ]), type = "response")

# combine data
rf_scores_all <- tibble(class_int = c(rep(0, nrow(rf_scores_tcga)), 
                                      rep(1, nrow(rf_scores_pancreas))), 
                        score = c(rf_scores_tcga_max, rf_scores_pancreas_max))

rf_scores_all <- rf_scores_all %>%
  mutate(class_char = ifelse(class_int == 0, "tcga", "pancreas"), 
         pred_class = c(as.character(rf_class_tcga), as.character(rf_class_pancreas)))

# plot score distributions
rf_scores_all %>% 
  ggplot(aes(score, fill = class_char)) +
  geom_density(alpha = 0.4) +
  theme_bw(20) + facet_wrap(facets = vars(pred_class))

# calculate AUC, plot ROC
rf_scores_roc <- with(rf_scores_all, roc(class_int, score))
rf_scores_roc
plot.roc(rf_scores_roc)


# implement outlier detection algorith ------------------

y <- rf_scores_all$class_int %>% as.vector()
x  <- rbind(rf_scores_tcga, rf_scores_pancreas)
split = sample(x = c("train", "val"), size = nrow(x), replace = TRUE, prob = c(0.7, 0.3))
x_train <- x[split == "train", ]
x_val <- x[split == "val", ]
y_train <- y[split == "train"]
y_val <- y[split == "val"]

fit_rf <- randomForest(x = x_train, y = as.factor(y_train))

pred_rf_train <- predict(fit_rf, newdata = x_train, type = "prob")
pred_rf_val <- predict(fit_rf, newdata = x_val, type = "prob")

pred_rf_class_train <- predict(fit_rf, newdata = x_train, type = "class")
pred_rf_class_val <- predict(fit_rf, newdata = x_val, type = "class")

pred_rf <- tibble(arrayId = c(rownames(x_train), rownames(x_val)), 
                  outlier_score = c(pred_rf_train[, 1], pred_rf_val[, 1]), 
                  outlier_class = c(pred_rf_class_train, pred_rf_class_val))
boxplot(pred_rf_val[, 1] ~ as.factor(y_val))

fit_rf_roc <- roc(y_val, pred_rf[, 1])
fit_rf_roc
plot.roc(fit_rf_roc)
with(fit_rf_roc, plot(thresholds, sensitivities))
with(fit_rf_roc, plot(thresholds, specificities))


anno <- left_join(anno, pred_rf)

anno %>% 
  filter(outlier_score < 0.4) %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType)) +
  geom_point() +
  theme_bw(base_size = 20)

# check umcu samples

betas_umc <- readRDS("./input/pancreas_betas_everything.rds")
betas_umc <- betas_umc[, anno_umc$arrayId]

anno_umc <- read_csv("./annotation/annotation_umcu.csv")

umc_scores_rf <- predict(object = rf_model, t(betas_umc[model_probes, ]), type = "prob")
umc_class_rf <- predict(object = rf_model, t(betas_umc[model_probes, ]), type = "class")

anno_umc <- anno_umc %>% 
  mutate(rf_class = as.vector(umc_class_rf), 
         rf_score = apply(umc_scores_rf, 1, max))
umc_outlier_class <- predict(fit_rf, newdata = umc_scores_rf)
umc_outlier_score <- predict(fit_rf, newdata = umc_scores_rf, type = "prob")

anno_umc <- anno_umc %>% 
  mutate(outlier = umc_outlier_class)

anno_umc