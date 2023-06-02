# Christoph Geisenberger
# github: @cgeisenberger
# last edited 04/01/2023 by AV Verschuur



# Load packages and source scripts ---------------------------------------------

library(tidyverse)
library(ggplot2)

library(keras)
library(randomForest)

library(Rtsne)
library(umap)

library(caret)
library(pROC)

source("./scripts/0_helpers.R")



# load data for pancreas data set ----------------------------------------------

anno <- readRDS("./output/sample_annotation_umap_purity.rds")
betas <- readRDS(file = "./input/betas_pancreas_everything.rds")
betas <- betas[, anno$arrayId]

# pick model probes
topvar_probes <- readRDS(file = "./output/pancreas_top_variable_probes_training_set.rds")
topvar_probes <- topvar_probes[1:5000]
betas <- betas[topvar_probes, ]



# load data for TCGA cases -----------------------------------------------------

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

# tcga data set stats
anno_tcga %>% 
  group_by(tumorType) %>% 
  summarise(n = n())



# load models ------------------------------------------------------------------

nn_model <- load_model_hdf5(filepath = "./output/nn_model.hdf5")
rf_model <- readRDS(file = "./output/rf_model_default.rds")



# run data through neural network ----------------------------------------------

# determine NN scores
nn_scores_tcga <- predict(object = nn_model, t(betas_tcga))
nn_scores_pancreas <- predict(object = nn_model, t(betas))

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
nn_classes_pancreas <- apply(nn_scores_pancreas, 1, which.max)
nn_classes_pancreas <- nn_labels[nn_classes_pancreas]

# collect data
nn_scores_all <- tibble(class_int = c(rep(0, nrow(nn_scores_tcga)), 
                                      rep(1, nrow(nn_scores_pancreas))), 
                        score = c(nn_scores_tcga_max, nn_scores_pancreas_max))

nn_scores_max <- nn_scores_all %>% 
  mutate(class_char = ifelse(class_int == 0, "tcga", "pancreas"),
         pred_class = c(as.character(nn_classes_tcga), as.character(nn_classes_pancreas))) %>% 
  select(class_char, class_int, score, pred_class)


# plot score distribution 
nn_scores_max %>% 
  ggplot(aes(class_char, score, fill = class_char)) +
  geom_boxplot(alpha = 0.5) +
  theme_bw(base_size = 20) +
  scale_fill_manual(values = branded_colors2) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  labs(x = NULL, y = "Neural Network Score")

# plot ROC curve
nn_scores_roc <- pROC::roc(nn_scores_max$class_int, nn_scores_max$score)
plot.roc(nn_scores_roc, grid = TRUE, )
nn_scores_roc



### run data through random forest ---------------------------------------------

# RF predictions
rf_scores_tcga <- predict(object = rf_model, t(betas_tcga), type = "prob")
rf_scores_pancreas <- predict(object = rf_model, t(betas), type = "prob")

# winning score
rf_scores_tcga_max <- apply(rf_scores_tcga, 1, max)
rf_scores_pancreas_max <- apply(rf_scores_pancreas, 1, max)

# class
rf_class_tcga <- predict(object = rf_model, t(betas_tcga), type = "raw")
rf_class_pancreas <- predict(object = rf_model, t(betas), type = "raw")

# combine data
rf_data <- bind_rows(as_tibble(rf_scores_pancreas, rownames = "arrayId"), 
                     as_tibble(rf_scores_tcga, rownames = "arrayId"))
rf_data <- rf_data %>% 
  mutate(winning_class = c(rf_class_pancreas, rf_class_tcga), 
         winning_score = c(rf_scores_pancreas_max, rf_scores_tcga_max), 
         tumorType = c(anno$tumorType, anno_tcga$tumorType),
         class_int = c(rep(1, length(rf_class_pancreas)), rep(0, length(rf_class_tcga))),
         class_char = ifelse(class_int == 0, "outlier", "pancreas"))


# plot distribution for TCGA classification results
rf_data %>% 
  filter(class_char == "outlier") %>% 
  group_by(winning_class) %>% 
  summarise(n = n()) %>% 
  mutate(prop = n/sum(n)*100) %>% 
  ggplot(aes(winning_class, prop, fill = winning_class)) +
  geom_col() +
  theme_bw(base_size = 24) +
  labs(x = "Predicted class", y = "Proportion of TCGA samples") +
  theme(legend.position = "none") +
  scale_fill_manual(values = branded_colors3)

# plot score distribution 
rf_data %>% 
  ggplot(aes(class_char, winning_score, fill = class_char)) +
  geom_boxplot(alpha = 0.5) +
  theme_bw(base_size = 20) +
  scale_fill_manual(values = branded_colors2) +
  labs(x = NULL, y = "Random Forest Score")

# plot score distribution with resolution per class
rf_data %>% 
  filter(winning_class %in% c("ACC", "PanNET", "PDAC")) %>% 
  ggplot(aes(class_char, winning_score, fill = winning_class)) +
  scale_fill_manual(values = branded_colors3) +
  geom_boxplot(alpha = 0.5) +
  theme_bw(base_size = 20) +
  facet_wrap(facets = vars(winning_class)) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Random Forest Score")

# calculate AUC, plot ROC
rf_data_roc <- with(rf_data, roc(class_int, winning_score))
rf_data_roc
plot.roc(rf_data_roc)



# train outlier detection model ------------------------------------------------

# split data into train and test cohort
set.seed(23432)
rf_data <- rf_data %>% 
  mutate(dataset = sample(c("train", "test"), size = nrow(rf_data), replace = TRUE))

od_model <- rf_data %>% 
  filter(dataset == "train") %>%
  glm(class_int ~ winning_score + winning_class, data = ., family = 'binomial')
summary(od_model)

od_pred <- rf_data %>%
  predict(object = od_model, newdata = ., type = "response")

rf_data <- rf_data %>% 
  mutate(od_prob = od_pred, 
         od_class = ifelse(od_prob > 0.5, "pancreas", "outlier"))

# plot outlier probability
rf_data %>% 
  ggplot(aes(class_char, od_prob, fill = class_char)) +
  geom_boxplot(alpha = 0.5) +
  geom_hline(yintercept = 0.5, lty = 2, col = "steelblue") +
  theme_bw(base_size = 24) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Outlier Probability") +
  scale_fill_manual(values = branded_colors2)

# plot ROC curve
od_prob_roc <- rf_data %>%
  pROC::roc(class_int, od_prob)
plot.roc(od_prob_roc)
od_prob_roc

rf_data %>% 
  filter(dataset == "test") %>%
  with(., table(class_char, od_class))



# use different cutoffs for outlier probability

cutoffs <- seq(from = 0, to = 0.95, length.out = 20)

real_class <- rf_data %>% 
  filter(dataset == "test") %>% 
  pull(class_char)
real_class <- ifelse(real_class == "pancreas", "pancreas", "outlier")

pred_class <- rf_data %>% 
  filter(dataset == "test") %>% 
  pull(od_class) %>% 
  as.character()

score <- rf_data %>% 
  filter(dataset == "test") %>% 
  pull(od_prob)

pass <- NULL
dropout <- NULL

for (i in 1:length(cutoffs)){
  pass_c <- score[real_class == "outlier"]
  pass_c <- sum(pass_c >= cutoffs[i])/length(pass_c)
  pass <- c(pass, pass_c)
  dropout_c <- score[real_class == "pancreas"]
  dropout_c <- sum(dropout_c < cutoffs[i])/length(dropout_c)
  dropout <- c(dropout, dropout_c)
}

cutoff_data <- tibble(cutoff = cutoffs, 
                         pass = pass * 100, 
                         dropout = dropout * 100)

cutoff_data %>% 
  pivot_longer(cols = -cutoff) %>% 
  filter(cutoff > 0) %>% 
  ggplot(aes(cutoff, value, col = name)) + 
    geom_line(lwd = 2) +
  geom_vline(xintercept = 0.5, col = "steelblue", lty = 2) +
  facet_wrap(facets = vars(name), scales = "free", ncol = 1) +
  theme_bw(base_size = 24) +
  theme(legend.position = "none", legend.title = element_blank()) +
  labs(x = "Cutoff (outlier prob)", y = "Samples (%)")


### save results to disk ------------------------
saveRDS(object = od_model, file = "./output/outlier_det_model.rds")
saveRDS(object = rf_data, file = "./output/outlier_det_results.rds")
