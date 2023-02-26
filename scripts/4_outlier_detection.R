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
  select(basename, tissue) %>% 
  rename("arrayId" = basename, 
         "tumorType" = tissue)

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

heatmap_anno <- as.data.frame(anno_tcga[, 2])
rownames(heatmap_anno) <- anno_tcga$arrayId
pheatmap::pheatmap(rf_scores_tcga, annotation_row = heatmap_anno, show_rownames = FALSE)

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
rf_scores_all <- tibble(arrayId = c(rownames(rf_scores_tcga), 
                                    rownames(rf_scores_pancreas)), 
                        class_int = c(rep(0, nrow(rf_scores_tcga)), 
                                      rep(1, nrow(rf_scores_pancreas))), 
                        score = c(rf_scores_tcga_max, rf_scores_pancreas_max))

rf_scores_all <- rf_scores_all %>%
  mutate(class_char = ifelse(class_int == 0, "tcga", "pancreas"), 
         pred_class = c(as.character(rf_class_tcga), as.character(rf_class_pancreas)))


# plot distribution for TCGA classification results
rf_scores_all %>% 
  filter(class_char == "tcga") %>% 
  group_by(pred_class) %>% 
  summarise(n = n()/nrow(rf_scores_all)) %>% 
  ggplot(aes(pred_class, n, fill = pred_class)) +
  geom_col() +
  theme_bw(base_size = 24) +
  labs(x = "Predicted class", y = "Proportion of TCGA samples") +
  theme(legend.position = "none") +
  scale_fill_manual(values = branded_colors3)

# plot score distribution 
rf_scores_all %>% 
  ggplot(aes(class_char, score, fill = class_char)) +
  geom_boxplot(alpha = 0.5) +
  theme_bw(base_size = 20) +
  scale_fill_manual(values = branded_colors2) +
  labs(x = NULL, y = "Random Forest Score")

# plot score distribution with resolution per class
rf_scores_all %>% 
  filter(pred_class %in% c("ACC", "PanNET", "PDAC")) %>% 
  ggplot(aes(class_char, score, fill = class_char)) +
  scale_fill_manual(values = branded_colors2) +
  geom_boxplot(alpha = 0.5) +
  theme_bw(base_size = 20) +
  facet_wrap(facets = vars(pred_class)) +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL)

# calculate AUC, plot ROC
rf_scores_roc <- with(rf_scores_all, roc(class_int, score))
rf_scores_roc
plot.roc(rf_scores_roc)



# train outlier detection model ------------------------------------------------

# split data
rf_scores_all <- rf_scores_all %>% 
  mutate(dataset = sample(c("train", "val"), size = nrow(rf_scores_all), replace = TRUE))

od_model <- rf_scores_all %>% 
  filter(dataset == "train") %>% 
  glm(class_int ~ score + pred_class , data = ., family = 'binomial')

od_pred <- rf_scores_all %>%
  predict(object = od_model, newdata = .)

rf_scores_all <- rf_scores_all %>% 
  mutate(od_logit = od_pred)

# convert logit to class
rf_scores_all <- rf_scores_all %>% 
  mutate(od_class = ifelse(od_logit > 0, "pancreas", "outlier"))

# plot logit scores
rf_scores_all %>% 
  ggplot(aes(class_char, od_logit, fill = class_char)) +
  geom_boxplot(alpha = 0.5) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = branded_colors2)

# plot ROC curve
logit_scores_roc <- rf_scores_all %>%
  pROC::roc(class_int, od_logit)
plot.roc(logit_scores_roc)
logit_scores_roc

# confusion matrix
rf_scores_all %>%  
  filter(dataset == "val") %>% 
  with(., table(class_char, od_class))

# add outlier information back to pancreas annotation 
anno <- left_join(anno, rf_scores_all[, c("arrayId", "od_logit", "od_class")])

with(rf_scores_all, table(class_char, od_class))

### plot UMAP with outlier detection results -----------------------------------

table(anno$od_class)

anno %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType, shape = od_class)) +
  geom_point(size = 5) +
  theme_classic(base_size = 10)
