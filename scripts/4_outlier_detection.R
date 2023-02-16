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

library(pROC)

source("./scripts/0_helpers.R")


# load data for pancreas data set ----------------------------------------------

anno <- readRDS("./output/sample_annotation_umap_purity.rds")
betas <- readRDS(file = "./input/pancreas_betas_everything.rds")
betas <- betas[, anno$arrayId]


# load data for TCGA cases -----------------------------------------------------

# processing was performed on server
anno_tcga <- read_csv("./input/annotation_tcga_validation.csv")

betas_tcga <- readRDS(file = "./input/betas_tcga_modelprobes.rds")
betas_tcga <- betas_tcga[, anno_tcga$basename]

# rename TCGA columnns
anno_tcga <- anno_tcga %>% 
  select(basename, tissue) %>% 
  rename("arrayId" = basename, 
         "tumorType" = tissue)

# tcga data set stats
anno_tcga %>% 
  group_by(tumorType) %>% 
  summarise(n = n())



# load models and model probes -------------------------------------------------

model_probes <- readRDS(file = "./output/pancreas_top_variable_probes.rds")[1:5000]
nn_model <- load_model_hdf5(filepath = "./output/nn_model.hdf5")
rf_model <- readRDS(file = "./output/rf_model_default.rds")


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
  geom_boxplot() +
  theme_bw(base_size = 20) +
  labs(x = NULL)

# plot score distributions
nn_scores_max %>% 
  ggplot(aes(class_char, score, fill = class_char)) +
  geom_boxplot() +
  theme_bw(20) + facet_wrap(facets = vars(pred_class)) +
  theme(legend.position = "none")

# plot ROC curve
nn_scores_roc <- pROC::roc(nn_scores_max$class_int, nn_scores_max$score)
plot.roc(nn_scores_roc)
nn_scores_roc


### run data through random forest ---------------------------------------------

# RF predictions
rf_scores_tcga <- predict(object = rf_model, t(betas_tcga[model_probes, ]), type = "prob")
rf_scores_pancreas <- predict(object = rf_model, t(betas[model_probes, ]), type = "prob")

# winning score
rf_scores_tcga_max <- apply(rf_scores_tcga, 1, max)
rf_scores_pancreas_max <- apply(rf_scores_pancreas, 1, max)

# class
rf_class_tcga <- predict(object = rf_model, t(betas_tcga[model_probes, ]), type = "raw")
rf_class_pancreas <- predict(object = rf_model, t(betas[model_probes, ]), type = "raw")

# combine data
rf_scores_all <- tibble(class_int = c(rep(0, nrow(rf_scores_tcga)), 
                                      rep(1, nrow(rf_scores_pancreas))), 
                        score = c(rf_scores_tcga_max, rf_scores_pancreas_max))

rf_scores_all <- rf_scores_all %>%
  mutate(class_char = ifelse(class_int == 0, "tcga", "pancreas"), 
         pred_class = c(as.character(rf_class_tcga), as.character(rf_class_pancreas)))

# plot score distribution 
rf_scores_all %>% 
  ggplot(aes(class_char, score, fill = class_char)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  labs(x = NULL)

# plot score distributions
rf_scores_all %>% 
  filter(pred_class != "NORM") %>% 
  ggplot(aes(score, fill = class_char)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = branded_colors2) +
  theme_bw(20) + facet_wrap(facets = vars(pred_class)) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  labs(x = "Random Forest Score", y = "Density")

# calculate AUC, plot ROC
rf_scores_roc <- with(rf_scores_all, roc(class_int, score))
rf_scores_roc
plot.roc(rf_scores_roc)



# train outlier detection model ------------------------------------------------

# split data
rf_scores_all <- rf_scores_all %>% 
  mutate(dataset = sample(c("train", "val"), size = nrow(rf_scores_all), replace = TRUE))


log.model <- rf_scores_all %>% 
  filter(dataset == "train") %>% 
  glm(class_int ~ score + pred_class , data = ., family = 'binomial')


log.pred <- rf_scores_all %>% 
  filter(dataset == "val") %>%
  predict(object = log.model, newdata = .)

log.pred <- rf_scores_all %>% 
  filter(dataset == "val") %>%
  mutate(logit = log.pred)

log.pred %>% 
  ggplot(aes(class_char, logit)) +
  geom_boxplot()

# plot ROC curve
logit_scores_roc <- log.pred %>%
  pROC::roc(class_int, logit)
plot.roc(logit_scores_roc)
logit_scores_roc




# check umcu samples
anno_umc <- read_csv("./annotation/annotation_umcu.csv")
anno_umc <- anno_umc %>% 
  filter(!sampleName == "ACC2")

betas_umc <- readRDS("./input/pancreas_betas_everything.rds")
betas_umc <- betas_umc[, anno_umc$arrayId]


umc_scores_rf <- predict(object = rf_model, t(betas_umc[model_probes, ]), type = "prob")
umc_class_rf <- predict(object = rf_model, t(betas_umc[model_probes, ]), type = "raw")

anno_umc <- anno_umc %>% 
  mutate(rf_class = as.vector(umc_class_rf), 
         rf_score = apply(umc_scores_rf, 1, max))
umc_outlier_score <- predict(log.model, newdata = umc_scores_rf)

anno_umc <- anno_umc %>% 
  mutate(outlier = umc_outlier_class)

anno_umc
