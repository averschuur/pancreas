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
# load sample betas and annotation and filter 

anno <- readRDS("./output/sample_annotation_umap_purity.rds")
betas <- readRDS(file = "./input/betas_pancreas_everything.rds")

betas <- betas[, anno$arrayId]

# obtain model performance 
rf_model <- readRDS(file = "./output/rf_model_default.rds")
#xgb_model <- readRDS(file = "./output/xgb_model_default.rds")
#nn_model <- load_model_hdf5(file = "./output/nn_model.hdf5")

# rf predictions
rf_pred_class <- predict(rf_model, newdata = t(betas))
rf_pred_scores <- predict(rf_model, newdata = t(betas), type = "prob")
rf_pred_scores_max <- apply(rf_pred_scores, 1, max)
rf_pred_class <- apply(rf_pred_scores, 1, function(x) colnames(rf_pred_scores)[which.max(x)]) %>% as.factor


# add performance to annotation
anno <- anno %>% 
  mutate(pred_rf = rf_pred_class, 
         pred_scores_rf = apply(rf_pred_scores, 1, max))


# (1) Calculate the correlations for the RF output matrix for all samples ----------------------------------------
anno %>% 
  ggplot(aes(pred_rf, pred_scores_rf)) +
  geom_boxplot()

cor(rf_pred_scores)
cor(rf_pred_scores$ACC, rf_pred_scores$PB)



# (2) Calculate the correlations for the RF output matrix for each class -----------------------------------------
names <- rownames(rf_pred_scores)
rownames(rf_pred_scores) <- NULL
data <- cbind(names,rf_pred_scores)
colnames(data) <- c("arrayId", "ACC", "NORM", "PanNEC", "PanNET", "PB", "PDAC", "SPN")

# take only ACC samples
ACC <- anno %>%
  filter(tumorType == "ACC") %>%
  select("arrayId")

ACC <- left_join(ACC, data, by = "arrayId")
ACC <- column_to_rownames(ACC, var = "arrayId")
cor(ACC)


# take all PB samples
PB <- anno %>%
  filter(tumorType == "PB") %>%
  select("arrayId")

PB <- left_join(PB, data, by = "arrayId")
PB <- column_to_rownames(PB, var = "arrayId")
cor(PB)

# take all SPN samples
SPN <- anno %>%
  filter(tumorType == "SPN") %>%
  select("arrayId")

SPN <- left_join(SPN, data, by = "arrayId")
SPN <- column_to_rownames(SPN, var = "arrayId")
cor(SPN)

# take all PanNET samples
PNET <- anno %>%
  filter(tumorType == "PanNET") %>%
  select("arrayId")

PNET <- left_join(PNET, data, by = "arrayId")
PNET <- column_to_rownames(PNET, var = "arrayId")
cor(PNET)

# take all PanNEC samples
PNEC <- anno %>%
  filter(tumorType == "PanNEC") %>%
  select("arrayId")

PNEC <- left_join(PNEC, data, by = "arrayId")
PNEC <- column_to_rownames(PNEC, var = "arrayId")
cor(PNEC)

# take all NORM samples
NORM <- anno %>%
  filter(tumorType == "NORM") %>%
  select("arrayId")

NORM <- left_join(NORM, data, by = "arrayId")
NORM <- column_to_rownames(NORM, var = "arrayId")
cor(NORM)


# (3) Calculate the correlations for the RF output matrix for the mixed tumors -----------------------------------
anno <- readRDS("./input/sample_annotation.rds")

# select Mixed tumors 
annoM <- anno %>%
  filter(tumorType %in% c("Mixed", "Mixed_ACC_DA", "Mixed_ACC_NEC", "Mixed_ACC_DA_NEC"))

# load unfiltered beta values
betas <- readRDS("./input/betas_pancreas_everything.rds")
betasM <- betas[, annoM$arrayId]

# obtain model performance 
rf_model <- readRDS(file = "./output/rf_model_default.rds")

# rf predictions
rf_pred_class <- predict(rf_model, newdata = t(betasM))
rf_pred_scores <- predict(rf_model, newdata = t(betasM), type = "prob")
rf_pred_scores_max <- apply(rf_pred_scores, 1, max)
rf_pred_class <- apply(rf_pred_scores, 1, function(x) colnames(rf_pred_scores)[which.max(x)]) %>% as.factor


# add performance to annotation
annoM <- annoM %>% 
  mutate(pred_rf = rf_pred_class, 
         pred_scores_rf = apply(rf_pred_scores, 1, max))

# get correlation
cor(rf_pred_scores)

# (3) calculate entropy ------------------------------------------------------------------------------------------
-sum(rf_pred_scores * log2(rf_pred_scores))

