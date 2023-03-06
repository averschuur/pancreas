# Christoph Geisenberger
# github: @cgeisenberger
# last edited 01/03/2023 by CG



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


# load annotation
anno <- readRDS(file = "./output/sample_annotation_umap_purity.rds")

# add outlier detection results
od <- readRDS(file = "./output/outlier_det_results.rds")
od <- od %>% 
  select(-tumorType)
anno <- left_join(anno, od, by = "arrayId")
rm(od)

# load data for cell type deconvolution
decomp <- read_csv("./test/meth_atlas/pancreas_data_deconv_output.csv")
decomp_rownames <- decomp$cell_type
decomp <- as.matrix(decomp[, -1])
rownames(decomp) <- decomp_rownames
rm(decomp_rownames)
decomp <- decomp[, anno$arrayId]

# first, look at misclassification rate per tumor type
anno <- anno %>% 
  mutate(class_correct = ifelse(winning_class == tumorType, 1, 0))
anno %>% 
  group_by(tumorType) %>% 
  summarise(no_correct = sum(class_correct == 0), 
            prop_correct = sum(class_correct)/n())


# look at data for specific tumor types
anno_sub <- anno %>% 
  filter(tumorType == "PanNET")

anno_sub_df <- anno_sub %>% 
  select(estimate, class_correct, od_class) %>% 
  mutate(class_correct = as.factor(class_correct)) %>% 
  as.data.frame
rownames(anno_sub_df) <- anno_sub$arrayId

decomp_sub <- decomp[, anno_sub$arrayId]

pheatmap::pheatmap(decomp_sub, annotation_col = anno_sub_df)
    





