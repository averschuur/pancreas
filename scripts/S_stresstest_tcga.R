### Load required packages and sources ----------------------------------------------------------------

library(tidyverse)
library(ggplot2)

library(keras)

library(Rtsne)
library(umap)

source("./scripts/functions.R")
source("./scripts/branded_colors.R")


### load data for TCGA cases -----------------------------------------------------

# processing was performed on server
anno_tcga <- read_csv("./annotation/cleaned/annotation_tcga_subsampled.csv")
betas_tcga <- readRDS(file = "./data/betas_tcga_subsampled.rds")

# filter out some useless columns from TCGA annotation 
anno_tcga <- anno_tcga %>% 
  select(basename, barcode, tissue)

anno_tcga <- anno_tcga %>% 
  mutate(squamous = ifelse(tissue %in% c("BLCA", "CESC", "ESCA", "HNSC", "LUSC"), "squamous", "other"))

# filter Pancreatic cancer samples
anno_tcga <- anno_tcga %>% 
  filter(tissue != "PAAD")
betas_tcga <- betas_tcga[, anno_tcga$basename]
