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

# load data for cel type deconvolution
decomp <- read_csv("./test/meth_atlas/pancreas_data_deconv_output.csv")

decomp_mat <- as.matrix(decomp[, -1])
rownames(decomp_mat) <- decomp$cell_type


decomp_umap <- umap(t(decomp_mat))




clusters_n <- 7
km <- kmeans(x = cor(decomp_mat), centers = clusters_n)
hc <- (1 - cor(decomp_mat)) %>% 
  as.dist %>% 
  hclust %>% 
  cutree(k = clusters_n)

decomp_anno <- tibble(arrayId = colnames(decomp_mat), 
                      umap_decomp_x = decomp_umap$layout[, 1], 
                      umap_decomp_y = decomp_umap$layout[, 2], 
                      clusters_km = km$cluster, 
                      clusters_hc = hc)

anno <- readRDS(file = "./output/sample_annotation_umap_purity.rds")

anno <- left_join(anno, decomp_anno)

# tumor type
anno %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType)) +
  geom_point(size = 4)

# kmeans
anno %>% 
  ggplot(aes(umap_x, umap_y, col = as.factor(clusters_km))) +
  geom_point(size = 4)

# hc
anno %>% 
  ggplot(aes(umap_x, umap_y, col = as.factor(clusters_hc))) +
  geom_point(size = 4)


decomp <- decomp %>% 
  pivot_longer(cols = 2:437, names_to = "arrayId", values_to = "prop")
decomp <- decomp %>% 
  select(arrayId, cell_type, prop)
colnames(decomp) <- decomp[, 1]

anno <- left_join(anno, decomp, multiple = "all")


test <- anno %>% 
  #group_by(tumorType, cell_type) %>% 
  #summarise(mean_prop = mean(prop)) %>% 
  filter(tumorType == "PanNET") %>% 
  select(arrayId, cell_type, prop) %>% 
  pivot_wider(names_from = cell_type, values_from = prop) %>% 
  select(-1) %>% 
  as.matrix() %>% 
  pheatmap::pheatmap()
ggplot(aes(arrayId, cell_type, size = prop, col = prop)) +
  geom_point()
theme(legend.position = "none") +
  facet_wrap(facets = vars(tumorType))
ungroup() %>% 
  select(2:26) %>% 
  pheatmap::pheatmap(cluster_rows = FALSE)
ggplot(aes(avg_beta, prop, col = cell_type)) +
  geom_col()
