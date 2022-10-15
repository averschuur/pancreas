# Christoph Geisenberger
# github: @cgeisenberger
# last edited 02/10/2022



### Preparation ----------------------------------------------------------------

library(tidyverse)
library(ggplot2)

library(keras)

library(Rtsne)
library(umap)

source("./code-christoph/0_functions.R")
source("./code-christoph/0_branded_colors.R")



# load data for TCGA cases -----------------------------------------------------

# processing was performed on server
anno_tcga <- read_csv("./input/annotation/cleaned/annotation_tcga_subsampled.csv")
betas_tcga <- readRDS(file = "./input/betas_tcga_subsampled.rds")

# filter out some useless columns from TCGA annotation 
anno_tcga <- anno_tcga %>% 
  select(basename, barcode, tissue)

anno_tcga <- anno_tcga %>% 
  mutate(squamous = ifelse(tissue %in% c("BLCA", "CESC", "ESCA", "HNSC", "LUSC"), "squamous", "other"))

# filter Pancreatic cancer samples
anno_tcga <- anno_tcga %>% 
  filter(tissue != "PAAD")
betas_tcga <- betas_tcga[, anno_tcga$basename]

# create UMAP for TCGA data
variance <- apply(betas_tcga, 1, var)
topvar <- order(variance, decreasing = TRUE)[1:5000]

# set custom settings
umap_settings = umap.defaults
umap_settings$n_neighbors = 15
umap_settings$min_dist = 0.3

umap_input <- betas_tcga[topvar, ]
umap <- umap(d = t(umap_input), config = umap_settings)

anno_tcga <- anno_tcga %>% 
  mutate(umap_x = umap$layout[, 1], 
         umap_y = umap$layout[, 2])

anno_tcga %>% 
  ggplot(aes(umap_x, umap_y, col = tissue)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_classic(base_size = 20) +
  labs(x = "Dim 1", y = "Dim 2")



# run data through neural network ----------------------------------------------


# load model
nn_model <- load_model_hdf5(filepath = "./input/model_nn.hdf5")
model_probes <- readRDS("./input/model_probes.rds")

# predict
input_model <- betas_tcga[model_probes, ]
tcga_scores <- predict(object = nn_model, t(input_model))

# add labels to score matrixs
tcga_labels <- c("ACC", "NORM", "PanNET", "PDAC", "SPN")
colnames(tcga_scores) <- tcga_labels

# calculate max score and label for all cases
tcga_max_score <- apply(tcga_scores, 1, max)
tcga_classes <- apply(tcga_scores, 1, function(x) tcga_labels[which.max(x)])

# add data to annotation
anno_tcga <- anno_tcga %>% 
  mutate(nn_class = tcga_classes, 
         nn_max = tcga_max_score,
         nn_score_acc = tcga_scores[, 1], 
         nn_score_norm = tcga_scores[, 2],
         nn_score_net = tcga_scores[, 3], 
         nn_score_pdac = tcga_scores[, 4], 
         nn_score_spn = tcga_scores[, 5])


anno_tcga <- anno_tcga %>% 
  mutate(mean_beta = apply(betas_tcga, 2, mean))


anno_tcga %>% 
  ggplot(aes(mean_beta, tissue, fill = tissue)) +
  geom_violin() +
  theme_classic(base_size = 20) +
  labs(x = "Avg. Beta", y = NULL) +
  theme(legend.position = "none")


anno_tcga %>% 
  ggplot(aes(nn_max)) +
  geom_histogram(fill = branded_colors2[3], alpha = 0.5, col =  branded_colors2[3]) +
  theme_classic(base_size = 24) +
  labs(x = "Highest Score", y = "Count")


anno_tcga %>% 
  ggplot(aes(tissue, nn_max)) +
  geom_jitter(width = 0.2, col = branded_colors2[3]) +
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = NULL, y = "Highest Score")

anno_tcga %>% 
  group_by(nn_class) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(n, nn_class, fill = nn_class)) +
  geom_col() +
  scale_fill_manual(values = branded_colors) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none") +
  labs(y = NULL, x = "Predicted class (n)")

anno_tcga %>% 
  group_by(nn_class) %>% 
  summarise(n = n()) %>% 
  ungroup %>% 
  mutate(n = n / sum(n))


test <- anno_tcga %>% 
  select(tissue, nn_class) %>% 
  group_by(tissue, nn_class) %>% 
  summarise(n = n()/25) %>% 
  ungroup %>% 
  pivot_wider(id_cols = tissue, names_from = nn_class, values_from = n)

test2 <- as.matrix(test[, 2:5])
test2 <- cbind(test2, rep(0, 16))
colnames(test2)[5] <- "SPN"
rownames(test2) <- test$tissue
test2[is.na(test2)] <- 0

superheat::superheat(test2, 
                     order.rows = 1:16,
                     bottom.label.text.angle = 90,
                     bottom.label.text.size = 6,
                     bottom.label.col = "white",
                     left.label.col = "white",
                     left.label.text.size = 6,
                     heat.pal = c("white",  branded_colors2[3]), 
                     X.text = round(test2, 2))    
  