# libraries
library(RFpurify)
library(minfi)
library(tidyverse)
library(Rtsne)
library(umap)

# load data and annotation 

betas <- readRDS("./00_christoph_data/methylation_data.rds")
anno <- readRDS("./00_christoph_annotation/sample_annotation.rds")


# keep only primaries of PNETs, ACCs and SPNs
anno <- anno %>% 
  filter(tumorType %in% c("PanNET", "ACC", "SPN")) %>% 
  filter(location == "primary")

# filter for samples where annotation is available
betas <- betas[, anno$arrayId]

# get purity estimates
absolute <- RFpurify::predict_purity_betas(betas = betas, method = "ABSOLUTE")
estimate <- RFpurify::predict_purity_betas(betas = betas, method = "ESTIMATE")

# collect data into tibble
anno <- anno %>% 
  mutate(absolute = absolute,
         estimate = estimate,
         avg_beta = apply(betas, 2, mean, na.rm = TRUE))

# make some basic plots
anno %>% 
  ggplot(aes(absolute, estimate)) +
  geom_point(size = 3) +
  geom_smooth() +
  theme_bw(base_size = 18) +
  labs(x = "Purity (ABSOLUTE)", y = "Purity (ESTIMATE)")

anno %>% 
  ggplot(aes(tumorType, estimate, col = tumorType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  theme_bw(base_size = 18)

anno %>% 
  ggplot(aes(tumorType, avg_beta, col = tumorType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  theme_bw(base_size = 18)



### unsupervised analysis: tSNE & UMAP ----------------------------------

# select most variable probes
most_var <- apply(betas, 1, mad)
most_var <- order(most_var, decreasing = TRUE)[1:5000]

# calculate correlation matrix, use 1-cor as dissimlarity measure
cor_matrix <- 1 - cor(betas[most_var, ])

# run tSNE
tsne <- Rtsne(X = cor_matrix)

# run UMAP
umap <- umap(d = t(betas[most_var, ]))

anno <- anno %>% 
  mutate(tsne_x = tsne$Y[, 1],
         tsne_y = tsne$Y[, 2], 
         umap_x = umap$layout[, 1], 
         umap_y = umap$layout[, 2])


# plot t-SNE map
anno %>% 
  ggplot(aes(tsne_x, tsne_y, col = tumorType)) + 
  geom_point(size = 4, alpha = 0.7) +
  #geom_label(aes(tsne_x, tsne_y, label = location, alpha = 0.5)) +
  labs(x = "tSNE 1", y = "tSNE 2") +
  theme_bw(base_size = 20)


# plot UMAP
anno %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType, size = 1/estimate)) + 
  geom_point(alpha = 0.7) +
  #geom_label(aes(tsne_x, tsne_y, label = location, alpha = 0.5)) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_bw(base_size = 20)










