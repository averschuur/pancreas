### preparation ----------------------------------------------------------------

# Load required packages and sources
library(RFpurify)
library(minfi)
library(minfiData)
library(tidyverse)
library(Rtsne)
library(umap)

source("./scripts/branded_colors.R")

# load sample annotation and filter --------------------------------------------

anno <- readRDS("./data/sample_annotation.rds")

# keep primaries (PanNET, ACC, SPN and PDAC) and normal pancreas, filter out UMC
anno <- anno %>% 
  filter(tumorType %in% c("PanNET", "ACC", "SPN", "PDAC", "normal", "acc normal", "PanNEC")) %>% 
  filter(location %in% c("primary", "pancreas")) %>% 
  filter(source != "UMCU")

### estimate tumor purity ------------------------------------------------------

# load unfiltered beta values
betas <- readRDS("./data/methylation_data.rds")
betas <- betas[,anno$arrayId]

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
  ggplot(aes(tumorType)) +
  geom_bar(aes(fill = source)) +
  scale_fill_manual(values = branded_colors) +
  theme_bw(base_size = 18)
ggsave("Figure 1A_count tumorType per source.png", path= "./output/")

anno %>% 
  ggplot(aes(absolute, estimate)) +
  geom_point(size = 3) +
  geom_smooth() +
  theme_bw(base_size = 18) +
  labs(x = "Purity (ABSOLUTE)", y = "Purity (ESTIMATE)")

anno %>% 
  ggplot(aes(tumorType, estimate, col = tumorType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(shape=tumorType), position = position_jitterdodge(0.5), alpha=0.9) +
  #geom_jitter(aes(fill=tumorType, shape=tumorType)) +
  scale_shape_manual(values = c(21, 25, 22, 23, 24, 8)) +
  scale_colour_manual(values = branded_colors) +
  theme_classic(base_size = 18)
ggsave("Figure 1D_tumor purity per tumorType_estimate.png", path= "./output/")

anno %>% 
  ggplot(aes(tumorType, absolute, col = tumorType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(shape=tumorType), position = position_jitterdodge(0.5), alpha=0.9) +
  #geom_jitter(aes(fill=tumorType, shape=tumorType)) +
  scale_shape_manual(values = c(21, 25, 22, 23, 24, 8)) +
  scale_colour_manual(values = branded_colors) +
  theme_classic(base_size = 18)
ggsave("Figure 1D_tumor purity per tumorType_absolute.png", path= "./output/")

anno %>% 
  ggplot(aes(tumorType, avg_beta, col = tumorType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_shape_manual(values = c(21, 25, 22, 23, 24, 8)) +
  scale_colour_manual(values = branded_colors) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("Figure 1C_average methylation.png", path= "./output/")

# clean up
rm(absolute, estimate, betas)



### unsupervised analysis: tSNE & UMAP ----------------------------------

# select most variable probes
most_var <- apply(betas, 1, var)
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
  #geom_text(data = anno %>% 
                    #filter(tumorType == "PanNET"), aes(tsne_x, tsne_y, label = sampleName)) +
  #geom_label(aes(tsne_x, tsne_y, label = sampleName, alpha = 0.5)) +
  scale_colour_manual(values = branded_colors) +
  labs(x = "tSNE 1", y = "tSNE 2") +
  theme_bw(base_size = 20)
ggsave("tSNE-all.png", path= "./output/")

# plot UMAP
anno %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType)) + 
  geom_point(alpha = 0.7) +
  #geom_label(aes(tsne_x, tsne_y, label = sampleName, alpha = 0.5)) +
  scale_colour_manual(values = branded_colors) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_bw(base_size = 20)
ggsave("UMAP-all.png", path= "./output/")



# perform dimensionlity reduction ----------------------------------------------


betas <- readRDS("./input/betas_filtered.rds")
betas <- betas[, anno$array_id]
any(is.na(betas))

# add avg. methylation to annotation
anno <- anno %>% 
  add_column(avg_beta = apply(betas, 2, mean, na.rm = TRUE))

# determine 5000 most variable probes across dataset and subset beta values
probe_var <- apply(betas, 1, var)
probes_topvar <- order(probe_var, decreasing = TRUE)[1:10000]
betas_topvar <- betas[probes_topvar, ]
saveRDS(object = betas_topvar, file = "./data/betas_filtered_topvar.rds")

# run UMAP
umap <- umap(d = t(betas_topvar), ret_model = TRUE)
#saveRDS(object = umap, file = "./input/umap_model.rds")

anno <- anno %>% 
  mutate(umap_x = umap$layout[, 1], 
         umap_y = umap$layout[, 2])
#saveRDS(object = anno, file = "./input/sample_annotation_extended.rds")


# plot UMAP
anno %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType)) +
  geom_point(size = 3) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = branded_colors) +
  labs(x = "UMAP 1", y = "UMAP 2")
ggsave("UMAP-all.png", path= "./output/")


### plot basic statistics for the data set -------------------------------------

anno %>% 
  group_by(tumorType) %>% 
  summarise(n = n()) %>% 
ggplot(aes(tumorType, n, fill = tumorType)) +
  geom_col() +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  scale_fill_manual(values = branded_colors) +
  labs(x = NULL, y = "Number of Cases")


# average methylation
anno %>% 
  ggplot(aes(tumorType, avg_beta, col = tumorType)) +
  geom_boxplot(lwd = 1, outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  scale_colour_manual(values = branded_colors) +
  labs(x = NULL, y = "Average Methylation (beta)")


# tumor purity
anno %>% 
  ggplot(aes(tumorType, absolute, col = tumorType)) +
  geom_boxplot(lwd = 1, outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  scale_colour_manual(values = branded_colors) +
  labs(x = NULL, y = "Tumor Purity")





### TEST GROUNDS --------------

umap2 <- umap(d = most_var)
plot(umap2$layout[, 1], umap2$layout[, 2])
test <- apply(most_var, 1, function(x) tapply(x, INDEX = anno$tumorType, FUN = mean, simplify = TRUE)) %>% t %>% as_tibble
test <- test %>% mutate(
  umap_x = umap2$layout[, 1], 
  umap_y = umap2$layout[, 2]
)

test %>% 
  pivot_longer(cols = 1:5, names_to = "label", values_to = "beta") %>% 
  filter(!label == "NORM") %>% 
  ggplot(aes(umap_x, umap_y, col = beta)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_classic(base_size = 20) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  facet_wrap(facets = vars(label))





