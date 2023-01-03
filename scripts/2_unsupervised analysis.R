# Christoph Geisenberger
# github: @cgeisenberger
# last edited 03/01/2023 by AV Verschuur


# libraries

# Load required packages and sources
library(RFpurify)
library(minfi)
library(minfiData)
library(tidyverse)
library(Rtsne)
library(umap)

source("./scripts/0_branded_colors.R")
source("./scripts/0_functions.R")

### load sample annotation and filter ------------------------------------------

anno <- readRDS("./data/sample_annotation.rds")

# keep primaries (PanNET, ACC, SPN and PDAC) and normal pancreas, filter out UMC
anno <- anno %>% 
  filter(tumorType %in% c("PanNET", "ACC", "SPN", "PDAC", "normal", "acc normal", "PanNEC", "PB")) %>% 
  filter(location %in% c("primary", "pancreas")) %>% 
  filter(source != "UMCU")


### estimate tumor purity ------------------------------------------------------

# load unfiltered beta values
betas <- readRDS("./input/pancreas_betas_everything.rds")
betas <- betas[,anno$arrayId]

# get purity estimates 
absolute <- RFpurify::predict_purity_betas(betas = betas, method = "ABSOLUTE")
estimate <- RFpurify::predict_purity_betas(betas = betas, method = "ESTIMATE")

# collect data into tibble
anno <- anno %>% 
  mutate(absolute = absolute,
         estimate = estimate,
         avg_beta = apply(betas, 2, mean, na.rm = TRUE))

# add bisulfite conversion scores

conv_scores <- read_csv(file = "./annotation/sample_annotation_conversion_scores.csv")
conv_scores <- conv_scores %>%
  select(arrayId, conversion)

anno <- left_join(anno, conv_scores)


### plot basic statistics for the data set -------------------------------------

# number of cases
anno %>% 
  group_by(tumorType) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(tumorType, n, fill = tumorType)) +
  geom_col() +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  scale_fill_manual(values = branded_colors) +
  labs(x = NULL, y = "Number of Cases")

# number of cases by source
anno %>% 
  ggplot(aes(tumorType)) +
  geom_bar(aes(fill = source)) +
  labs(title="Tumor Type By Source",x = "", y = "No. of cases") +
  scale_fill_manual(values = branded_colors) +
  theme_bw(base_size = 18) +
ggsave("Figure 1A_count tumorType per source.png", path= "./output/")

# tumor purity
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
  scale_shape_manual(values = c(21, 25, 22, 23, 24, 8, 10)) +
  scale_colour_manual(values = branded_colors) +
  theme_classic(base_size = 18)
ggsave("Figure 1D_tumor purity per tumorType_estimate.png", path= "./output/")

anno %>% 
  ggplot(aes(tumorType, absolute, col = tumorType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(shape=tumorType), position = position_jitterdodge(0.5), alpha=0.9) +
  #geom_jitter(aes(fill=tumorType, shape=tumorType)) +
  scale_shape_manual(values = c(21, 25, 22, 23, 24, 8, 1)) +
  scale_colour_manual(values = branded_colors) +
  theme_classic(base_size = 18)
ggsave("Figure 1D_tumor purity per tumorType_absolute.png", path= "./output/")

# average methylation
anno %>% 
  ggplot(aes(tumorType, avg_beta, col = tumorType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_shape_manual(values = c(21, 25, 22, 23, 24, 8, 1)) +
  scale_colour_manual(values = branded_colors) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("Figure 1C_average methylation.png", path= "./output/")

# clean up
rm(absolute, estimate, betas)


### dimensionality reduction ---------------------------------------------------

betas <- readRDS("./input/pancreas_betas_filtered.rds")
betas <- betas[, anno$arrayId]
any(is.na(betas))

# add avg. methylation to annotation
anno <- anno %>% 
  add_column(avg_beta = apply(betas, 2, mean, na.rm = TRUE))

# determine most variable probes across dataset and subset beta values
probe_var <- apply(betas, 1, var)
probes_topvar <- rownames(betas)[order(probe_var, decreasing = TRUE)[1:10000]]
betas_topvar <- betas[probes_topvar, ]

saveRDS(object = probes_topvar, file = "./input/pancreas_top_variable_probes.rds")

# plot UMAP
anno %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType)) +
  geom_point(size = 4) +
  theme_classic(base_size = 24) +
  scale_color_manual(values = branded_colors1) +
  theme(legend.direction = "horizontal", legend.position = "bottom") +
  labs(x = "UMAP 1", y = "UMAP 2")


anno %>% 
  ggplot(aes(umap_x, umap_y, col = conversion)) +
  geom_point(size = 3) +
  paletteer::scale_color_paletteer_c("grDevices::Blue-Red 2") +
  theme_classic(base_size = 24) +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks = element_blank(), 
        legend.direction = "horizontal", 
        legend.position = "bottom", 
        legend.key.width=unit(3,"cm")) +
  labs(x = NULL, y = NULL)

anno %>% 
  ggplot(aes(umap_x, umap_y, col = absolute)) +
  geom_point(size = 3) +
  paletteer::scale_color_paletteer_c("grDevices::Blue-Red 2") +
  theme_classic(base_size = 24) +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks = element_blank(), 
        legend.direction = "horizontal", 
        legend.position = "bottom", 
        legend.key.width=unit(3,"cm")) +
  labs(x = NULL, y = NULL)





