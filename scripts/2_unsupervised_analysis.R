# Christoph Geisenberger
# github: @cgeisenberger
# last edited 26/02/2023 by CG


# libraries

# Load required packages and sources
library(minfi)
library(RFpurify)

library(tidyverse)

library(Rtsne)
library(umap)

library(caret)

source("./scripts/0_helpers.R")



# load sample annotation and filter --------------------------------------------

anno <- readRDS("./input/sample_annotation.rds")

# keep primaries (PNET, PNEC, ACC, SPN, PB and PDAC) and normal pancreas (NORM)
anno <- anno %>% 
  filter(tumorType %in% c("PanNET", "ACC", "SPN", "PDAC", "NORM", "PanNEC", "PB")) %>% 
  filter(location %in% c("primary", "pancreas"))

# exclude MiNEN and biopsy samples
anno <- anno[c(1:269,271, 273, 275, 277, 279:326),]


# split into training and testing cohort ---------------------------------------

set.seed(12345678)
train_index <- createDataPartition(as.factor(anno$tumorType), 
                                   times = 1, list = TRUE)
train_groups <- rep("train", nrow(anno))
train_groups[train_index$Resample1] <- "test"

anno <- anno %>% 
  mutate(cohort = train_groups)
rm(train_index, train_groups)



# estimate tumor purity --------------------------------------------------------

# load unfiltered beta values
betas <- readRDS("./input/betas_pancreas_everything.rds")
betas <- betas[, anno$arrayId]

# get purity estimates 
absolute <- RFpurify::predict_purity_betas(betas = betas, method = "ABSOLUTE")
estimate <- RFpurify::predict_purity_betas(betas = betas, method = "ESTIMATE")

# collect data into tibble
anno <- anno %>% 
  mutate(absolute = absolute,
         estimate = estimate,
         avg_beta_unfiltered = apply(betas, 2, mean, na.rm = TRUE))

# add bisulfite conversion scores
conv_scores <- readRDS(file = "./annotation/sample_annotation_conversion_scores.csv")
anno <- left_join(anno, conv_scores)


# plot basic statistics for the data set ---------------------------------------

# overview of all variables
anno %>% 
  select(-c(arrayId, source, sampleName, location)) %>% 
  GGally::ggpairs()


# number of cases
anno %>% 
  group_by(tumorType) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(tumorType, n, fill = tumorType)) +
  geom_col() +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Number of Cases")

# number of cases by source
anno %>% 
  ggplot(aes(tumorType)) +
  geom_bar(aes(fill = source)) +
  labs(title="Tumor Type By Source",x = "", y = "No. of cases") +
  theme_bw(base_size = 18)

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
  geom_point(position = position_jitterdodge(0.5), alpha=0.7) +
  #geom_jitter(aes(fill=tumorType, shape=tumorType)) +
  scale_shape_manual(values = c(21, 25, 22, 23, 24, 8, 10)) +
  theme_classic(base_size = 18)

anno %>% 
  ggplot(aes(tumorType, absolute, col = tumorType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(0.5), alpha=0.9) +
  #geom_jitter(aes(fill=tumorType, shape=tumorType)) +
  scale_shape_manual(values = c(21, 25, 22, 23, 24, 8, 1)) +
  theme_classic(base_size = 18)

# average methylation
anno %>% 
  ggplot(aes(tumorType, avg_beta_unfiltered, col = tumorType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.7) +
  scale_shape_manual(values = c(21, 25, 22, 23, 24, 8, 1)) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# clean up
rm(absolute, estimate, betas)



### dimensionality reduction ---------------------------------------------------

# load filtered beta values
betas <- readRDS("./input/betas_pancreas_filtered.rds")
betas <- betas[, anno$arrayId]
any(is.na(betas))

# add avg. methylation to annotation
anno <- anno %>% 
  add_column(avg_beta_filtered = apply(betas, 2, mean, na.rm = TRUE))

# determine most variable probes across dataset and subset beta values
probe_var <- apply(betas[, anno$cohort == "train"], 1, var)
probes_topvar <- rownames(betas)[order(probe_var, decreasing = TRUE)]
#saveRDS(object = probes_topvar, file = "./output/pancreas_top_variable_probes_training_set_0112023.rds")

# pick betas for 5,000 top variable probes
betas_topvar <- betas[probes_topvar[1:5000], ]

# run UMAP
set.seed(45098)
umap_settings <- umap.defaults
umap_settings$n_neighbors = 15
umap_settings$min_dist = 0.2

umap <- umap(d = t(betas_topvar), config = umap_settings, ret_model = TRUE)

anno <- anno %>% 
  mutate(umap_x = umap$layout[, 1], 
         umap_y = umap$layout[, 2])

saveRDS(object = umap, file = "./output/umap_model_01112023.rds")
saveRDS(object = anno, file = "./output/sample_annotation_umap_purity_01112023.rds")


# plot UMAP
anno %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType)) +
  geom_point(size = 4) +
  #geom_text(aes(label = sampleName), size = 4) +
  theme_classic(base_size = 24) +
  theme(legend.direction = "horizontal", legend.position = "bottom") +
  labs(x = "UMAP 1", y = "UMAP 2")


anno %>% 
  ggplot(aes(umap_x, umap_y, col = conversion)) +
  geom_point(size = 5) +
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
  geom_point(size = 5) +
  paletteer::scale_color_paletteer_c("grDevices::Blue-Red 2") +
  theme_classic(base_size = 24) +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks = element_blank(), 
        legend.direction = "horizontal", 
        legend.position = "bottom", 
        legend.key.width=unit(3,"cm")) +
  labs(x = NULL, y = NULL)

