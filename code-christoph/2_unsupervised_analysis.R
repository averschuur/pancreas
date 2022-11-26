# Christoph Geisenberger
# github: @cgeisenberger
# last edited 21/11/2022


# libraries
library(RFpurify)
library(minfi)
library(tidyverse)
library(Rtsne)
library(umap)

source("./code-christoph/0_branded_colors.R")
source("./code-christoph/0_functions.R")


# load sample annotation and filter --------------------------------------------

anno <- readRDS("./input/sample_annotation.rds")

# keep primaries (PanNEC, PanNET, ACC, SPN and PDAC) and normal pancreas w?o UMC
anno <- anno %>% 
  filter(label %in% c("PanNEC", "PanNET", "PB", "ACC", "SPN", "PDAC", "NORM")) %>% 
  filter(location %in% c("primary", "pancreas")) %>% 
  filter(source != "UMCU")


anno %>% 
  group_by(label, source) %>% 
  summarise(n = n())



### estimate tumor purity ------------------------------------------------------

# load unfiltered beta values
betas <- readRDS("./input/pancreas_betas_everything.rds")
betas <- betas[, anno$array_id]

# get purity estimates 
absolute <- RFpurify::predict_purity_betas(betas = betas, method = "ABSOLUTE")
estimate <- RFpurify::predict_purity_betas(betas = betas, method = "ESTIMATE")

# collect data into tibble
anno <- anno %>% 
  mutate(absolute = absolute,
         estimate = estimate)

# clean up
rm(absolute, estimate, betas)



# dimensionality reduction -----------------------------------------------------

betas <- readRDS("./input/pancreas_betas_filtered.rds")
betas <- betas[, anno$array_id]
any(is.na(betas))

# add avg. methylation to annotation
anno <- anno %>% 
  add_column(avg_beta = apply(betas, 2, mean, na.rm = TRUE))

# determine 5000 most variable probes across dataset and subset beta values
probe_var <- apply(betas, 1, var)
probes_topvar <- rownames(betas)[order(probe_var, decreasing = TRUE)[1:5000]]
betas_topvar <- betas[probes_topvar, ]

saveRDS(object = probes_topvar, file = "./input/pancreas_5k_top_variable_probes.rds")

# run UMAP
umap <- umap(d = t(betas_topvar), ret_model = TRUE)
saveRDS(object = umap, file = "./input/umap_model.rds")

anno <- anno %>% 
  mutate(umap_x = umap$layout[, 1], 
         umap_y = umap$layout[, 2])
saveRDS(object = anno, file = "./input/sample_annotation_umap_purity.rds")


# plot UMAP
anno %>% 
  ggplot(aes(umap_x, umap_y, col = label)) +
  geom_point(size = 4) +
  theme_classic(base_size = 24) +
  scale_color_manual(values = branded_colors1) +
  theme(legend.direction = "horizontal") +
  labs(x = "UMAP 1", y = "UMAP 2")


anno %>% 
  ggplot(aes(umap_x, umap_y, col = avg_beta)) +
  geom_point(size = 8) +
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
  geom_point(size = 8) +
  paletteer::scale_color_paletteer_c("grDevices::Blue-Red 2") +
  theme_classic(base_size = 24) +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks = element_blank(), 
        legend.direction = "horizontal", 
        legend.position = "bottom", 
        legend.key.width=unit(3,"cm")) +
  labs(x = NULL, y = NULL)



### plot basic statistics for the data set -------------------------------------

anno %>% 
  group_by(label) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(label, n, fill = label)) +
  geom_col() +
  theme_classic(base_size = 24) +
  theme(legend.position = "none") +
  scale_fill_manual(values = branded_colors1) +
  labs(x = NULL, y = "Number of Cases")
  

# average methylation
anno %>% 
  ggplot(aes(label, avg_beta, col = label)) +
  geom_boxplot(lwd = 1, outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none") +
  scale_colour_manual(values = branded_colors1) +
  labs(x = NULL, y = "Average Methylation (beta)")


# tumor purity
anno %>% 
  ggplot(aes(label, absolute, col = label)) +
  geom_boxplot(lwd = 1, outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none") +
  scale_colour_manual(values = branded_colors1) +
  labs(x = NULL, y = "Tumor Purity")





### TEST GROUNDS --------------

umap2 <- umap(d = betas_topvar)
plot(umap2$layout[, 1], umap2$layout[, 2])
test <- apply(betas_topvar, 1, function(x) tapply(x, INDEX = anno$label, FUN = mean, simplify = TRUE)) %>% t %>% as_tibble
test <- test %>% mutate(
  umap_x = umap2$layout[, 1], 
  umap_y = umap2$layout[, 2]
)

test %>% 
  pivot_longer(cols = 1:7, names_to = "label", values_to = "beta") %>% 
  filter(!label == "NORM") %>% 
  ggplot(aes(umap_x, umap_y, col = beta)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_classic(base_size = 20) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  facet_wrap(facets = vars(label))

