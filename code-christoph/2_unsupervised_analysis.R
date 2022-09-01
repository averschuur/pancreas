# Christoph Geisenberger
# github: @cgeisenberger
# last edited 28/08/2022


# libraries
library(RFpurify)
library(minfi)
library(minfiData)
library(tidyverse)
library(Rtsne)
library(umap)

source("./code/branded_colors.R")


# load sample annotation and filter --------------------------------------------

anno <- readRDS("./input/sample_annotation.rds")

# keep primaries (PanNET, ACC, SPN and PDAC) and normal pancreas, filter out UMC
anno <- anno %>% 
  filter(label %in% c("PanNET", "ACC", "SPN", "PDAC", "NORM")) %>% 
  filter(location %in% c("primary", "pancreas")) %>% 
  filter(source != "UMCU")



### estimate tumor purity ------------------------------------------------------

# load unfiltered beta values
betas <- readRDS("./input/betas.rds")
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



# perform dimensionlity reduction ----------------------------------------------


betas <- readRDS("./input/betas_probes_filtered.rds")
betas <- betas[, anno$array_id]
any(is.na(betas))

# add avg. methylation to annotation
anno <- anno %>% 
  add_column(avg_beta = apply(betas, 2, mean, na.rm = TRUE))

# determine 5000 most variable probes across dataset and subset beta values
probe_var <- apply(betas, 1, var)
probes_topvar <- order(probe_var, decreasing = TRUE)[1:10000]
betas_topvar <- betas[probes_topvar, ]
saveRDS(object = betas_topvar, file = "./input/betas_filtered_topvar.rds")

# run UMAP
umap <- umap(d = t(betas_topvar))

anno <- anno %>% 
  mutate(umap_x = umap$layout[, 1], 
         umap_y = umap$layout[, 2])

# plot UMAP
anno %>% 
  ggplot(aes(umap_x, umap_y, col = label)) +
  geom_point(size = 3) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = branded_colors) +
  labs(x = "UMAP 1", y = "UMAP 2")


### plot basic statistics for the data set -------------------------------------

anno %>% 
  group_by(label) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(label, n, fill = label)) +
  geom_col() +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  scale_fill_manual(values = branded_colors) +
  labs(x = NULL, y = "Number of Cases")
  

# average methylation
anno %>% 
  ggplot(aes(label, avg_beta, col = label)) +
  geom_boxplot(lwd = 1, outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  scale_colour_manual(values = branded_colors) +
  labs(x = NULL, y = "Average Methylation (beta)")


# tumor purity
anno %>% 
  ggplot(aes(label, absolute, col = label)) +
  geom_boxplot(lwd = 1, outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  scale_colour_manual(values = branded_colors) +
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
  pivot_longer(cols = 1:5, names_to = "label", values_to = "beta") %>% 
  filter(!label == "NORM") %>% 
  ggplot(aes(umap_x, umap_y, col = beta)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_classic(base_size = 20) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  facet_wrap(facets = vars(label))

