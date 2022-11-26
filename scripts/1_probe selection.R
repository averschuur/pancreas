Selection MVPs

### preparation ----------------------------------------------------------------
# libraries
library(RFpurify)
library(minfi)
library(minfiData)
library(tidyverse)
library(Rtsne)
library(umap)

branded_colors <- c("#2e4057", "#da4167", "#3cdbd3", "#8d96a3", "#f4d35e", "#f6d8ae", "#fff8ae")

# load sample annotation and filter --------------------------------------------
betas <- readRDS("./data/methylation_data_filtered.rds")
betas <- betas[, anno$array_id]
any(is.na(betas))

anno <- readRDS("./data/sample_annotation.rds")

anno <- anno %>% 
  filter(tumorType %in% c("PanNET", "ACC", "SPN", "PDAC", "normal", "acc normal")) %>% 
  filter(location %in% c("primary", "pancreas")) %>% 
  filter(source != "UMCU")

betas <- betas[,anno$arrayId]

# add avg. methylation to annotation
anno <- anno %>% 
  add_column(avg_beta = apply(betas, 2, mean, na.rm = TRUE))

# determine 100, 200, 500, 1000, 2000, 5000, 10000 most variable probes across dataset and subset beta values
probe_var <- apply(betas, 1, var)
probes_topvar <- order(probe_var, decreasing = TRUE)[1:20000]
betas_topvar <- betas[probes_topvar, ]
saveRDS(object = betas_topvar, file = "./data/betas_filtered_20000topvar.rds")

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
