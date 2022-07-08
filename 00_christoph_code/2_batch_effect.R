# libraries & sources

library(tidyverse)
library(minfi)

source("./00_christoph_code/functions.R")

# load data and annotation 
betas <- readRDS("./00_christoph_data/methylation_data.rds")
anno <- readRDS("./00_christoph_annotation/sample_annotation.rds")

# check overlapping sample types between studies
anno %>% 
  group_by(source, tumorType) %>% 
  summarise(n = n()) %>% 
  arrange(tumorType)



### Investigation 1: Normal pancreatic tissue, DiDomenico vs. Jakel ------------

# keep only normal pancreatic tissue
anno_normal <- anno %>% 
  filter(tumorType == "normal" & location == "pancreas")

anno_normal %>% 
  group_by(source) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(source, n, fill = source)) +
  geom_col(width = 0.4) +
  theme_bw(base_size = 28) +
  labs(x = NULL, y = "Normal tissues (n)") +
  theme(legend.position = "none")

# list idats
idat <- list.files(path = "./00_christoph_data/idat/", pattern = "_Grn.idat")
file <- list.files(path = "./00_christoph_data/idat/", pattern = "_Grn.idat", full.names = TRUE)

file <- file[match(paste0(anno_normal$arrayId, "_Grn.idat"), idat)]
raw_normal <- read.metharray(basenames = file)
rm(idat, file)


# look into QC measures 1: bisulfite conversion efficiency 
conv_normal <- getControlBeta(raw_normal, controls = "BISULFITE CONVERSION II") %>% 
  as_tibble() %>% 
  group_by(arrayId) %>% 
  summarise(conversion = min(value[channel == "Red"]) / max(value[channel == "Green"]))

anno_normal <- anno_normal %>% 
  left_join(conv_normal)

# plot 
anno_normal %>% 
  ggplot(aes(source, conversion)) +
  geom_jitter(col = "magenta", size = 2, width = 0.1) +
  theme_bw(base_size = 24) +
  labs(x = NULL, y = "Control Score") +
  theme(legend.position = "none")

pp_normal <- preprocessNoob(raw_normal)
test <- getQ(pp_normal)



rm(raw_normal, conv_normal)

# compare methylation

betas_normal <- betas[, anno_normal$arrayId]

betas_normal_avg <- tibble(
  DiDomenico = apply(betas_normal[, anno_normal$source == "DiDomenico"], 1, mean, na.rm = TRUE), 
  Jakel = apply(betas_normal[, anno_normal$source == "Jakel"], 1, mean, na.rm = TRUE)
)
  
betas_normal_avg %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(value, fill = name)) +
  geom_density(alpha = 0.4) +
  theme_bw(base_size = 24) +
  labs(x = "Beta", y = "Density")

ttest <- genefilter::rowttests(x = betas_normal, fac = as.factor(anno_normal$source)) %>% 
  as_tibble(rownames = "probe")

ttest <- ttest %>% 
  mutate(d = betas_normal_avg$DiDomenico, 
         j = betas_normal_avg$Jakel)

ttest %>% 
  ggplot(aes(dm)) +
  geom_density(fill = "grey85", alpha = 0.6) +
  theme_bw(base_size = 24) +
  labs(x = "Delta (beta)", y = "Density")


ttest %>% 
  sample_n(size = 30000) %>% 
  ggplot(aes(d, -log10(p.value))) +
  geom_point(alpha = 0.2) +
  theme_bw(base_size = 24) +
  labs(x = "Avg. methylation (beta)", y = "-log10(p-value)")

ttest %>% 
  sample_n(size = 30000) %>% 
  ggplot(aes(d, dm)) +
  geom_point(fill = "grey85", alpha = 0.2) +
  theme_bw(base_size = 24) +
  labs(x = "Delta (beta)", y = "-log10(p-value")



### look at sample correlation heatmap ---------------------------------------

heatmap_anno <- anno_normal %>% 
  select(source, conversion) %>% 
  as.data.frame()
rownames(heatmap_anno) <- anno_normal$arrayId

cor(betas_normal) %>% 
  pheatmap::pheatmap(annotation_col = heatmap_anno,
                     labels_row = heatmap_anno$source, 
                     show_colnames = FALSE)



library(Rtsne)
