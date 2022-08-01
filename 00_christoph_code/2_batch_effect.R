# libraries & sources

library(tidyverse)
library(minfi)
library(paletteer)

source("./00_christoph_code/functions.R")


# create custom color scale ------------------------------------
branded_colors <- c("#f6d8ae", "#2e4057", "#da4167", "#3cdbd3", "#f4d35e", "#8d96a3")


# load data and annotation 
betas <- readRDS("./00_christoph_data/methylation_data.rds")
anno <- readRDS("./00_christoph_annotation/sample_annotation.rds")

# check overlapping sample types between studies
anno %>% 
  filter(source != "UMCU") %>%
  filter(tumorType != "MACNEC") %>%
  #filter(location %in% c("primary", "acc normal", "pancreas")) %>% 
  group_by(tumorType, source, location) %>% 
  summarise(n = n()) %>% 
  arrange(tumorType) #%>% 
  pull(n) %>% sum

  
  
# n = 233 incl UMCU
# n = 215 excl. UMCU

### Investigation 1: Normal pancreatic tissue, DiDomenico vs. Jakel ------------

# keep only normal pancreatic tissue
anno_normal <- anno %>% 
  filter(tumorType == "normal" & location == "pancreas")

anno_normal %>% 
  group_by(source) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(source, n, fill = source)) +
  geom_col(width = 0.4) +
  labs(x = NULL, y = "Normal tissues (n)") +
  theme_bw(base_size = 28) +
  theme(legend.position = "none") +
  scale_fill_manual(values = branded_colors) +
  ylim(c(0, 15))

# list idats
idat <- list.files(path = "./00_christoph_data/idat/", pattern = "_Grn.idat")
file <- list.files(path = "./00_christoph_data/idat/", pattern = "_Grn.idat", full.names = TRUE)

file <- file[match(paste0(anno_normal$arrayId, "_Grn.idat"), idat)]
raw_normal <- read.metharray(basenames = file)
rm(idat, file)

# investigate QC measures 1: bisulfite conversion efficiency 
conv_normal <- getControlBeta(raw_normal, controls = "BISULFITE CONVERSION II") %>% 
  as_tibble() %>% 
  group_by(arrayId) %>% 
  summarise(conversion = min(value[channel == "Red"]) / max(value[channel == "Green"]))

anno_normal <- anno_normal %>% 
  left_join(conv_normal)

# plot 
anno_normal %>% 
  ggplot(aes(source, conversion, fill = source)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(size = 2, width = 0.1) +
  theme_bw(base_size = 24) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Control Score") +
  scale_fill_manual(values = branded_colors)

rm(raw_normal, conv_normal)


# compare methylation

betas_normal <- betas[, anno_normal$arrayId]

betas_normal_avg <- tibble(
  DiDomenico = apply(betas_normal[, anno_normal$source == "DiDomenico"], 1, mean, na.rm = TRUE), 
  Jakel = apply(betas_normal[, anno_normal$source == "Jakel"], 1, mean, na.rm = TRUE)
)
  
betas_normal_avg %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(name, value, fill = name)) +
  geom_violin(alpha = 0.8) +
  theme_bw(base_size = 24) +
  labs(x = "Beta", y = "Density") +
  scale_fill_manual(values = branded_colors) +
  theme(legend.position = "none")


betas_normal_avg %>% 
  sample_n(size = 30000) %>% 
  ggplot(aes(DiDomenico, Jakel)) +
  geom_bin2d(bins = 100) +
  paletteer::scale_fill_paletteer_c("scico::tokyo") +
  theme_bw(base_size = 28) +
  geom_abline(slope = 1, intercept = 0, col = "grey", lty = 2) +
  theme(legend.position = "none")




# perform t-test

ttest <- genefilter::rowttests(x = betas_normal, fac = as.factor(anno_normal$source)) %>% 
  as_tibble(rownames = "probe")

ttest <- ttest %>% 
  mutate(d = betas_normal_avg$DiDomenico, 
         j = betas_normal_avg$Jakel)


# estimate FDR for p < 0.001
ttest_perm <- genefilter::rowttests(x = betas_normal, fac = as.factor(sample(anno_normal$source))) %>% 
  as_tibble(rownames = "probe")
p_cutoff <- 0.001
sum(ttest_perm$p.value < p_cutoff)/sum(ttest$p.value < p_cutoff) # 0.01806834



# volcano plot

ttest %>%
  sample_n(size = 30000) %>% 
  ggplot(aes(dm, -log10(p.value))) +
  geom_point(alpha = 0.4, col = branded_colors[6]) +
  theme_bw(base_size = 24) +
  labs(x = "Methylation difference (beta)", y = "-log10(p-value)")


# show that significant beta values are enriched for low methylation 
  
ttest %>% 
  filter(p.value < 0.001) %>% 
  sample_n(size = 30000) %>% 
  ggplot(aes(j)) +
  geom_density(alpha = 0.2, fill = branded_colors[3]) +
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



### Investigation 2: PanNets, DiDomenico vs. Jakel vs. Chan --------------------


# select (1) primary (2) PanNETs which are (3) not from UMCU
anno_pnet <- anno %>% 
  filter(tumorType == "PanNET" & location == "primary" & source != "UMCU")

# plot sample numbers per study
anno_pnet %>% 
  group_by(source) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(source, n, fill = source)) +
  geom_col(width = 0.4) +
  labs(x = NULL, y = "Normal tissues (n)") +
  theme_bw(base_size = 28) +
  theme(legend.position = "none") +
  scale_fill_manual(values = branded_colors)

# collect beta values, add avg. beta per sample to annotation
betas_pnet <- betas[, anno_pnet$arrayId]
anno_pnet <- anno_pnet %>% 
  mutate(mean_beta = apply(betas_pnet, 2, mean, na.rm = TRUE))

anno_pnet %>% 
  ggplot(aes(source, mean_beta, fill = source)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(values = branded_colors) +
  theme_bw(base_size = 24)

# t-SNE of correlation matrix
cor_pnet <- cor(betas_pnet)

umap_pnet <- umap::umap(d = cor_pnet)
umap_pnet <- umap::umap(d = t(betas_pnet[sample(1:nrow(betas_pnet), size = 5000), ]))


plot(umap_pnet$layout)

anno_pnet <- anno_pnet %>% 
  mutate(umap_x = umap_pnet$layout[, 1], 
         umap_y = umap_pnet$layout[, 2])

anno_pnet %>% 
  ggplot(aes(umap_x, umap_y, col = source, alpha = conversion)) +
  #geom_point(size = 3) +
  geom_label(label = anno_pnet$sampleName) +
  scale_color_manual(values = branded_colors) +
  theme_bw(base_size = 24) +
  labs(x = "t-SNE 1", y = "t-SNE 2")



pheatmap::pheatmap(cor_pnet, labels_row = anno_pnet$source)
tsne_pnet <- Rtsne::Rtsne(t(betas_pnet[sample(1:nrow(betas_pnet), size = 5000), ]))

anno_pnet <- anno_pnet %>% 
  mutate(tsne_x = tsne_pnet$Y[, 1], 
         tsne_y = tsne_pnet$Y[, 2])

anno_pnet %>% 
  ggplot(aes(tsne_x, tsne_y, col = source, alpha = conversion)) +
  geom_point(size = 3) +
  scale_color_manual(values = branded_colors) +
  theme_bw(base_size = 24) +
  labs(x = "t-SNE 1", y = "t-SNE 2")


# list idats, read raw data
idat_pnet <- list.files(path = "./00_christoph_data/idat/", pattern = "_Grn.idat")
file_pnet <- list.files(path = "./00_christoph_data/idat/", pattern = "_Grn.idat", full.names = TRUE)
file_pnet <- file_pnet[match(paste0(anno_pnet$arrayId, "_Grn.idat"), idat_pnet)]
raw_pnet <- read.metharray(basenames = file_pnet) # 95 samples, all 450k
rm(idat_pnet, file_pnet)


mean_beta_pnet <- tibble(
  Chan = apply(betas_pnet[, anno_pnet$source == "Chan"], 1, mean, na.rm = TRUE),
  DiDomenico = apply(betas_pnet[, anno_pnet$source == "DiDomenico"], 1, mean, na.rm = TRUE), 
  Jakel = apply(betas_pnet[, anno_pnet$source == "Jakel"], 1, mean, na.rm = TRUE)
)


conv_pnet <- getControlBeta(raw_pnet, controls = "BISULFITE CONVERSION II") %>% 
  as_tibble() %>% 
  group_by(arrayId) %>% 
  summarise(conversion = min(value[channel == "Red"]) / max(value[channel == "Green"]))

anno_pnet <- anno_pnet %>% 
  left_join(conv_pnet)


anno_pnet %>% 
  ggplot(aes(source, conversion )) +
  geom_jitter()

# look at beta value distribution

betas_pnet <- betas[, anno_pnet$arrayId]

anno_pnet <- anno_pnet %>% 
  mutate(mean_beta = apply(betas_pnet, 2, mean, na.rm = TRUE))

anno_pnet %>% 
  ggplot(aes(source, mean_beta)) +
  geom_boxplot()

mean_beta_pnet <- tibble(
  Chan = apply(betas_pnet[, anno_pnet$source == "Chan"], 1, mean, na.rm = TRUE),
  DiDomenico = apply(betas_pnet[, anno_pnet$source == "DiDomenico"], 1, mean, na.rm = TRUE), 
  Jakel = apply(betas_pnet[, anno_pnet$source == "Jakel"], 1, mean, na.rm = TRUE)
)

# density plot of beta value distribution
mean_beta_pnet %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(value, fill = name)) +
  geom_density(alpha = 0.8) +
  theme_bw(base_size = 24) +
  labs(x = "Beta", y = "Density") +
  scale_fill_manual(values = branded_colors) +
  theme(legend.position = "none") +
  facet_grid(rows = vars(name))


mean_beta_pnet %>% 
  sample_n(size = 30000) %>% 
  ggplot(aes(Chan, DiDomenico)) +
  geom_point() +
  theme_bw(base_size = 28) +
  geom_abline(slope = 1, intercept = 0, col = "grey", lty = 2) +
  theme(legend.position = "none")


# investigate QC measures

oob <- getOOB(raw_pnet)
str(oob)

test <- cor(oob$Red)
test %>% pheatmap::pheatmap(labels_row = anno_pnet$source)

getQC( mapToGenome(raw_pnet))

