### Load required packages and sources  ----------------------------------------------

library(minfi)
library(RFpurify)

library(tidyverse)

library(Rtsne)
library(umap)

source("./scripts/0_helpers.R")



#### load data and annotation 
anno <- readRDS("./input/sample_annotation.rds")
betas <- readRDS(file = "./input/betas_pancreas_everything.rds")


# check overlapping sample types between studies
anno <- anno %>% 
  filter(tumorType %in% c("PB", "ACC", "PanNEC", "PanNET", "NORM", "PDAC", "SPN", "Mixed")) %>%
  filter(location %in% c("primary", "pancreas")) 

anno %>%
  group_by(tumorType, source, location) %>% 
  summarise(n = n()) %>% 
  arrange(tumorType) %>% 
  pull(n) %>% sum

#n= 204 excl UMCU and MACNEC
  

### Investigation 1: Normal pancreatic tissue, DiDomenico vs. Jakel ------------

# keep only normal pancreatic tissue
anno_normal <- anno %>% 
  filter(tumorType == "NORM")

anno_normal %>% 
  group_by(source) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(source, n, fill = source)) +
  geom_col(width = 0.4) +
  labs(x = NULL, y = "Normal tissues (n)") +
  theme_bw(base_size = 28) +
  theme(legend.position = "none") +
  scale_fill_manual(values = branded_colors1) +
  ylim(c(0, 30))
ggsave("BatchEffect_Normal Pancreas_No.png", path= "./plots/", dpi = 500)

# list idats
idat <- list.files(path = "./data/ALL IDATS", pattern = "_Grn.idat")
file <- list.files(path = "./data/ALL IDATS", pattern = "_Grn.idat", full.names = TRUE)

file <- file[match(paste0(anno_normal$arrayId, "_Grn.idat"), idat)]
raw_normal <- read.metharray(basenames = file, force = TRUE)
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
  scale_fill_manual(values = branded_colors1)
ggsave("BatchEffect_Normal Pancreas_Controle Score.png", path= "./output/")

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
  scale_fill_manual(values = branded_colors1) +
  theme(legend.position = "none")
ggsave("BatchEffect_Normal Pancreas_Violin_Beta's avg.png", path= "./output/")

betas_normal_avg %>% 
  sample_n(size = 30000) %>% 
  ggplot(aes(DiDomenico, Jakel)) +
  geom_bin2d(bins = 100) +
  paletteer::scale_fill_paletteer_c("scico::tokyo") +
  theme_bw(base_size = 28) +
  geom_abline(slope = 1, intercept = 0, col = "grey", lty = 2) +
  theme(legend.position = "none")
ggsave("BatchEffect_Normal Pancreas_Beta's avg.png", path= "./output/")


### perform t-test

ttest <- genefilter::rowttests(x = betas_normal, fac = as.factor(anno_normal$source)) %>% 
  as_tibble(rownames = "probe")

ttest <- ttest %>% 
  mutate(d = betas_normal_avg$DiDomenico, 
         j = betas_normal_avg$Jakel)


# estimate FDR for p < 0.001
ttest_perm <- genefilter::rowttests(x = betas_normal, fac = as.factor(sample(anno_normal$source))) %>% 
  as_tibble(rownames = "probe")
p_cutoff <- 0.001
sum(ttest_perm$p.value < p_cutoff)/sum(ttest$p.value < p_cutoff)
# 0.0006957974


# volcano plot

ttest %>%
  sample_n(size = 30000) %>% 
  ggplot(aes(dm, -log10(p.value))) +
  geom_point(alpha = 0.4, col = branded_colors1[6]) +
  theme_bw(base_size = 24) +
  labs(x = "Methylation difference (beta)", y = "-log10(p-value)")
ggsave("BatchEffect_Normal Pancreas_Volcano_T-Testg.png", path= "./output/")

# show that significant beta values are enriched for low methylation 
  
ttest %>% 
  filter(p.value < 0.001) %>% 
  sample_n(size = 30000) %>% 
  ggplot(aes(j)) +
  geom_density(alpha = 0.2, fill = branded_colors1[3]) +
  theme_bw(base_size = 24) +
  labs(x = "Avg. methylation (beta)", y = "-log10(p-value)")
ggsave("BatchEffect_Normal Pancreas_DensityPlot_T-Test.png", path= "./output/")


ttest %>% 
  sample_n(size = 30000) %>% 
  ggplot(aes(d, dm)) +
  geom_point(fill = "grey85", alpha = 0.2) +
  theme_bw(base_size = 24) +
  labs(x = "Delta (beta)", y = "-log10(p-value")
ggsave("BatchEffect_Normal Pancreas_T-Test.png", path= "./output/")


### look at sample correlation heatmap ---------------------------------------

heatmap_anno <- anno_normal %>% 
  select(source, conversion) %>% 
  as.data.frame()
rownames(heatmap_anno) <- anno_normal$arrayId

heat <- cor(betas_normal) %>% 
  pheatmap::pheatmap(annotation_col = heatmap_anno,
                     labels_row = heatmap_anno$source, 
                     show_colnames = FALSE)
save_pheatmap_pdf(heat, "./output/BatchEffect_Normal Pancreas_Heatmap cor.pdf")

# run UMAP umap
set.seed(45098)
umap_settings <- umap.defaults
umap_settings$n_neighbors = 15
umap_settings$min_dist = 0.2

#umap <- umap(d = sample_cor)
umap <- umap(d = t(betas_normal), config = umap_settings, ret_model = TRUE)

anno_normal <- anno_normal %>% 
  mutate(umap_x = umap$layout[, 1], 
         umap_y = umap$layout[, 2])

anno_normal %>% 
  ggplot(aes(umap_x, umap_y, col = source)) +
  geom_point(size = 4)


### Investigation 2: PanNets, DiDomenico vs. Jakel vs. Chan vs. Yachida vs. UMCU --------------------

# select PanNETs 
anno_pnet <- anno %>% 
  filter(tumorType == "PanNET")


# plot sample numbers per study
anno_pnet %>% 
  group_by(source) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(source, n, fill = source)) +
  geom_col(width = 0.4) +
  labs(x = NULL, y = "Normal tissues (n)") +
  #scale_x_discrete(limits = c("Chan", "DiDomenico", "Jakel", "UMCU", "yachida"),
                   #breaks = c("Chan", "DiDomenico", "Jakel", "UMCU", "yachida"),
                   #labels = c("Chan", "DiDomenico", "Jakel", "UMCU", "Yachida")) +
  theme_bw(base_size = 28) +
  theme(axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 16),
        axis.text = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = "none")+
  scale_fill_manual(values = branded_colors1)
ggsave("BatchEffect_PNET_No.png", path= "./output/")
 

# list idats, read raw data
idat <- list.files(path = "./data/ALL IDATS", pattern = "_Grn.idat")
file <- list.files(path = "./data/ALL IDATS", pattern = "_Grn.idat", full.names = TRUE)

file <- file[match(paste0(anno_pnet$arrayId, "_Grn.idat"), idat)]

file_size <- file.info(file)$size

idats_epic <- file[file_size > 10000000]
idats_450k <- file[file_size < 10000000]

raw_epic <- read.metharray(basenames = idats_epic, force = TRUE)
raw_450k <- read.metharray(basenames = idats_450k)

raw_pnet <- read.metharray(basenames = file, force = TRUE)
rm(idat, file)



rm(idat_pnet, file_pnet)


# investigate QC measures 1: bisulfite conversion efficiency 

conv_pnet2 <- getControlBeta(raw_pnet2, controls = "BISULFITE CONVERSION II") %>% 
  as_tibble() %>% 
  group_by(arrayId) %>% 
  summarise(conversion = min(value[channel == "Red"]) / max(value[channel == "Green"]))

anno_pnet2 <- anno_pnet2 %>% 
  left_join(conv_pnet2)

anno_pnet <- rbind(anno_pnet1, anno_pnet2)

#plot
anno_pnet %>% 
  ggplot(aes(source, conversion, fill = source)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(size = 2, width = 0.1) +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 16),
        axis.text = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = "none")+
  labs(x = NULL, y = "Control Score") +
  scale_x_discrete(limits = c("Chan", "DiDomenico", "Jakel", "yachida"),
                   breaks = c("Chan", "DiDomenico", "Jakel", "yachida"),
                   labels = c("Chan", "DiDomenico", "Jakel", "Yachida")) +
  scale_fill_manual(values = branded_colors)
ggsave("BatchEffect_PNET_Controle Score.png", path= "./output/")

# collect beta values, add avg. beta per sample to annotation
betas_pnet <- betas[, anno_pnet$arrayId]

betas_pnet_avg <- tibble(
  Chan = apply(betas_pnet[, anno_pnet$source == "Chan"], 1, mean, na.rm = TRUE),
  DiDomenico = apply(betas_pnet[, anno_pnet$source == "DiDomenico"], 1, mean, na.rm = TRUE), 
  Jakel = apply(betas_pnet[, anno_pnet$source == "Jakel"], 1, mean, na.rm = TRUE),
  Yachida = apply(betas_pnet[, anno_pnet$source == "yachida"], 1, mean, na.rm = TRUE)
)

#violin plot
betas_pnet_avg %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(name, value, fill = name)) +
  geom_violin(alpha = 0.8) +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 16),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=1, size = 16),
        axis.text = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = "none")+
  labs(x = "Beta values", y = "Density") +
  scale_fill_manual(values = branded_colors)
ggsave("BatchEffect_PNET_Violin_Beta's avg.png", path= "./output/")


betas_pnet_avg %>% 
  sample_n(size = 30000) %>% 
  ggplot(aes(Chan, Jakel)) +
  geom_bin2d(bins = 100) +
  paletteer::scale_fill_paletteer_c("scico::tokyo") +
  theme_bw(base_size = 28) +
  theme(axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 16),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=1, size = 16),
        axis.text = element_text(face = "bold", size = 10),
        legend.position = "none")+
  geom_abline(slope = 1, intercept = 0, col = "grey", lty = 2) +
  theme(legend.position = "none")
ggsave("Beta's avg_Chan-Jakel_pnet.png", path= "./output/")


anno_pnet <- anno_pnet %>% 
  mutate(mean_beta = apply(betas_pnet, 2, mean, na.rm = TRUE))

anno_pnet %>% 
  ggplot(aes(source, mean_beta, fill = source)) +
  geom_boxplot(alpha = 0.8) +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 16),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=1, size = 16),
        axis.text = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = "none") +
  labs(x = "", y = "Mean beta values") +
  scale_fill_manual(values = branded_colors)
ggsave("BatchEffect_PNET_AVG betas.png", path= "./output/")


# density plot of beta value distribution
betas_pnet_avg %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(value, fill = name)) +
  geom_density(alpha = 0.8) +
  theme_bw(base_size = 24) +
  labs(x = "Beta", y = "Density") +
  scale_fill_manual(values = branded_colors) +
  theme(legend.position = "none",
        strip.text = element_text(size = 10)) +
  facet_grid(rows = vars(name)) 
ggsave("BatchEffect_PNET_Density plots.png", path= "./output/")


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
  geom_point(size = 3) +
  #geom_label(label = anno_pnet$sampleName) +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 16),
        axis.title = element_text(face = "bold", hjust = 0.5, vjust=1, size = 16),
        axis.text = element_text(face = "bold", size = 10),
        legend.title=element_text(face = "bold", size=10), 
        legend.text=element_text(size=9)) +
  labs(x = "UMAP 1", y = "UMAP 2", alpha = "Conversion", col = "Source") +
  scale_color_manual(labels = c("Chan", "DiDomenico", "Jakel", "Yachida"), values = branded_colors)
ggsave("BatchEffect_PNET_UMAP.png", path= "./output/")


# t-SNE:
tsne_pnet <- Rtsne(t(betas_pnet[sample(1:nrow(betas_pnet), size = 5000), ]), perplexity = 15)

anno_pnet <- anno_pnet %>% 
  mutate(tsne_x = tsne_pnet$Y[, 1], 
         tsne_y = tsne_pnet$Y[, 2])

anno_pnet %>% 
  ggplot(aes(tsne_x, tsne_y, col = source, alpha = conversion)) +
  geom_point(size = 3) +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 16),
        axis.title = element_text(face = "bold", hjust = 0.5, vjust=1, size = 16),
        axis.text = element_text(face = "bold", size = 10),
        legend.title=element_text(face = "bold", size=10), 
        legend.text=element_text(size=9)) +
  labs(x = "t-SNE 1", y = "t-SNE 2", alpha = "Conversion", col = "Source") +
  scale_color_manual(labels = c("Chan", "DiDomenico", "Jakel", "Yachida"), values = branded_colors)
ggsave("BatchEffect_PNET_tSNE.png", path= "./output/")

### look at sample correlation heatmap ---------------------------------------

heatmap_anno <- anno_pnet %>% 
  select(source, conversion) %>% 
  as.data.frame()
rownames(heatmap_anno) <- anno_pnet$arrayId

heat <- cor(betas_pnet) %>% 
  pheatmap::pheatmap(annotation_col = heatmap_anno,
                     labels_row = heatmap_anno$source, 
                     show_colnames = FALSE)
save_pheatmap_pdf(heat, "./output/BatchEffect_PNET_Heatmap cor.pdf")
