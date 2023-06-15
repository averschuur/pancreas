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

# plot control score
anno_normal %>% 
  ggplot(aes(source, conversion, fill = source)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(size = 2, width = 0.1) +
  theme_bw(base_size = 24) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Control Score") +
  scale_fill_manual(values = branded_colors1)

#stastiek
stat <- anno_normal[,c(2,7)]
attach(stat)
t.test(conversion~source)



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

#
betas_normal_avg %>% 
  sample_n(size = 30000) %>% 
  ggplot(aes(DiDomenico, Jakel)) +
  geom_bin2d(bins = 100) +
  paletteer::scale_fill_paletteer_c("scico::tokyo") +
  theme_bw(base_size = 28) +
  geom_abline(slope = 1, intercept = 0, col = "grey", lty = 2) +
  theme(legend.position = "none")
ggsave("BatchEffect_Normal Pancreas_Beta's avg.png", path= "./output/")

#stastiek
Jakel <- betas_normal_avg[,2]
jakel = "Jakel"
Jakel <- cbind(Jakel, jakel)
colnames(Jakel) <- c("avg_beta", "source")
DiDominico <- betas_normal_avg[,1]
didomenico = "DiDomenico"
DiDominico <- cbind(DiDominico, didomenico)
colnames(DiDominico) <- c("avg_beta", "source")
stat <- rbind(Jakel, DiDominico)
attach(stat)
t.test(avg_beta~source)


#mean beta values
anno_normal <- anno_normal %>% 
  mutate(mean_beta = apply(betas_normal, 2, mean, na.rm = TRUE))

anno_normal %>% 
  ggplot(aes(source, mean_beta, fill = source)) +
  geom_boxplot(alpha = 0.8) +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 16),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=1, size = 16),
        axis.text = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = "none") +
  labs(x = "", y = "Mean beta values") +
  scale_fill_manual(values = branded_colors1)

#stastiek
stat <- anno_normal[,c(2,8)]
attach(stat)
t.test(mean_beta~source)

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


### Investigation 2: pNET, DiDomenico vs. Jakel vs. Chan vs. Yachida vs. UMCU --------------------

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

conv_pnet <- getControlBeta(raw_epic, controls = "BISULFITE CONVERSION II") %>% 
  as_tibble() %>% 
  group_by(arrayId) %>% 
  summarise(conversion = min(value[channel == "Red"]) / max(value[channel == "Green"]))

conv_pnet_2 <- getControlBeta(raw_450k, controls = "BISULFITE CONVERSION II") %>% 
  as_tibble() %>% 
  group_by(arrayId) %>% 
  summarise(conversion = min(value[channel == "Red"]) / max(value[channel == "Green"]))

conv_pnet <- bind_rows(conv_pnet, conv_pnet_2)

anno_pnet <- anno_pnet %>% 
  left_join(conv_pnet)


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
  #scale_x_discrete(limits = c("Chan", "DiDomenico", "Jakel", "yachida"),
   #                breaks = c("Chan", "DiDomenico", "Jakel", "yachida"),
    #               labels = c("Chan", "DiDomenico", "Jakel", "Yachida")) +
  scale_fill_manual(values = branded_colors1)

#stastiek
stat <- anno_pnet[,c(2,7)]
attach(stat)
boxplot(conversion~source)

pg<-lm(conversion~source)
anova(pg)

mfit2 <- aov(conversion~factor(source))
TukeyHSD(mfit2)

#                          diff          lwr         upr     p adj
#DiDomenico-Chan    -0.19300185 -0.261566205 -0.12443750 0.0000000
#Jakel-Chan         -0.11532088 -0.201194208 -0.02944755 0.0027709
#UMCU-Chan          -0.23692291 -0.338620279 -0.13522554 0.0000000
#Jakel-DiDomenico    0.07768097  0.001475953  0.15388599 0.0434010
#yachida-DiDomenico  0.14091682  0.077910524  0.20392311 0.0000001
#UMCU-Jakel         -0.12160203 -0.228599642 -0.01460441 0.0173903
#yachida-UMCU        0.18483787  0.086801769  0.28287398 0.0000074

# collect beta values, add avg. beta per sample to annotation
betas_pnet <- betas[, anno_pnet$arrayId]

betas_pnet_avg <- tibble(
  Chan = apply(betas_pnet[, anno_pnet$source == "Chan"], 1, mean, na.rm = TRUE),
  DiDomenico = apply(betas_pnet[, anno_pnet$source == "DiDomenico"], 1, mean, na.rm = TRUE), 
  Jakel = apply(betas_pnet[, anno_pnet$source == "Jakel"], 1, mean, na.rm = TRUE),
  UMCU = apply(betas_pnet[, anno_pnet$source == "UMCU"], 1, mean, na.rm = TRUE),
  Yachida = apply(betas_pnet[, anno_pnet$source == "yachida"], 1, mean, na.rm = TRUE)
)

# density plot of beta value distribution
betas_pnet_avg %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(value, fill = name)) +
  geom_density(alpha = 0.8) +
  theme_bw(base_size = 24) +
  labs(x = "Beta", y = "Density") +
  scale_fill_manual(values = branded_colors1) +
  theme(legend.position = "none",
        strip.text = element_text(size = 10)) +
  facet_grid(rows = vars(name)) 


#stastiek
Chan <- betas_pnet_avg[,1]
chan = "Chan"
Chan <- cbind(Chan, chan)
colnames(Chan) <- c("avg_beta", "source")
DiDominico <- betas_pnet_avg[,2]
didomenico = "DiDomenico"
DiDominico <- cbind(DiDominico, didomenico)
colnames(DiDominico) <- c("avg_beta", "source")
Jakel <- betas_pnet_avg[,3]
jakel = "Jakel"
Jakel <- cbind(Jakel, jakel)
colnames(Jakel) <- c("avg_beta", "source")
UMCU <- betas_pnet_avg[,4]
umcu = "UMCU"
UMCU <- cbind(UMCU, umcu)
colnames(UMCU) <- c("avg_beta", "source")
Yachida <- betas_pnet_avg[,5]
yachida = "Yachida"
Yachida <- cbind(Yachida, yachida)
colnames(Yachida) <- c("avg_beta", "source")
stat <- rbind(Chan, DiDominico)
stat <- rbind(stat, Jakel)
stat <- rbind(stat, UMCU)
stat <- rbind(stat, Yachida)

stat<-stat[complete.cases(stat),]

attach(stat)
boxplot(avg_beta~source)

pg<-lm(avg_beta~source)
anova(pg)

mfit2 <- aov(avg_beta~factor(source))
TukeyHSD(mfit2)

# obtain mean betas
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
  scale_fill_manual(values = branded_colors1)

#stastiek
stat <- anno_pnet[,c(2,8)]
attach(stat)
boxplot(mean_beta~source)

pg<-lm(mean_beta~source)
anova(pg)

mfit2 <- aov(mean_beta~factor(source))
TukeyHSD(mfit2)

#                           diff           lwr         upr     p adj
#yachida-DiDomenico  0.024114304  0.0004842730 0.047744335 0.0430417



# UMAP
umap_pnet <- umap::umap(d = t(betas_pnet[sample(1:nrow(betas_pnet), size = 5000), ]))


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





### Investigation 3: ACC tissue, Benhamida vs UMCU vs vs. Jakel ------------

# keep only ACC tissue
anno_acc <- anno %>% 
  filter(tumorType == "ACC")

anno_acc %>% 
  group_by(source) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(source, n, fill = source)) +
  geom_col(width = 0.4) +
  labs(x = NULL, y = "Normal tissues (n)") +
  theme_bw(base_size = 28) +
  theme(legend.position = "none") +
  scale_fill_manual(values = branded_colors1) +
  ylim(c(0, 30))

## obtain conversionscore
#list idats, read raw data
idat <- list.files(path = "./data/ALL IDATS", pattern = "_Grn.idat")
file <- list.files(path = "./data/ALL IDATS", pattern = "_Grn.idat", full.names = TRUE)

file <- file[match(paste0(anno_acc$arrayId, "_Grn.idat"), idat)]

file_size <- file.info(file)$size

idats_epic <- file[file_size > 10000000]
idats_450k <- file[file_size < 10000000]

raw_epic <- read.metharray(basenames = idats_epic, force = TRUE)
raw_450k <- read.metharray(basenames = idats_450k)

rm(idat, file, idats_epic, idats_450k)


# investigate QC measures 1: bisulfite conversion efficiency 

conv_acc <- getControlBeta(raw_epic, controls = "BISULFITE CONVERSION II") %>% 
  as_tibble() %>% 
  group_by(arrayId) %>% 
  summarise(conversion = min(value[channel == "Red"]) / max(value[channel == "Green"]))

conv_acc_2 <- getControlBeta(raw_450k, controls = "BISULFITE CONVERSION II") %>% 
  as_tibble() %>% 
  group_by(arrayId) %>% 
  summarise(conversion = min(value[channel == "Red"]) / max(value[channel == "Green"]))

conv_acc <- bind_rows(conv_acc, conv_acc_2)

anno_acc <- anno_acc %>% 
  left_join(conv_acc)

#stastiek
stat <- anno_acc[,c(2,7)]
attach(stat)
boxplot(conversion~source)

pg<-lm(conversion~source)
anova(pg)

mfit2 <- aov(conversion~factor(source))
TukeyHSD(mfit2)

#                      diff        lwr          upr     p adj
#Jakel-Benhamida -0.1027644 -0.1966770 -0.008851712 0.0290909
#UMCU-Benhamida  -0.2160587 -0.3588708 -0.073246600 0.0018377

#plot conversion score
anno_acc %>% 
  ggplot(aes(source, conversion, fill = source)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(size = 2, width = 0.1) +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 16),
        axis.text = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = "none")+
  labs(x = NULL, y = "Control Score") +
  scale_fill_manual(values = branded_colors1)

# collect beta values, add avg. beta per sample to annotation
betas_acc <- betas[, anno_acc$arrayId]

betas_acc_avg <- tibble(
  Benhamida = apply(betas_pnet[, anno_acc$source == "Benhamida"], 1, mean, na.rm = TRUE), 
  Jakel = apply(betas_pnet[, anno_acc$source == "Jakel"], 1, mean, na.rm = TRUE),
  UMCU = apply(betas_pnet[, anno_acc$source == "UMCU"], 1, mean, na.rm = TRUE)
)

# density plot of beta value distribution
betas_acc_avg %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(value, fill = name)) +
  geom_density(alpha = 0.8) +
  theme_bw(base_size = 24) +
  labs(x = "Beta", y = "Density") +
  scale_fill_manual(values = branded_colors1) +
  theme(legend.position = "none",
        strip.text = element_text(size = 10)) +
  facet_grid(rows = vars(name)) 

#stastiek
Benhamida <- betas_acc_avg[,1]
benhamida = "Benhamida"
Benhamida <- cbind(Benhamida, benhamida)
colnames(Benhamida) <- c("avg_beta", "source")
Jakel <- betas_acc_avg[,2]
jakel = "Jakel"
Jakel <- cbind(Jakel, jakel)
colnames(Jakel) <- c("avg_beta", "source")
UMCU <- betas_acc_avg[,3]
umcu = "UMCU"
UMCU <- cbind(UMCU, umcu)
colnames(UMCU) <- c("avg_beta", "source")
stat <- rbind(Benhamida, Jakel)
stat <- rbind(stat, UMCU)

attach(stat)
boxplot(avg_beta~source)

pg<-lm(avg_beta~source)
anova(pg)

mfit2 <- aov(avg_beta~factor(source))
TukeyHSD(mfit2)

# obtain mean beta values
anno_acc <- anno_acc %>% 
  mutate(mean_beta = apply(betas_acc, 2, mean, na.rm = TRUE))

anno_acc %>% 
  ggplot(aes(source, mean_beta, fill = source)) +
  geom_boxplot(alpha = 0.8) +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 16),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=1, size = 16),
        axis.text = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = "none") +
  labs(x = "", y = "Mean beta values") +
  scale_fill_manual(values = branded_colors1)

#stastiek
stat <- anno_acc[,c(2,8)]
attach(stat)
boxplot(mean_beta~source)

pg<-lm(mean_beta~source)
anova(pg)

mfit2 <- aov(mean_beta~factor(source))
TukeyHSD(mfit2)

# UMAP
umap_acc <- umap::umap(d = t(betas_acc[sample(1:nrow(betas_acc), size = 5000), ]))


anno_acc <- anno_acc %>% 
  mutate(umap_x = umap_acc$layout[, 1], 
         umap_y = umap_acc$layout[, 2])

anno_acc %>% 
  ggplot(aes(umap_x, umap_y, col = source, alpha = conversion)) +
  geom_point(size = 3) +
  #geom_label(label = anno_pnet$sampleName) +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 16),
        axis.title = element_text(face = "bold", hjust = 0.5, vjust=1, size = 16),
        axis.text = element_text(face = "bold", size = 10),
        legend.title=element_text(face = "bold", size=10), 
        legend.text=element_text(size=9)) +
  labs(x = "UMAP 1", y = "UMAP 2", alpha = "Conversion", col = "Source") 


### Investigation 4: SPN tissue, Selenica vs UMCU ------------

# keep only spn tissue
anno_spn <- anno %>% 
  filter(tumorType == "SPN")

anno_spn %>% 
  group_by(source) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(source, n, fill = source)) +
  geom_col(width = 0.4) +
  labs(x = NULL, y = "Normal tissues (n)") +
  theme_bw(base_size = 28) +
  theme(legend.position = "none") +
  scale_fill_manual(values = branded_colors1) +
  ylim(c(0, 30))

## obtain conversionscore
#list idats, read raw data
idat <- list.files(path = "./data/ALL IDATS", pattern = "_Grn.idat")
file <- list.files(path = "./data/ALL IDATS", pattern = "_Grn.idat", full.names = TRUE)

file <- file[match(paste0(anno_spn$arrayId, "_Grn.idat"), idat)]

file_size <- file.info(file)$size

idats_epic <- file[file_size > 10000000]
idats_450k <- file[file_size < 10000000]

raw_epic <- read.metharray(basenames = idats_epic, force = TRUE)

rm(idat, file, idats_epic, idats_450k)


# investigate QC measures 1: bisulfite conversion efficiency 

conv_spn <- getControlBeta(raw_epic, controls = "BISULFITE CONVERSION II") %>% 
  as_tibble() %>% 
  group_by(arrayId) %>% 
  summarise(conversion = min(value[channel == "Red"]) / max(value[channel == "Green"]))


anno_spn <- anno_spn %>% 
  left_join(conv_spn)


#plot conversion score
anno_spn %>% 
  ggplot(aes(source, conversion, fill = source)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(size = 2, width = 0.1) +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 16),
        axis.text = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = "none")+
  labs(x = NULL, y = "Control Score") +
  scale_fill_manual(values = branded_colors1)

#stastiek
stat <- anno_spn[,c(2,7)]
attach(stat)
t.test(conversion~source)

# collect beta values, add avg. beta per sample to annotation
betas_spn <- betas[, anno_spn$arrayId]

betas_spn_avg <- tibble(
  Selenica = apply(betas_spn[, anno_spn$source == "Selenica"], 1, mean, na.rm = TRUE),
  UMCU = apply(betas_spn[, anno_spn$source == "UMCU"], 1, mean, na.rm = TRUE)
)

# density plot of beta value distribution
betas_spn_avg %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(value, fill = name)) +
  geom_density(alpha = 0.8) +
  theme_bw(base_size = 24) +
  labs(x = "Beta", y = "Density") +
  scale_fill_manual(values = branded_colors1) +
  theme(legend.position = "none",
        strip.text = element_text(size = 10)) +
  facet_grid(rows = vars(name)) 

#stastiek
Selenica <- betas_spn_avg[,1]
selenica = "Selenica"
Selenica <- cbind(Selenica, jakel)
colnames(Selenica) <- c("avg_beta", "source")
UMCU <- betas_normal_avg[,2]
umcu = "UMCU"
UMCU <- cbind(UMCU, umcu)
colnames(UMCU) <- c("avg_beta", "source")
stat <- rbind(Selenica, UMCU)
attach(stat)
t.test(avg_beta~source)


anno_spn <- anno_spn %>% 
  mutate(mean_beta = apply(betas_spn, 2, mean, na.rm = TRUE))

anno_spn %>% 
  ggplot(aes(source, mean_beta, fill = source)) +
  geom_boxplot(alpha = 0.8) +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 16),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=1, size = 16),
        axis.text = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = "none") +
  labs(x = "", y = "Mean beta values") +
  scale_fill_manual(values = branded_colors1)

#stastiek
stat <- anno_spn[,c(2,8)]
attach(stat)
t.test(mean_beta~source)

# UMAP
umap_spn <- umap::umap(d = t(betas_spn[sample(1:nrow(betas_spn), size = 5000), ]))


anno_spn <- anno_spn %>% 
  mutate(umap_x = umap_spn$layout[, 1], 
         umap_y = umap_spn$layout[, 2])

anno_spn %>% 
  ggplot(aes(umap_x, umap_y, col = source, alpha = conversion)) +
  geom_point(size = 3) +
  #geom_label(label = anno_pnet$sampleName) +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 16),
        axis.title = element_text(face = "bold", hjust = 0.5, vjust=1, size = 16),
        axis.text = element_text(face = "bold", size = 10),
        legend.title=element_text(face = "bold", size=10), 
        legend.text=element_text(size=9)) +
  labs(x = "UMAP 1", y = "UMAP 2", alpha = "Conversion", col = "Source") 