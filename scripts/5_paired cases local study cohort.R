# Christoph Geisenberger
# github: @cgeisenberger
# last edited 04/01/2023 by AV Verschuur


### Load required packages and sources ----------------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(minfi)

source("./scripts/0_helpers.R")


# load annotation, start from raw data

anno_umcu <- read_csv(file = "./annotation/annotation_umcu.csv")
anno_paired <- read_csv(file = "./annotation/annotation_umcu_paired_samples.csv")

idats <- detect_idats(dir = "./data/UMCU/Idat Files_All/")
idats <- idats[match(anno_umcu$arrayId, idats$array_id), ]
#idats <- idats[1:22,]

# infer platform from filesize, apparantly different generations of EPIC were used
file_size <- file.info(idats$path)$size

# load raw data
raw <- minfi::read.metharray(idats$path, force = TRUE)

# extract SNP values to make sure samples are matched
snp_beta <- minfi::getSnpBeta(raw)

# extract beta values
betas <- minfi::preprocessNoob(raw, dyeMethod= "single") %>% 
  getBeta()
rm(raw)

### Analysis of all matched cases together ---------------------------------------
i = 9


test <- betas[, c(anno_paired$array_id_1[i], anno_paired$array_id_2[i])]
colnames(test) <- c("sample1", "sample2")
dim(test)

test %>% 
  as_tibble() %>% 
  slice_sample(n = 50000) %>% 
  ggplot(aes(sample1, sample2)) +
  #geom_bin2d(bins = 50) +
  geom_point(pch = ".", col = branded_colors2[2], alpha = 0.4) +
  scale_fill_distiller(palette = "PuRd", direction = 1) +
  theme_classic(base_size = 20)

test <- cor(betas[, anno_paired$array_id_1], betas[, anno_paired$array_id_2])
#test[lower.tri(test)] <- NA
test0 <- test %>% 
  as_tibble(rownames = "source") %>% 
  pivot_longer(cols = 2:10, names_to = "target", values_to = "cor")

test2 <- anno_paired %>% 
  select(array_id_1, array_id_2) %>% 
  rename(array_id_1 = "source", array_id_2 = "matched")

test <- test %>% 
  left_join(test2)

test <- test %>% 
  mutate(matched = ifelse(matched == target, 1, 0))

test %>% 
  ggplot(aes(source, cor, col = as.factor(matched))) +
  geom_jitter(width = 0.1) +
  theme_classic(base_size = 20)
ggsave("matched cases.pdf", path= "./plots/", dpi=500)

test %>% 
  ggplot(aes(as.factor(matched), cor)) +
  geom_boxplot() +
  theme_classic(base_size = 20)

### Matched cases biopsy -------------------------------------------------------
anno_paired_b <- anno_paired %>%
  filter(location_2 == 'biopsy')

i = 4


test <- betas[, c(anno_paired_b$array_id_2[i], anno_paired_b$array_id_1[i])]
colnames(test) <- c("sample1", "sample2")
dim(test)

test %>% 
  as_tibble() %>% 
  slice_sample(n = 50000) %>% 
  ggplot(aes(sample1, sample2)) +
  #geom_bin2d(bins = 50) +
  geom_point(pch = ".", col = branded_colors2[2], alpha = 0.4) +
  scale_fill_distiller(palette = "PuRd", direction = 1) +
  theme_classic(base_size = 20)

test <- cor(betas[, anno_paired_b$array_id_2], betas[, anno_paired_b$array_id_1])
#test[lower.tri(test)] <- NA
test <- test %>% 
  as_tibble(rownames = "source") %>% 
  pivot_longer(cols = 2:5, names_to = "target", values_to = "cor")

label <- anno_paired_b[,c(1,3)]
colnames(label) <- c("target", "label")
label2 <- anno_paired_b[,c(2,4)]
colnames(label2) <- c("source", "label_s")


test2 <- anno_paired_b %>% 
  select(array_id_2, array_id_1) %>% 
  rename(array_id_2 = "source", array_id_1 = "matched")

test <- test %>% 
  left_join(test2)

test <- test %>% 
  mutate(matched = ifelse(matched == target, 1, 0))

test <- full_join(test, label, by = "target")
test <- full_join(test, label2, by = "source")

test %>% 
  ggplot(aes(source, cor, col = as.factor(matched))) +
 # geom_jitter(width = 0.1) +
  geom_jitter(width = 0.1, height = 0.25, size = 1) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 30) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank())

test %>% 
  ggplot(aes(as.factor(matched), cor, fill = matched)) +
  geom_boxplot() +
  theme_bw(base_size = 18) 

# statistics
my_data <- test[,c(4,3)]

group_by(my_data, matched) %>%
  summarise(
    count = n(),
    mean = mean(cor, na.rm = TRUE),
    sd = sd(cor, na.rm = TRUE)
  )

res <- t.test(cor ~ as.factor(matched), data = my_data, paired = FALSE)


### Matched cases metastases----------------------------------------------------
anno_paired_m <- anno_paired %>%
  filter(location_2 == 'metastasis')

i = 4


test <- betas[, c(anno_paired_m$array_id_2[i], anno_paired_m$array_id_1[i])]
colnames(test) <- c("sample1", "sample2")
dim(test)

test %>% 
  as_tibble() %>% 
  slice_sample(n = 50000) %>% 
  ggplot(aes(sample1, sample2)) +
  #geom_bin2d(bins = 50) +
  geom_point(pch = ".", col = branded_colors2[2], alpha = 0.4) +
  scale_fill_distiller(palette = "PuRd", direction = 1) +
  theme_classic(base_size = 20)

test <- cor(betas[, anno_paired_m$array_id_2], betas[, anno_paired_m$array_id_1])
#test[lower.tri(test)] <- NA
test <- test %>% 
  as_tibble(rownames = "source") %>% 
  pivot_longer(cols = 2:6, names_to = "target", values_to = "cor")

test2 <- anno_paired_m %>% 
  select(array_id_2, array_id_1) %>% 
  rename(array_id_2 = "source", array_id_1 = "matched")

test <- test %>% 
  left_join(test2)

test <- test %>% 
  mutate(matched = ifelse(matched == target, 1, 0))

label <- anno_paired_m[,c(1,3)]
colnames(label) <- c("target", "label")
label2 <- anno_paired_m[,c(2,4)]
colnames(label2) <- c("source", "label_s")

test[test == "V4"] <- "201533570026_R01C01"
test <- full_join(test, label, by = "target")
test <- full_join(test, label2, by = "source")





test <- test[c(1,2,3,7,8,9,10,14,15,16,17,21,22,23,24,28, 29, 30, 31, 35),]

test %>% 
  ggplot(aes(source, cor, col = as.factor(matched))) +
  # geom_jitter(width = 0.1) +
  geom_jitter(width = 0.1, height = 0.25, size = 1) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 30) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank())

test %>% 
  ggplot(aes(as.factor(matched), cor, fill = matched)) +
  geom_boxplot() +
  theme_bw(base_size = 18) 

# statistics
my_data <- test[,c(4,3)]

group_by(my_data, matched) %>%
  summarise(
    count = n(),
    mean = mean(cor, na.rm = TRUE),
    sd = sd(cor, na.rm = TRUE)
  )

t.test(cor ~ as.factor(matched), data = my_data, paired = FALSE)

