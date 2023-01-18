# Christoph Geisenberger
# github: @cgeisenberger
# last edited 04/01/2023 by AV Verschuur


### Load required packages and sources ----------------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(minfi)

source("./scripts/0_functions.R")
source("./scripts/0_branded_colors.R")


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
test <- test %>% 
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

test %>% 
  ggplot(aes(as.factor(matched), cor)) +
  geom_boxplot() +
  theme_classic(base_size = 20)