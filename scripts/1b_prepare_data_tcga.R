# libraries

library(tidyverse)
library(minfi)
library(doParallel)

source("./0_helpers.R")


# load annotation
anno_tcga <- read_csv(file = "./annotation/TCGA/annotation_tcga_primaries.csv")

# get model probes
top_variable_probes <- readRDS(file = "./output/pancreas_top_variable_probes_training_set.rds")
top_variable_probes <- top_variable_probes[1:5000]
  
# scan idat files, select those from training and methylation data
idat_dir <- "/media/Yoda/General/idats_tcga_cgeisenberger/"
idats <- detect_idats(dir = idat_dir)

# remove cases where RED file is missing (n = 67 out of 15170)
idats <- idats %>% 
  mutate(red_exists = file.exists(str_replace(string = path, 
                                              pattern = "_Grn.idat", 
                                              replacement = "_Red.idat")))

idats <- idats %>% 
  filter(red_exists) %>% 
  select(-red_exists)


# check for which samples idats are available
anno_tcga <- anno_tcga %>% 
  mutate(available = (basename %in% idats$array_id)) %>% 
  filter(available) %>% 
  group_by(tissue) %>% 
  mutate(n_total = sum(available))

# remove pancreatic cancer (PAAD) and brain cancers (GBM, LGG)
anno_tcga <- anno_tcga %>% 
  filter(!tissue %in% c("PAAD", "GBM", "LGG"))

anno_tcga %>% pull(tissue) %>% table

# select corresponding idats
idats <- idats[match(anno_tcga$basename, idats$array_id), ]

# load training data (takes 17.4 minutes)
betas <- mclapply(X = idats$path,
                          FUN = function(p){
                            data <- minfi::read.metharray(p) %>% 
                              preprocessNoob %>% 
                              getBeta()
                            sub <- data[top_variable_probes, ]
                            return(data)
                          }, 
                  mc.cores = 96)

betas <- Reduce(f = cbind, x = betas)
betas <- as.matrix(betas)
colnames(betas) <- anno_tcga$basename
rownames(betas) <- top_variable_probes

# save data to files
saveRDS(object = anno_tcga[, 1:9], file = "./output/sample_anno_tcga.rds")
saveRDS(object = betas, file = "./input/betas_tcga_modelprobes.rds")



