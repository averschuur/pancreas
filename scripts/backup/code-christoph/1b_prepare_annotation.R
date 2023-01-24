# Christoph Geisenberger
# github: @cgeisenberger
# last edited 21/11/2022



### load required packages -----------------------------------------------------

library(tidyverse)
library(minfi)


### load annotation ------------------------------------------------------------

anno_files <- list.files(path = "./annotation/",
                         pattern = ".csv",
                         full.names = TRUE)
anno_files <- anno_files[!grepl(x = anno_files, pattern = "*conversion_scores*")]
anno_files <- anno_files[!grepl(x = anno_files, pattern = "*annotation_umcu_paired_samples*")]


anno <- lapply(as.list(anno_files), read_csv)
anno <- Reduce(f = bind_rows, x = anno)
anno

saveRDS(object = anno, file = "./input/sample_annotation.rds")


# double check annotation and beta values correspond

betas <- readRDS(file = "./input/pancreas_betas_everything.rds")
betas_filtered <- readRDS(file = "./input/pancreas_betas_filtered.rds")

all(colnames(betas) %in% anno$arrayId)
all(anno$arrayId %in% colnames(betas))

all(colnames(betas_filtered) %in% anno$arrayId)
all(anno$arrayId %in% colnames(betas_filtered))

all(colnames(betas) == colnames(betas_filtered))
