# Christoph Geisenberger
# github: @cgeisenberger
# last edited 25/01/2023 by GEI



### load required packages -----------------------------------------------------

library(tidyverse)
library(minfi)
library(doParallel)



# load annotation --------------------------------------------------------------

anno_files <- list.files(path = "./annotation/",
                         pattern = ".csv",
                         full.names = TRUE)
anno_files <- anno_files[!grepl(x = anno_files, pattern = "*tcga*")]
anno_files <- anno_files[!grepl(x = anno_files, pattern = "*annotation_umcu_paired_samples*")]
anno_files <- anno_files[!grepl(x = anno_files, pattern = "*sample_annotation_conversion_scores*")]


anno <- lapply(as.list(anno_files), read_csv)
anno <- Reduce(f = bind_rows, x = anno)
anno



# prepare methylation data -----------------------------------------------------

idats <- list.files(path = "./idat-pancreas/",
                    recursive = TRUE,
                    full.names = TRUE, 
                    pattern = "_Grn.idat")


# infer platform from filesize
file_size <- file.info(idats)$size

idats_epic <- idats[file_size > 10000000]
idats_450k <- idats[file_size < 10000000]

raw_epic <- read.metharray(basenames = idats_epic, force = TRUE)
raw_450k <- read.metharray(basenames = idats_450k)

# preprocess data
preprocessed_epic <- raw_epic %>%
  preprocessNoob(dyeMethod= "single")
preprocessed_450k <- raw_450k %>%
  preprocessNoob(dyeMethod= "single")

# combine 450k and EPIC data, save full (combined) data to disk
betas_epic <- getBeta(preprocessed_epic)
betas_450k <- getBeta(preprocessed_450k)

common_probes <- intersect(rownames(betas_epic), rownames(betas_450k))
betas_combined <- cbind(betas_450k[common_probes, ], betas_epic[common_probes, ])

# clean up
rm(raw_epic, raw_450k)
rm(idats, idats_450k, idats_epic, file_size)

# identify probes (1) on X/Y chr, (2) SNPs/multimappers or (3) are crossreactive

# sex chromosomes
platform_anno <- getAnnotation(preprocessed_450k)
filtered_probes <- list()
filtered_probes$xy <- which(platform_anno$chr %in% c("chrX","chrY"), arr.ind = TRUE)
filtered_probes$xy <- rownames(preprocessed_450k)[filtered_probes$xy]

# SNPs and multimappers
betas_450k_filtered <- preprocessed_450k %>% 
  mapToGenome %>% 
  dropLociWithSnps %>% 
  getBeta()
filtered_probes$snp <- intersect(rownames(betas_450k), rownames(betas_450k_filtered))
filtered_probes$snp <- which(rownames(betas_450k) %in% filtered_probes$snp, arr.ind = TRUE)
filtered_probes$snp <- rownames(betas_450k)[-filtered_probes$snp]

# cross reactive
filtered_probes$xr <- read_csv("./input/48639-non-specific-probes-Illumina450k.csv", 
                       col_names = "probe_id")
filtered_probes$xr <- filtered_probes$xr$probe_id[-1] %>% as.character()

# clean up and save to file
saveRDS(object = filtered_probes, file = "./input/filtered_probes_list.rds")
rm(platform_anno, betas_450k_filtered)

# filter data 
probes_remove <- Reduce(f = union, x = filtered_probes)
probes_remove <- intersect(probes_remove, rownames(betas_combined))
probes_remove <- match(probes_remove, rownames(betas_combined))

betas_combined_filtered <- betas_combined[-probes_remove, ]

# clean up
rm(betas_450k, betas_epic, filtered_probes, probes_remove, common_probes)
rm(preprocessed_450k, preprocessed_epic)



# double check annotation vs. beta values, save files to disk ------------------

dim(betas_combined)

all(colnames(betas_combined) == colnames(betas_combined_filtered))
all(colnames(betas_combined) %in% anno$arrayId)
all(anno$arrayId %in% colnames(betas_combined))


# anno
saveRDS(object = anno, file = "./input/sample_annotation.rds")

# all betas
saveRDS(object = betas_combined, file = "./betas_pancreas_everything.rds")

# filtered betas
saveRDS(object = betas_combined_filtered, 
        file = "./betas_pancreas_filtered.rds")
