# Christoph Geisenberger
# github: @cgeisenberger
# last edited 28/08/2022



### load required packages -----------------------------------------------------

library(tidyverse)
library(minfi)



### load annotation ------------------------------------------------------------

anno_files <- list.files(path = "./input/annotation/cleaned/",
                        pattern = ".csv",
                        full.names = TRUE)

anno <- lapply(as.list(anno_files), read_csv)

anno <- Reduce(f = bind_rows, x = anno)

saveRDS(object = anno, file = "./input/sample_annotation.rds")



### prepare methylation data ---------------------------------------------------

idats <- list.files(path = "./input/idat/",
                    full.names = TRUE, 
                    pattern = "_Grn.idat")


# infer platform from filesize
file_size <- file.info(idats)$size

idats_epic <- idats[file_size > 10000000]
idats_450k <- idats[file_size < 10000000]

raw_epic <- read.metharray(basenames = idats_epic, force = TRUE)
raw_450k <- read.metharray(basenames = idats_450k)



### perform preprocessing and extract beta values ------------------------------

preprocessed_epic <- raw_epic %>%
  preprocessNoob(dyeMethod= "single")

preprocessed_450k <- raw_450k %>%
  preprocessNoob(dyeMethod= "single")



### Perform filtering: Remove chrX/chrY, probes with SNPs and multimappers -----


# remove X and Y chromosomes 

platform_anno_450k <- getAnnotation(preprocessed_450k)
filtered_450K <- preprocessed_450k[!platform_anno_450k$chr %in% c("chrX","chrY"), ]
rm(platform_anno_450k)

platform_anno_epic <- getAnnotation(preprocessed_epic)
filtered_epic <- preprocessed_epic[!platform_anno_epic$chr %in% c("chrX","chrY"), ]
rm(platform_anno_epic)



# filter out SNPs

filtered_450K <- filtered_450K %>% 
  mapToGenome %>% 
  dropLociWithSnps
filtered_epic <- filtered_epic %>% 
  mapToGenome %>% 
  dropLociWithSnps


# extract Beta values, combine and clean up

filtered_450K <- getBeta(filtered_450K)
filtered_epic <- getBeta(filtered_epic)

common_probes <- intersect(rownames(filtered_450K), rownames(filtered_epic))
betas_filtered <- cbind(filtered_450K[common_probes, ], filtered_epic[common_probes, ])

rm(raw_450k, raw_epic)
rm(filtered_450K, filtered_epic)


# filter out X reactive probes

x_reactive <- read_csv("./input/48639-non-specific-probes-Illumina450k.csv", 
                       col_names = "probe_id")

betas_filtered <- betas_filtered[!rownames(betas_filtered) %in% x_reactive$probe_id, ]
rm(x_reactive)



### create unfiltered beta value matrix ----------------------------------------

betas_450k <- getBeta(preprocessed_450k)
betas_epic <- getBeta(preprocessed_epic)

common_probes <- intersect(rownames(betas_450k), rownames(betas_epic))
betas <- cbind(betas_450k[common_probes, ], betas_epic[common_probes, ])

rm(betas_450k, betas_epic)
rm(preprocessed_450k, preprocessed_epic)


### match sample annotation and beta value matrix ------------------------------

samples_common <- intersect(anno$array_id, colnames(betas))

length(samples_common) == nrow(anno)
length(samples_common) == ncol(betas)

anno <- anno %>% 
  filter(array_id %in% samples_common)

betas <- betas[, anno$array_id]
betas_filtered <- betas_filtered[, anno$array_id]

# save everything to disk
saveRDS(object = betas, file = "./input/betas.rds")
saveRDS(object = betas_filtered, file = "./input/betas_filtered.rds")
saveRDS(object = anno, file = "./input/sample_annotation_incomplete.rds")



