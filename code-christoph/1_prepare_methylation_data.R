# Christoph Geisenberger
# github: @cgeisenberger
# last edited 21/11/2022



### load required packages -----------------------------------------------------

library(tidyverse)
library(minfi)



### prepare beta value matrices ------------------------------------------------

idats <- list.files(path = "./input/idat/",
                    recursive = TRUE,
                    full.names = TRUE, 
                    pattern = "_Grn.idat")


# infer platform from filesize
file_size <- file.info(idats)$size

idats_epic <- idats[file_size > 10000000]
idats_450k <- idats[file_size < 10000000]

raw_epic <- read.metharray(basenames = idats_epic, force = TRUE)
raw_450k <- read.metharray(basenames = idats_450k)


# preprocess and extract beta values

preprocessed_epic <- raw_epic %>%
  preprocessNoob(dyeMethod= "single")
rm(raw_epic)
saveRDS(object = preprocessed_epic, file = "./input/preprocessed_epic.rds")
rm(preprocessed_epic)

preprocessed_450k <- raw_450k %>%
  preprocessNoob(dyeMethod= "single")
rm(raw_450k)
saveRDS(object = preprocessed_450k, file = "./input/preprocessed_450k.rds")
rm(preprocessed_450k)

rm(idats, idats_450k, idats_epic, file_size)


# load preprocessed data

preprocessed_450k <- readRDS(file = "./input/preprocessed_450k.rds")
preprocessed_epic <- readRDS(file = "./input/preprocessed_epic.rds")


# remove probes on sex chromosomes

platform_anno_450k <- getAnnotation(preprocessed_450k)
filtered_450K <- preprocessed_450k[!platform_anno_450k$chr %in% c("chrX","chrY"), ]
rm(platform_anno_450k)

platform_anno_epic <- getAnnotation(preprocessed_epic)
filtered_epic <- preprocessed_epic[!platform_anno_epic$chr %in% c("chrX","chrY"), ]
rm(platform_anno_epic)


# remove probes containing SNPs and multimappers

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

rm(filtered_450K, filtered_epic)


# filter out X reactive probes

x_reactive <- read_csv("./input/48639-non-specific-probes-Illumina450k.csv", 
                       col_names = "probe_id")

betas_filtered <- betas_filtered[!rownames(betas_filtered) %in% x_reactive$probe_id, ]
rm(x_reactive)



# create unfiltered beta value matrix

betas_450k <- getBeta(preprocessed_450k)
betas_epic <- getBeta(preprocessed_epic)

common_probes <- intersect(rownames(betas_450k), rownames(betas_epic))
betas <- cbind(betas_450k[common_probes, ], betas_epic[common_probes, ])

rm(betas_450k, betas_epic)
rm(preprocessed_450k, preprocessed_epic)
rm(common_probes)


# save beta values data to disk 
saveRDS(object = betas, file = "./input/pancreas_betas_everything.rds")
saveRDS(object = betas_filtered, file = "./input/pancreas_betas_filtered.rds")
