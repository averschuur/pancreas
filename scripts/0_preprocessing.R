
### Load required packages -----------------------------------------
library(tidyverse)
library(minfi)



### Prepare methylation data ---------------------------------------

idats <- list.files(path = "./input/ALL IDATS",
                    full.names = TRUE, 
                    pattern = "_Grn.idat")

# infer platform from filesize

fileSizes <- file.info(idats)$size

idatsEpic <- idats[fileSizes > 10000000]
idats450k <- idats[fileSizes < 10000000]

arrayRawEpic <- read.metharray(basenames = idatsEpic, force = TRUE)
arrayRaw450k <- read.metharray(basenames = idats450k)


### Normalize and filter probes ------------------------------------

filtered450K <- arrayRaw450k %>%
  preprocessNoob(dyeMethod= "single") %>%
  ratioConvert(what = "both", keepCN = TRUE) %>%
  mapToGenome

filteredEPIC <- arrayRawEpic %>%
  preprocessNoob(dyeMethod= "single") %>%
  ratioConvert(what = "both", keepCN = TRUE) %>%
  mapToGenome

# remove any probes that have failed in one or more samples

detP <- arrayRaw450k %>%
  detectionP
detP <- detP[match(featureNames(filtered450K),rownames(detP)),] 

keep <- rowSums(detP < 0.01) == ncol(filtered450K)
table(keep)
filtered450K <- filtered450K[keep,]


detP <- arrayRawEpic %>%
  detectionP
detP <- detP[match(featureNames(filteredEPIC),rownames(detP)),]

keep <- rowSums(detP < 0.01) == ncol(filteredEPIC)
table(keep)
filteredEPIC <- filteredEPIC[keep,]

# filter out X and Y chromosomes

annotation <- getAnnotation(filtered450K)
keep <- !(featureNames(filtered450K) %in% 
            annotation$Name[annotation$chr %in% 
                              c("chrX","chrY")])
table(keep)
filtered450K <- filtered450K[keep,]

annotation <- getAnnotation(filteredEPIC)
keep <- !(featureNames(filteredEPIC) %in% 
            annotation$Name[annotation$chr %in% 
                              c("chrX","chrY")])
table(keep)
filteredEPIC <- filteredEPIC[keep,]

# filter out SNPs

filtered450K <- dropLociWithSnps(filtered450K)

filteredEPIC <- dropLociWithSnps(filteredEPIC)

# filter out X reactive probes

xReactiveProbes <- read.csv("https://github.com/sirselim/illumina450k_filtering/blob/master/48639-non-specific-probes-Illumina450k.csv", stringsAsFactors=FALSE)

keep <- !(featureNames(filtered450K) %in% xReactiveProbes$TargetID)
table(keep)
filtered450K <- filtered450K[keep,]

keep <- !(featureNames(filteredEPIC) %in% xReactiveProbes$TargetID)
table(keep)
filteredEPIC <- filteredEPIC[keep,]

# Extract betas

betas450k <- getBeta(filtered450K)
betasEpic <- getBeta(filteredEPIC)

# merge datasets, keep probes available in both platfroms

commonProbes <- intersect(rownames(betasEpic), rownames(betas450k))

betas <- cbind(betasEpic[commonProbes, ], betas450k[commonProbes, ])

# clean column names (remove GSM_xxx prefix)

colnames(betas) <-  colnames(betas) %>%
  str_extract(pattern = "[0-9]*_R[0-9]{2}C[0-9]{2}$")

saveRDS(object = betas, file = "./data/methylation_data_filtered.rds")

### load annotation ----------------------------------------------------

annoFiles <- list.files(path = "./00_christoph_annotation/cleaned/",
                        pattern = ".csv",
                        full.names = TRUE)

anno <- lapply(as.list(annoFiles), read_csv)

anno <- Reduce(f = bind_rows, x = anno)

saveRDS(object = anno, file = "./data/sample_annotation.rds")


