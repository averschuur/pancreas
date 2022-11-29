

### Prepare methylation data ChanJakelSelenicaDiDomenico ------------------------------------
idats <- list.files(path = "./input/ALL IDATS",
                    full.names = TRUE, 
                    pattern = "_Grn.idat")


# infer platform from filesize
fileSizes <- file.info(idats)$size

idatsEpic <- idats[fileSizes > 10000000]
idats450k <- idats[fileSizes < 10000000]

arrayRawEpic <- read.metharray(basenames = idatsEpic, force = TRUE)
saveRDS(object = arrayRawEpic, file = "./data/preprocessing/arrayRawEpic_ChanJakelSelenicaDiDomenico.rds")
arrayRaw450k <- read.metharray(basenames = idats450k)
saveRDS(object = arrayRaw450k, file = "./data/preprocessing/arrayRaw450k_ChanJakelSelenicaDiDomenico.rds")


# Normalize probes 
filtered450k <- arrayRaw450k %>%
  preprocessNoob(dyeMethod= "single") %>%
  ratioConvert(what = "both", keepCN = TRUE) %>%
  mapToGenome
saveRDS(object = filtered450k, file = "./data/preprocessing/filtered450K_ChanJakelSelenicaDiDomenico.rds")

filteredEPIC <- arrayRawEpic %>%
  preprocessNoob(dyeMethod= "single") %>%
  ratioConvert(what = "both", keepCN = TRUE) %>%
  mapToGenome
saveRDS(object = filteredEPIC, file = "./data/preprocessing/filteredEPIC_ChanJakelSelenicaDiDomenico.rds")

# Extract betas
betas450k <- getBeta(filtered450k)
saveRDS(object = betas450k, file = "./data/preprocessing/betas450k_ChanJakelSelenicaDiDomenico.rds")
betasEpic <- getBeta(filteredEPIC)
saveRDS(object = betasEpic, file = "./data/preprocessing/betasEpic_ChanJakelSelenicaDiDomenico.rds")


###  Filter probes
# remove any probes that have failed in one or more samples
detP <- arrayRaw450k %>%
  detectionP
detP <- detP[match(featureNames(filtered450k),rownames(detP)),] 

keep <- rowSums(detP < 0.01) == ncol(filtered450k)
table(keep)
filtered450k <- filtered450k[keep,]

detP <- arrayRawEpic %>%
  detectionP
detP <- detP[match(featureNames(filteredEPIC),rownames(detP)),]

keep <- rowSums(detP < 0.01) == ncol(filteredEPIC)
table(keep)
filteredEPIC <- filteredEPIC[keep,]


# filter out X and Y chromosomes
annotation <- getAnnotation(filtered450k)
keep <- !(featureNames(filtered450k) %in% 
            annotation$Name[annotation$chr %in% 
                              c("chrX","chrY")])
table(keep)
filtered450K <- filtered450k[keep,]

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
betasfiltered450k <- getBeta(filtered450K)
saveRDS(object = betasfiltered450k, file = "./data/preprocessing/betasfiltered450k_ChanJakelSelenicaDiDomenico.rds")
betasfilteredEpic <- getBeta(filteredEPIC)
saveRDS(object = betasfilteredEpic, file = "./data/preprocessing/betasfilteredEpic_ChanJakelSelenicaDiDomenico.rds")

rm(list=ls())


### Prepare methylation data Endo ------------------------------------
idatsEndo <- list.files(path = "./input/2021_Endo_PDAC, Normal Pancreas/Idat Files_PDAC",
                     full.names = TRUE, 
                     pattern = "_Grn.idat")

arrayRawEndo <- read.metharray(basenames = idatsEndo)
saveRDS(object = arrayRawEndo, file = "./data/preprocessing/arrayRawEndo.rds")

filteredEndo <- arrayRawEndo %>%
  preprocessNoob(dyeMethod= "single") %>%
  ratioConvert(what = "both", keepCN = TRUE) %>%
  mapToGenome
saveRDS(object = filteredEndo, file = "./data/preprocessing/filteredEndo.rds")

betasEndo <- getBeta(filteredEndo)
saveRDS(object = betasEndo, file = "./data/preprocessing/betasEndo.rds")


# remove any probes that have failed in one or more samples
detP <- arrayRawEndo %>%
  detectionP
detP <- detP[match(featureNames(filteredEndo),rownames(detP)),]

keep <- rowSums(detP < 0.01) == ncol(filteredEndo)
table(keep)
filteredEndo <- filteredEndo[keep,]


# filter out X and Y chromosomes
annotation <- getAnnotation(filteredEndo)
keep <- !(featureNames(filteredEndo) %in% 
            annotation$Name[annotation$chr %in% 
                              c("chrX","chrY")])
table(keep)
filteredEndo <- filteredEndo[keep,]


# filter out SNPs
filteredEndo <- dropLociWithSnps(filteredEndo)


# filter out X reactive probes
xReactiveProbes <- read.csv("https://github.com/sirselim/illumina450k_filtering/blob/master/48639-non-specific-probes-Illumina450k.csv", stringsAsFactors=FALSE)

keep <- !(featureNames(filteredEndo) %in% xReactiveProbes$TargetID)
table(keep)
filteredEndo <- filteredEndo[keep,]


# Extract betas
betasfilteredEndo <- getBeta(filteredEndo)
saveRDS(object = betasfilteredEndo, file = "./data/preprocessing/betasfilteredEndo.rds")

rm(list=ls())


### Prepare methylation data Yachida ------------------------------------
idatsYachida <- list.files(path = "./input/Yachida_EPIC_NETNEC/Idat Files",
                             full.names = TRUE, 
                             pattern = "_Grn.idat")

arrayRawYachida <- read.metharray(basenames = idatsYachida, force = TRUE)
saveRDS(object = arrayRawYachida, file = "./data/preprocessing/arrayRawYachida.rds")

filteredYachida <- arrayRawYachida %>%
  preprocessNoob(dyeMethod= "single") %>%
  ratioConvert(what = "both", keepCN = TRUE) %>%
  mapToGenome
saveRDS(object = filteredYachida, file = "./data/preprocessing/filteredYachida.rds")

betasYachida <- getBeta(filteredYachida)
saveRDS(object = betasYachida, file = "./data/preprocessing/betasYachida.rds")

# remove any probes that have failed in one or more samples
detP <- arrayRawYachida %>%
  detectionP
detP <- detP[match(featureNames(filteredYachida),rownames(detP)),]

keep <- rowSums(detP < 0.01) == ncol(filteredYachida)
table(keep)
filteredYachida <- filteredYachida[keep,]

# filter out X and Y chromosomes
annotation <- getAnnotation(filteredYachida)
keep <- !(featureNames(filteredYachida) %in% 
            annotation$Name[annotation$chr %in% 
                              c("chrX","chrY")])
table(keep)
filteredYachida <- filteredYachida[keep,]


# filter out SNPs
filteredYachida <- dropLociWithSnps(filteredYachida)


# filter out X reactive probes
xReactiveProbes <- read.csv("https://github.com/sirselim/illumina450k_filtering/blob/master/48639-non-specific-probes-Illumina450k.csv", stringsAsFactors=FALSE)

keep <- !(featureNames(filteredYachida) %in% xReactiveProbes$TargetID)
table(keep)
filteredYachida <- filteredYachida[keep,]


# Extract betas
betasfilteredYachida <- getBeta(filteredYachida)
saveRDS(object = betasfilteredYachida, file = "./data/preprocessing/betasfilteredYachida.rds")

rm(list=ls())

### Prepare methylation data Benhamida ------------------------------------
idatsBenhamida <- list.files(path = "./input/2022_Benhamida_EPIC_PC_ACC",
                        full.names = TRUE, 
                        pattern = "_Grn.idat")

arrayRawBenhamida <- read.metharray(basenames = idatsBenhamida)
saveRDS(object = arrayRawBenhamida, file = "./data/preprocessing/arrayRawBenhamida.rds")

filteredBenhamida <- arrayRawBenhamida %>%
  preprocessNoob(dyeMethod= "single") %>%
  ratioConvert(what = "both", keepCN = TRUE) %>%
  mapToGenome
saveRDS(object = filteredBenhamida, file = "./data/preprocessing/filteredBenhamida.rds")

betasBenhamida <- getBeta(filteredBenhamida)
saveRDS(object = betasBenhamida, file = "./data/preprocessing/betasBenhamida.rds")


# remove any probes that have failed in one or more samples
detP <- arrayRawBenhamida %>%
  detectionP
detP <- detP[match(featureNames(filteredBenhamida),rownames(detP)),]

keep <- rowSums(detP < 0.01) == ncol(filteredBenhamida)
table(keep)
filteredBenhamida <- filteredBenhamida[keep,]

# filter out X and Y chromosomes
annotation <- getAnnotation(filteredBenhamida)
keep <- !(featureNames(filteredBenhamida) %in% 
            annotation$Name[annotation$chr %in% 
                              c("chrX","chrY")])
table(keep)
filteredBenhamida <- filteredBenhamida[keep,]


# filter out SNPs
filteredBenhamida <- dropLociWithSnps(filteredBenhamida)


# filter out X reactive probes
xReactiveProbes <- read.csv("https://github.com/sirselim/illumina450k_filtering/blob/master/48639-non-specific-probes-Illumina450k.csv", stringsAsFactors=FALSE)

keep <- !(featureNames(filteredBenhamida) %in% xReactiveProbes$TargetID)
table(keep)
filteredBenhamida <- filteredBenhamida[keep,]


# Extract betas
betasfilteredBenhamida <- getBeta(filteredBenhamida)
saveRDS(object = betasfilteredBenhamida, file = "./data/preprocessing/betasfilteredBenhamida.rds")

rm(list=ls())

################################################
### merge unfiltered datasets                ###
### keep probes available in both platforms  ###
################################################

#load data
betasEPIC <- readRDS("./data/preprocessing/betasEpic_ChanJakelSelenicaDiDomenico.rds")
betas450k <- readRDS("./data/preprocessing/betas450k_ChanJakelSelenicaDiDomenico.rds")
betasEndo <- readRDS("./data/preprocessing/betasEndo.rds")
betasYachida <- readRDS("./data/preprocessing/betasYachida.rds")
betasBenhamida <- readRDS("./data/preprocessing/betasBenhamida.rds")

#merge betas
commonProbes <- intersect(rownames(betasEPIC), rownames(betas450k))
betas <- cbind(betasEPIC[commonProbes, ], betas450k[commonProbes, ])
#### 237 samples

# add Endo data
commonProbes <- intersect(rownames(betas), rownames(betasEndo))
betas <- cbind(betas[commonProbes, ], betasEndo[commonProbes, ])
#### 319 samples

# add Yachida data
commonProbes <- intersect(rownames(betas), rownames(betasYachida))
betas <- cbind(betas[commonProbes, ], betasYachida[commonProbes, ])
#### 367 samples

# add Benhamida data
commonProbes <- intersect(rownames(betas), rownames(betasBenhamida))
betas <- cbind(betas[commonProbes, ], betasBenhamida[commonProbes, ])
#### 396 samples

# clean column names (remove GSM_xxx prefix)
colnames(betas) <- colnames(betas) %>%
  str_extract(pattern = "[0-9]*_R[0-9]{2}C[0-9]{2}$")

saveRDS(object = betas, file = "./data/methylation_data.rds")

rm(list=ls())



################################################
### merge filtered datasets                  ###
### keep probes available in both platforms  ###
################################################

#load data
filteredbetasEPIC <- readRDS("./data/preprocessing/betasfilteredEpic_ChanJakelSelenicaDiDomenico.rds")
filteredbetas450k <- readRDS("./data/preprocessing/betasfiltered450k_ChanJakelSelenicaDiDomenico.rds")
filteredbetasEndo <- readRDS("./data/preprocessing/betasfilteredEndo.rds")
filteredbetasYachida <- readRDS("./data/preprocessing/betasfilteredYachida.rds")
filteredbetasBenhamida <- readRDS("./data/preprocessing/betasfilteredBenhamida.rds")

# merge betas
commonProbes <- intersect(rownames(filteredbetasEPIC), rownames(filteredbetas450k))
betas <- cbind(filteredbetasEPIC[commonProbes, ], filteredbetas450k[commonProbes, ])

# add Endo data
commonProbes <- intersect(rownames(betas), rownames(filteredbetasEndo))
betas <- cbind(betas[commonProbes, ], filteredbetasEndo[commonProbes, ])

# add Yachida data
commonProbes <- intersect(rownames(betas), rownames(filteredbetasYachida))
betas <- cbind(betas[commonProbes, ], filteredbetasYachida[commonProbes, ])

# add Benhamida data
commonProbes <- intersect(rownames(betas), rownames(filteredbetasBenhamida))
betas <- cbind(betas[commonProbes, ], filteredbetasBenhamida[commonProbes, ])

# clean column names (remove GSM_xxx prefix)
colnames(betas) <-  colnames(betas) %>%
  str_extract(pattern = "[0-9]*_R[0-9]{2}C[0-9]{2}$")

saveRDS(object = betas, file = "./data/methylation_data_filtered.rds")

rm(list=ls())


