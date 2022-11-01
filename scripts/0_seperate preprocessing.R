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

# add Endo data
idatsEndo <- list.files(path = "./input/2021_Endo_PDAC, Normal Pancreas/Idat Files_PDAC",
                     full.names = TRUE, 
                     pattern = "_Grn.idat")

arrayRawEndo <- read.metharray(basenames = idatsEndo)

filteredEndo <- arrayRawEndo %>%
  preprocessNoob(dyeMethod= "single") %>%
  ratioConvert(what = "both", keepCN = TRUE) %>%
  mapToGenome

betasEndo <- getBeta(filteredEndo)

### Normalize probes -----------------------------------------------

filtered450K <- arrayRaw450k %>%
  preprocessNoob(dyeMethod= "single") %>%
  ratioConvert(what = "both", keepCN = TRUE) %>%
  mapToGenome

filteredEPIC <- arrayRawEpic %>%
  preprocessNoob(dyeMethod= "single") %>%
  ratioConvert(what = "both", keepCN = TRUE) %>%
  mapToGenome

# Extract betas

betas450k <- getBeta(filtered450K)
betasEpic <- getBeta(filteredEPIC)


# merge datasets, keep probes available in both platfroms

commonProbes <- intersect(rownames(betasEpic), rownames(betas450k))

betas <- cbind(betasEpic[commonProbes, ], betas450k[commonProbes, ])

#### 237 samples

# add Endo data

commonProbes <- intersect(rownames(betas), rownames(betasEndo))
betas <- cbind(betas[commonProbes, ], betasEndo[commonProbes, ])

#### 319 samples


# clean column names (remove GSM_xxx prefix)

colnames(betas) <- colnames(betas) %>%
  str_extract(pattern = "[0-9]*_R[0-9]{2}C[0-9]{2}$")

saveRDS(object = betas, file = "./data/methylation_data.rds")

### Filer probes ----------------------------------------------------
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

detP <- arrayRawEndo %>%
  detectionP
detP <- detP[match(featureNames(filteredEndo),rownames(detP)),]

keep <- rowSums(detP < 0.01) == ncol(filteredEndo)
table(keep)
filteredEndo <- filteredEndo[keep,]

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

annotation <- getAnnotation(filteredEndo)
keep <- !(featureNames(filteredEndo) %in% 
            annotation$Name[annotation$chr %in% 
                              c("chrX","chrY")])
table(keep)
filteredEndo <- filteredEndo[keep,]

# filter out SNPs

filtered450K <- dropLociWithSnps(filtered450K)

filteredEPIC <- dropLociWithSnps(filteredEPIC)

filteredEndo <- dropLociWithSnps(filteredEndo)

# filter out X reactive probes

xReactiveProbes <- read.csv("https://github.com/sirselim/illumina450k_filtering/blob/master/48639-non-specific-probes-Illumina450k.csv", stringsAsFactors=FALSE)

keep <- !(featureNames(filtered450K) %in% xReactiveProbes$TargetID)
table(keep)
filtered450K <- filtered450K[keep,]

keep <- !(featureNames(filteredEPIC) %in% xReactiveProbes$TargetID)
table(keep)
filteredEPIC <- filteredEPIC[keep,]

keep <- !(featureNames(filteredEndo) %in% xReactiveProbes$TargetID)
table(keep)
filteredEndo <- filteredEndo[keep,]

# Extract betas

betas450k <- getBeta(filtered450K)
betasEpic <- getBeta(filteredEPIC)
betasEndo <- getBeta(filteredEndo)

# merge datasets, keep probes available in both platfroms

commonProbes <- intersect(rownames(betasEpic), rownames(betas450k))

betas <- cbind(betasEpic[commonProbes, ], betas450k[commonProbes, ])

# add Endo data

commonProbes <- intersect(rownames(betas), rownames(betasEndo))
betas <- cbind(betas[commonProbes, ], betasEndo[commonProbes, ])

# clean column names (remove GSM_xxx prefix)

colnames(betas) <-  colnames(betas) %>%
  str_extract(pattern = "[0-9]*_R[0-9]{2}C[0-9]{2}$")

saveRDS(object = betas, file = "./data/methylation_data_filtered.rds")

