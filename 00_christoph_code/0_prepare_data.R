

library(tidyverse)
library(minfi)



### prepare methylation data 

idats <- list.files(path = "./00_christoph_data/idat/",
                    full.names = TRUE, 
                    pattern = "_Grn.idat")


# infer platform from filesize

fileSizes <- file.info(idats)$size

idatsEpic <- idats[fileSizes > 10000000]
idats450k <- idats[fileSizes < 10000000]

arrayRawEpic <- read.metharray(basenames = idatsEpic, force = TRUE)
arrayRaw450k <- read.metharray(basenames = idats450k)


# normalize and extract beta values

betasEpic <- arrayRawEpic %>% 
  preprocessNoob %>% 
  getBeta()

betas450k <- arrayRaw450k %>% 
  preprocessNoob %>% 
  getBeta()


# merge datasets, keep probes available in both platfroms

commonProbes <- intersect(rownames(betasEpic), rownames(betas450k))

betas <- cbind(betasEpic[commonProbes, ], betas450k[commonProbes, ])

# clean column names (remove GSM_xxx prefix)

colnames(betas) <-  colnames(betas) %>%
  str_extract(pattern = "[0-9]*_R[0-9]{2}C[0-9]{2}$")

saveRDS(object = betas, file = "./00_christoph_data/methylation_data.rds")


### load annotation ------------------------------------


annoFiles <- list.files(path = "./00_christoph_annotation/cleaned/",
                        pattern = ".csv",
                        full.names = TRUE)

anno <- lapply(as.list(annoFiles), read_csv)

anno <- Reduce(f = bind_rows, x = anno)

saveRDS(object = anno, file = "./00_christoph_annotation/sample_annotation.rds")


