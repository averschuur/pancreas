# Christoph Geisenberger
# github: @cgeisenberger
# last edit 24/11/2023


# libraries

library(minfi)
library(tidyverse)
library(umap)
library(clusteval)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

source("./scripts/0_helpers.R")


# import sample annotation
sample_anno <- readRDS(file = "./output/sample_annotation_umap_purity_01112023.rds")

# import pancreas data 
betas <- readRDS("./input/betas_pancreas_filtered.rds")
betas <- betas[, sample_anno$arrayId]

# import most variable probes
topvar_probes <- readRDS("./output/pancreas_top_variable_probes_training_set_01112023.rds")
topvar_probes <- topvar_probes[1:5000]
betas_topvar <- betas[topvar_probes, ]
#betas_topvar <- betas[sample(1:nrow(betas), size = 5000), ]

# import 450K annotation
anno_array <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno_array <- anno_array[rownames(betas), ] %>% 
  as_tibble(rownames = "probe_id")

# clean up the array annotation a little bit
anno_array <- anno_array %>% 
  mutate(topvar = ifelse(probe_id %in% topvar_probes, "topvar", "other"), 
         Relation_to_Island = str_replace(Relation_to_Island, 
                                          pattern = "[N,S]{1}_", replacement = ""))



# Jaccard Index of clustering all vs. topvar probes ----------------------------

# no of clusters (equals no. of tumor types)
cl_k <- 7
set.seed(2315235)

# run clustering for most variable probes
km1 <- kmeans(t(betas_topvar), centers = cl_k)
km1 <- km1$cluster

hc1 <- fastcluster::hclust(dist(t(betas_topvar)))
stump1 <- cutree(hc1, k = cl_k)

# now for all probes
km2 <- kmeans(t(betas), centers = cl_k)
km2 <- km2$cluster

hc2 <- fastcluster::hclust(dist(t(betas)))
stump2 <- cutree(hc2, k = cl_k)

#visualize on umap
sample_anno %>% 
  mutate(cl_km_var = as.factor(km1), 
         cl_km_all = as.factor(km2), 
         cl_hc_var = as.factor(stump1), 
         cl_hc_all = as.factor(stump2)) %>% 
  pivot_longer(cols = starts_with("cl"), values_to = "Cluster") %>% 
  ggplot(aes(umap_x, umap_y, col = Cluster)) +
  geom_point(size = 3) +
  facet_wrap(facets = vars(name)) +
  theme_minimal(base_size = 18) +
  labs(x = "Umap 1", y = "Umap 2")
ggsave(filename = "./plots/umap_clustering_topvar_vs_all.pdf", width = 10, height = 10)


# measure overlap between clustering and tumor types
cs_hc1 <- cluster_similarity(stump1, sample_anno$tumorType, similarity = "jaccard")
cs_hc2 <- cluster_similarity(stump2, sample_anno$tumorType, similarity = "jaccard")
cluster_similarity(stump1, stump2, similarity = "jaccard")

cs_km1 <- cluster_similarity(km1, sample_anno$tumorType, similarity = "jaccard")
cs_km2 <- cluster_similarity(km2, sample_anno$tumorType, similarity = "jaccard")
cluster_similarity(km1, km2, similarity = "jaccard")

# combine data
cl_eval <- tibble(probes = c("topvar", "all"), hc = c(cs_hc1, cs_hc2), km = c(cs_km1, cs_km2))
cl_eval
write_csv(x = cl_eval, file = "./output/cluster_eval_topvar_vs_all.csv")

table(sample_anno$tumorType, stump1) %>% apply(., 2, function(x) max(x)) %>% sum/321
table(stump1, stump2)# %>% apply(., 2, function(x) max(x)) %>% sum/321
table(sample_anno$tumorType, km1) %>% apply(., 2, function(x) max(x)) %>% sum/321
table(sample_anno$tumorType, km2) %>% apply(., 2, function(x) max(x)) %>% sum/321



# investigate distribution of most variable probes -----------------------------

# extract probe locations
topvar_location <- anno_array %>% 
  select(topvar, Relation_to_Island) %>% 
  group_by(topvar, Relation_to_Island) %>% 
  summarise(n = n()) %>%
  group_by(topvar) %>% 
  mutate(prop = n / sum(n) *100) %>% 
  ungroup
#saveRDS(object = topvar_location, file = "./output/revision_reviewer2_topvar_locations.rds")

# plot
topvar_location %>% 
  ggplot(aes(Relation_to_Island, prop, fill = topvar)) +
  geom_col(position = "dodge") +
  theme_classic(base_size = 24) +
  labs(x = NULL, y = "Proportion of probes (%)")

# Chi-Square test (p-val = 4.281e-08)
topvar_location %>% 
  pivot_wider(id_cols = Relation_to_Island, names_from = topvar, values_from = n) %>% 
  select(other, topvar) %>% 
  chisq.test()


# what about gene bodies ? -----------------------------------------------------

# extract probes which are unequivocally associated with one gene
gb <- anno_array %>% 
  filter(UCSC_RefGene_Accession != "") %>% 
  filter(str_detect(UCSC_RefGene_Accession, pattern = ";", negate = TRUE))

# count and normalize w.r.t. to gene position
gb_counts <- gb %>% 
  mutate(cpg_island = ifelse(Relation_to_Island == "Island", "CpG Island", "Non-Island")) %>% 
  select(UCSC_RefGene_Group, cpg_island, topvar) %>% 
  dplyr::rename("location" = UCSC_RefGene_Group) %>% 
  group_by(location, cpg_island, topvar) %>% 
  summarise(n = n()) %>% 
  group_by(topvar) %>% 
  mutate(prop = n/sum(n) * 100)

gb_counts %>% 
  ggplot(aes(location, prop, fill = topvar)) +
  geom_col(position = "dodge") +
  theme_classic(base_size = 24) +
  facet_wrap(nrow = 2, facets = vars(cpg_island), scales = "free_y") +
  labs(x = NULL, y = "Proportion of probes (%)")



# look at differential methylation ---------------------------------------------

prom_anno <- gb %>% 
  filter(Relation_to_Island == "Island", 
         str_detect(UCSC_RefGene_Group, pattern = "TSS")) %>% 
  group_by(UCSC_RefGene_Accession) %>% 
  mutate(n = n()) %>% 
  filter(n > 3)

betas_prom <- betas[prom_anno$probe_id, sample_anno$arrayId]
betas_prom_avg <- apply(betas_prom, 2, FUN = function(x) 
  tapply(x, INDEX = prom_anno$UCSC_RefGene_Name, FUN = mean, na.rm = TRUE))

prom_umap <- umap(t(betas_prom_avg))

sample_anno <- sample_anno %>% 
  mutate(umap_prom_x = prom_umap$layout[, 1], 
         umap_prom_y = prom_umap$layout[, 2])

sample_anno %>% 
  ggplot(aes(umap_prom_x, umap_prom_y, col = tumorType)) +
  geom_point(size = 3) +
  theme_minimal(base_size = 18) +
  labs(x = "Umap 1", y = "Umap 2")

diff_meth_groups <- c("PanNET", "NORM")

anno_comp <- sample_anno %>% 
  filter(tumorType %in% diff_meth_groups)

diff_meth <- genefilter::rowttests(betas_prom_avg[, anno_comp$arrayId], 
                                        fac = as.factor(anno_comp$tumorType))
diff_meth$mean_pnet = apply(betas_prom_avg[, sample_anno$tumorType == diff_meth_groups[1]], 1, mean)
diff_meth$mean_norm = apply(betas_prom_avg[, sample_anno$tumorType == diff_meth_groups[2]], 1, mean)

diff_meth <- diff_meth %>% 
  as_tibble(rownames = "gene_symbol") %>% 
  mutate(pval_log = -log10(p.value))

sum(diff_meth$p.value < 0.001)

diff_meth %>% 
  mutate(filtered_label = ifelse(pval_log > 10, gene_symbol, "")) %>% 
  ggplot(aes(dm, pval_log)) +
  geom_point() +
  geom_text(aes(label = filtered_label), nudge_y = 0.5) +
  theme_bw(base_size = 18) +
  labs(x = "Methylation change (diff. beta)", y = "Significance (-log10 p-value)")

diff_meth %>% 
  ggplot(aes(dm)) +
  geom_histogram(bins = 200) +
  theme_bw(base_size = 18) +
  labs(x = "Methylation change (diff. beta)", y = "Significance (-log10 p-value)")

write_csv(x = diff_meth, file = "./output/diff_meth_PanNET_vs_Norm.csv")
