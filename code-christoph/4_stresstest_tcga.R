# Christoph Geisenberger
# github: @cgeisenberger
# last edited 02/10/2022



### Preparation ----------------------------------------------------------------

library(tidyverse)
library(ggplot2)

library(keras)

library(Rtsne)
library(umap)

source("./code-christoph/0_functions.R")
source("./code-christoph/0_branded_colors.R")



# load data for pancreas data set ----------------------------------------------

anno <- readRDS("./input/sample_annotation_umap_purity.rds")
betas <- readRDS(file = "./input/pancreas_betas_everything.rds")
betas <- betas[, anno$array_id]


# load data for TCGA cases -----------------------------------------------------

# processing was performed on server
anno_tcga <- read_csv("./input/annotation/cleaned/annotation_tcga_subsampled.csv")
betas_tcga <- readRDS(file = "./input/betas_tcga_subsampled.rds")
betas_tcga <- betas_tcga[, anno_tcga$basename]


# remove some columns from TCGA annotation 
anno_tcga <- anno_tcga %>% 
  select(basename, barcode, tissue)

anno_tcga %>% 
  group_by(tissue) %>% 
  summarise(n = n())

anno_tcga$tissue[anno_tcga$tissue == "ACC"] <- "AdC"


# load models and model probes -------------------------------------------------

model_probes <- readRDS("./output/model_probes.rds")
nn_model <- load_model_hdf5(filepath = "./output/model_nn.hdf5")
rf_model <- readRDS(file = "./output/model_rf.rds")


# run data through neural network ----------------------------------------------

# predict
input_model <- betas_tcga[model_probes, ]
tcga_scores_nn <- predict(object = nn_model, t(input_model))

# add labels to score matrixs
tcga_labels <- c("ACC", "NORM", "PanNEC", "PanNET", "PB", "PDAC", "SPN")
colnames(tcga_scores) <- tcga_labels

# calculate max score and label for all cases
tcga_max_score <- apply(tcga_scores, 1, max)
tcga_classes <- apply(tcga_scores, 1, function(x) tcga_labels[which.max(x)])

# add data to annotation
anno_tcga <- anno_tcga %>% 
  mutate(nn_class = tcga_classes, 
         nn_max = tcga_max_score,
         nn_score_acc = tcga_scores[, 1], 
         nn_score_norm = tcga_scores[, 2],
         nn_score_nec = tcga_scores[, 3],
         nn_score_net = tcga_scores[, 4], 
         nn_score_pb = tcga_scores[, 5], 
         nn_score_pdac = tcga_scores[, 6], 
         nn_score_spn = tcga_scores[, 7])


anno_tcga <- anno_tcga %>% 
  mutate(mean_beta = apply(betas_tcga, 2, mean))


anno_tcga %>% 
  ggplot(aes(mean_beta, tissue, fill = tissue)) +
  geom_violin() +
  theme_classic(base_size = 20) +
  labs(x = "Avg. Beta", y = NULL) +
  theme(legend.position = "none")


anno_tcga %>% 
  ggplot(aes(nn_max)) +
  geom_histogram(fill = branded_colors2[3], alpha = 0.5, col =  branded_colors2[3]) +
  theme_classic(base_size = 24) +
  labs(x = "Highest Score", y = "Count")


anno_tcga %>% 
  ggplot(aes(tissue, nn_max)) +
  geom_jitter(width = 0.2, col = branded_colors2[3]) +
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = NULL, y = "Highest Score")

anno_tcga %>% 
  group_by(nn_class) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(n, nn_class, fill = nn_class)) +
  geom_col() +
  scale_fill_manual(values = branded_colors1) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none") +
  labs(y = NULL, x = "Predicted class (n)")

anno_tcga %>% 
  group_by(nn_class) %>% 
  summarise(n = n()) %>% 
  ungroup %>% 
  mutate(n = n / sum(n))


test <- anno_tcga %>% 
  select(tissue, nn_class) %>% 
  group_by(tissue, nn_class) %>% 
  summarise(n = n()/25) %>% 
  ungroup %>% 
  pivot_wider(id_cols = tissue, names_from = nn_class, values_from = n)

test2 <- as.matrix(test[, 2:7])
rownames(test2) <- test$tissue
test2[is.na(test2)] <- 0

superheat::superheat(test2, 
                     #order.rows = 1:28,
                     pretty.order.rows = TRUE,
                     bottom.label.text.angle = 90,
                     bottom.label.text.size = 6,
                     bottom.label.col = "white",
                     left.label.col = "white",
                     left.label.text.size = 6,
                     #X.text = round(test2, 2),
                     heat.pal = c("white",  branded_colors2[3]))
pheatmap::pheatmap(test2)



### run data through random forest ---------------------------------------------

# predict
tcga_scores_rf <- predict(object = rf_model, t(betas_tcga[model_probes, ]), type = "prob")
pancreas_scores_rf <- predict(object = rf_model, t(betas[model_probes, ]), type = "prob")

tcga_scores_max <- apply(tcga_scores_rf, 1, max)
pancreas_scores_max <- apply(pancreas_scores_rf, 1, max)


tcga_class_rf <- predict(object = rf_model, t(betas_tcga[model_probes, ]), type = "response")
pancreas_class_rf <- predict(object = rf_model, t(betas[model_probes, ]), type = "response")

dim(tcga_scores_rf)
dim(pancreas_scores_rf)


test_anno <- tibble(dataset = c(rep("TCGA", nrow(anno_tcga)), rep("validation", nrow(anno))), 
                    label = c(anno_tcga$tissue, anno$label), 
                    rf_class = c(as.vector(tcga_class_rf), as.vector(pancreas_class_rf)), 
                    rf_score = c(tcga_scores_max, pancreas_scores_max))

test_anno <- test_anno %>% 
  mutate(correct = ifelse(label == rf_class, "yes", "no"))

test_anno %>% 
  ggplot(aes(rf_score, colour = dataset)) +
  geom_density(linewidth = 2) +
  theme_classic(base_size = 24)
  #facet_wrap(~ rf_class)

cutoff <- seq(from = 0, to = 1, by = 0.01)

cutoff_data <- matrix(data = NaN, nrow = length(cutoff), ncol = 4)
for (i in 1:length(cutoff)){
  cutoff_data[i, 1] <- cutoff[i]
  cutoff_data[i, 2] <- sum(test_anno$rf_score[test_anno$dataset == "TCGA"] > cutoff[i])
  cutoff_data[i, 3] <- sum(test_anno$rf_score[test_anno$dataset == "validation" & test_anno$correct == "yes"] > cutoff[i])
  cutoff_data[i, 4] <- sum(test_anno$rf_score[test_anno$dataset == "validation" & test_anno$correct == "no"] > cutoff[i])
}

colnames(cutoff_data) <-  c("cutoff", "outliers", "correct", "misclassified")

cutoff_data <- as_tibble(cutoff_data)

cutoff_data <- cutoff_data %>% 
  mutate(accuracy = correct/(correct + misclassified))
cutoff_data <- cutoff_data %>% 
  mutate(outliers = outliers/max(outliers), 
         predictable = correct / max(correct))


cutoff_data %>% 
  select(cutoff, outliers, accuracy, predictable) %>% 
  pivot_longer(cols = 2:4, names_to = "statistic", values_to = "value") %>% 
  ggplot(aes(cutoff, value, col = statistic)) +
  geom_line(linewidth = 2) +
  theme_classic(base_size = 24)

test <- randomForest(y = as.factor(c(rep("outlier", ncol(betas_tcga)), rep("real", ncol(betas)))),
                     x = rbind(tcga_scores_rf, pancreas_scores_rf))
test
table(test$predicted[test_anno$dataset == "validation"], test_anno$label[test_anno$dataset == "validation"])

# mahalanobis distance

rf_cov <- cov(pancreas_scores_rf)
rf_means <- colMeans(pancreas_scores_rf)

md_tcga <- mahalanobis(x = tcga_scores_rf , center = rf_means , cov = rf_cov, tol = 1e-20)
md_pancreas <- mahalanobis(x = pancreas_scores_rf , center = rf_means , cov = rf_cov, tol = 1e-20)

test_anno <- test_anno %>% 
  mutate(md = c(md_tcga, md_pancreas))

test_anno %>% 
  ggplot(aes(dataset, md, colour = dataset)) +
  geom_boxplot(linewidth = 1) +
  theme_classic(base_size = 24)

test_anno %>% 
  ggplot(aes(rf_score, colour = correct)) +
  geom_density() +
  facet_wrap(~ rf_class)
  
test_anno %>% 
  ggplot(aes(rf_score, colour = dataset)) +
  geom_density(linewidth = 2) +
  theme_classic(base_size = 24)
