# Validation


### preparation ----------------------------------------------------------------
# libraries
library(tidyverse)
library(Rtsne)
library(umap)
library(RFpurify)

# load data and annotation 
betasUMCU <- readRDS("./data/methylation_data.rds")
annoUMCU <- readRDS("./00_annotation/sample_annotation.rds")

# keep only primaries of PNETs, ACCs and SPNs
annoUMCU <- annoUMCU %>% 
  filter(!tumorType %in% c("MACNEC", "Mixed")) %>% 
  filter(source == "UMCU")

# filter for samples where annotation is available
betasUMCU <- betasUMCU[, annoUMCU$arrayId]


### estimate tumor purity ------------------------------------------------------

# get purity estimates 
absolute <- RFpurify::predict_purity_betas(betas = betasUMCU, method = "ABSOLUTE")
estimate <- RFpurify::predict_purity_betas(betas = betasUMCU, method = "ESTIMATE")

# collect data into tibble
annoUMCU <- annoUMCU %>% 
  mutate(absolute = absolute,
         estimate = estimate,
         avg_beta = apply(betasUMCU, 2, mean, na.rm = TRUE))

# make some basic plots
annoUMCU %>% 
  ggplot(aes(absolute, estimate)) +
  geom_point(size = 3) +
  geom_smooth() +
  theme_bw(base_size = 18) +
  labs(x = "Purity (ABSOLUTE)", y = "Purity (ESTIMATE)")

annoUMCU %>% 
  ggplot(aes(tumorType, estimate, col = tumorType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(shape=tumorType), position = position_jitterdodge(0.5), alpha=0.9) +
  #geom_jitter(aes(fill=tumorType, shape=tumorType)) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  scale_colour_manual(values = branded_colors) +
  theme_bw(base_size = 18)
ggsave("Figure 3_tumor purity per tumorType.png", path= "./output/")


### some descriptive statistics ------------------------------------------------
annoUMCU %>% 
  ggplot(aes(tumorType)) +
  geom_bar(aes(fill = location)) +
  scale_fill_manual(values = branded_colors) +
  theme_bw(base_size = 18)
ggsave("Figure 3_count tumorType per location.png", path= "./output/")

annoUMCU %>% 
  ggplot(aes(tumorType, avg_beta, col = tumorType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_colour_manual(values = branded_colors) +
  theme_bw(base_size = 18)
ggsave("Figure 3_average methylation.png", path= "./output/")



### unsupervised analysis: tSNE & UMAP ----------------------------------

# use most variable probes obtained from trainingset


# calculate correlation matrix, use 1-cor as dissimlarity measure
cor_matrix <- 1 - cor(betasUMCU[top_var_probes, ])

# run tSNE
tsne <- Rtsne(X = cor_matrix, perplexity = 5)

# run UMAP
umap <- umap(d = t(betasUMCU[top_var_probes, ]))

annoUMCU <- annoUMCU %>% 
  mutate(tsne_x = tsne$Y[, 1],
         tsne_y = tsne$Y[, 2], 
         umap_x = umap$layout[, 1], 
         umap_y = umap$layout[, 2])


# plot t-SNE map
annoUMCU %>% 
  ggplot(aes(tsne_x, tsne_y, col = tumorType)) + 
  geom_point(size = 4, alpha = 0.7) +
  #geom_text(data = anno %>% 
  #filter(tumorType == "PanNET"), aes(tsne_x, tsne_y, label = sampleName)) +
  #geom_label(aes(tsne_x, tsne_y, label = sampleName, alpha = 0.5)) +
  geom_label_repel(data=annotation2, aes(x=x, y=y, label=label), ylim = -75, xlim = 25, color= "#da4167", size = 4, alpha = 0.7) +
  scale_colour_manual(values = branded_colors) +
  labs(x = "tSNE 1", y = "tSNE 2") +
  theme_bw(base_size = 20)
ggsave("Figure5_tSNE-UMCU_SPNcases.png", path= "./output/")

# plot UMAP
annoUMCU %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType, size = 1/estimate)) + 
  geom_point(alpha = 0.7) +
  #geom_label(aes(tsne_x, tsne_y, label = sampleName, alpha = 0.5)) +
  scale_colour_manual(values = branded_colors) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_bw(base_size = 20)
ggsave("Figure3_UMAP-UMCU.png", path= "./output/")

### unsupervised analysis: tSNE & UMAP for cases reports    -------------------
# set annotation
annotation <- data.frame(
  x= c(-79.985050, -43.875610),
  y = c(-18.166953, -42.981255),
  label = c("primary", "metastasis")
)

# plot t-SNE map
annoUMCU %>% 
  ggplot(aes(tsne_x, tsne_y, col = tumorType)) + 
  geom_point(size = 4, alpha = 0.7) +
  #geom_text(data = anno %>% 
  #filter(tumorType == "PanNET"), aes(tsne_x, tsne_y, label = sampleName)) +
  #geom_label(aes(tsne_x, tsne_y, label = sampleName, alpha = 0.5)) +
  geom_label_repel(data=annotation2, aes(x=x, y=y, label=label), ylim = -50, color= "#f6d8ae", size = 4, alpha = 0.7) +
  scale_colour_manual(values = branded_colors) +
  labs(x = "tSNE 1", y = "tSNE 2") +
  theme_bw(base_size = 20)
ggsave("Figure5_tSNE-UMCU_ACCcases.png", path= "./output/")

# set annotation
annotation2 <- data.frame(
  x= -50.430093,
  y = -66.676524,
  label = "primary"
)

# plot t-SNE map
annoUMCU %>% 
  ggplot(aes(tsne_x, tsne_y, col = tumorType)) + 
  geom_point(size = 4, alpha = 0.7) +
  #geom_text(data = anno %>% 
  #filter(tumorType == "PanNET"), aes(tsne_x, tsne_y, label = sampleName)) +
  #geom_label(aes(tsne_x, tsne_y, label = sampleName, alpha = 0.5)) +
  geom_label_repel(data=annotation2, aes(x=x, y=y, label=label), ylim = -75, xlim = 25, color= "#da4167", size = 4, alpha = 0.7) +
  scale_colour_manual(values = branded_colors) +
  labs(x = "tSNE 1", y = "tSNE 2") +
  theme_bw(base_size = 20)
ggsave("Figure5_tSNE-UMCU_SPNcases.png", path= "./output/")


### heatmap --------------------------------------------------------------------

heat <- pheatmap::pheatmap(cor_matrix, labels_row = annoUMCU$sampleName)
save_pheatmap_pdf(heat, "./output/heatmap_UMCU.pdf")


### calculate model performance in validation set -------------------------------------

validation_data <- betasUMCU[top_var_probes, annoUMCU$arrayId]

# define target values, i.e. labels
validation_labels <- annoUMCU$tumorType

# convert to one-hot encoding for neural network
validation_labels_onehot <- to_one_hot(validation_labels)

# NEURAL NETWORK 

# scores
pred_nn_scores <- predict(object = nn_model, t(validation_data))
colnames(pred_nn_scores) <- colnames(validation_labels_onehot)

# classes
pred_nn_classes <- apply(pred_nn_scores, 1, function(x){
  colnames(pred_nn_scores)[which.max(x)]
})

#merge levels of training and test dataset together and print the confusionMatrix
u <- union(pred_nn_scores, pred_nn_classes)
t <- table(prediction = pred_nn_classes, actual = annoUMCU$tumorType)
confusionMatrix(t)
