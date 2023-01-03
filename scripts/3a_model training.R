######################
### Model training ###
######################

### Load required packages and sources ------------------------------------------------------------
library(tidyverse)
library(ggplot2)

library(keras)
library(randomForest)
library(xgboost)

library(Rtsne)
library(umap)

source("./scripts/functions.R")
source("./scripts/branded_colors.R")

# load annotation and data
anno <- readRDS("./data/sample_annotation.rds")
betas <- readRDS(file = "./data/methylation_data_filtered.rds")

# remove MACNECs and Mixed tumors
anno <- anno %>% 
  filter(!tumorType %in% c("MACNEC", "Mixed", "?", "Mixed_ACC_NEC", "Mixed_ACC_DA", "Mixed_ACC_DA_NEC")) %>% 
  filter(!location %in% c("acinar cells", "alpha cells", "beta cells", 
                          "ductal cells", "MACNEC normal", "?")) %>%
  filter(!sampleName == "UMCU_ACC2")

betas <- betas[, anno$arrayId]

### split sample into training and test cohort -----------------------

# use n = 10 samples per tumor type and study for the training cohort
train_anno <- anno %>% 
  filter(source != "UMCU") %>%
  filter(location %in% c("primary", "pancreas")) %>%
  filter(source != "yachida" | tumorType != "PanNET") %>%
  group_by(source, tumorType, location) %>% 
  slice_head(n = 10) %>% 
  ungroup()

# training set stats
train_anno %>% 
  group_by(tumorType, source, location) %>%
  summarise(n = n())

# use rest of the data set as test cohort
test_anno <- anno %>% 
  anti_join(y = train_anno)

# select train data
train_data <- betas[, train_anno$arrayId]

# select most variable probes 
top_var_probes <- apply(train_data, 1, var) %>% 
  order(decreasing = TRUE)
top_var_probes <- rownames(train_data)[top_var_probes[1:5000]]
top_var_probes_names <- rownames(train_data)[top_var_probes]
saveRDS(object = top_var_probes_names, file = "./data/model_probes.RDS")

# subset training and test data
train_data <- betas[top_var_probes, train_anno$arrayId]
test_data <- betas[top_var_probes, test_anno$arrayId]

# define target values, i.e. labels
train_labels <- train_anno$tumorType
test_labels <- test_anno$tumorType

# convert to one-hot encoding for neural network
train_labels_onehot <- to_one_hot(train_labels)
test_labels_onehot <- to_one_hot(test_labels)

# convert to 0-based numeric vector for XGBoost
train_labels_0based <- as.numeric(as.factor(train_anno$tumorType)) - 1
test_labels_0based <- as.numeric(as.factor(test_anno$tumorType)) - 1


# determine groups for 5-fold cross-validation ---------------------------------

set.seed(2341324)
k <- 5
indices <- sample(1:ncol(train_data))
folds <- cut(1:length(indices), breaks = k, labels = FALSE)


################
### Training ###
################

### train neural network ---------------------------------------------------------

# set parameters
num_epochs <- 35
all_accuracy_histories <- NULL

# helper function for building the NN
build_model <- function() {
  model <- keras_model_sequential() %>% 
    layer_dense(units = 64, activation = "relu", 
                input_shape = nrow(train_data)) %>% 
    layer_dense(units = 64, activation = "relu") %>% 
    layer_dense(units = 64, activation = "relu") %>% 
    layer_dense(units = ncol(train_labels_onehot), activation = 'softmax')
  
  model %>% compile(
    optimizer = "rmsprop", 
    loss = "categorical_crossentropy", 
    metrics = c("accuracy")
  )
}

# run k-fold crossvalidation 
for (i in 1:k) {
  cat("Processing fold #", i , "\n")
  
  val_indices <- indices[which(folds == i, arr.ind = TRUE)]
  val_data <- t(train_data[, val_indices])
  val_labels <- train_labels_onehot[val_indices, ]
  
  partial_train_data <- t(train_data[, -val_indices])
  partial_train_labels <- train_labels_onehot[-val_indices, ]
  
  model <- build_model()
  
  history <- model %>% fit(
    partial_train_data, 
    partial_train_labels, 
    validation_data = list(val_data, val_labels),
    epochs = num_epochs, 
    batch_size = ncol(train_labels_onehot), 
    verbose = 1
  )
  
  accuracy_history <- c(history$metrics$accuracy, history$metrics$val_accuracy)
  
  all_accuracy_histories <- cbind(all_accuracy_histories, accuracy_history)
}

# assemble performance metrics into dataframe
nn_stats <- tibble(
  epochs = rep(1:num_epochs, 2), 
  measure = rep(c("accuracy", "val_accuracy"), each = 35),
  mean = apply(all_accuracy_histories, 1, mean), 
  sd = apply(all_accuracy_histories, 1, sd)
)

# check overall accuracy:
model %>% evaluate(val_data, val_labels)


nn_stats %>%
  ggplot(aes(epochs, mean, col = measure)) +
  geom_line(lwd = 2) +
  theme_classic(base_size = 20) +
  geom_ribbon(aes(x = epochs, ymin = mean-sd, ymax = mean+sd, fill = measure), alpha = 0.2) +
  scale_colour_manual(values = branded_colors) +
  scale_fill_manual(values = branded_colors)

#### train final model ---------------------------------------------------------
nn_model <- keras_model_sequential() %>% 
  layer_dense(units = 64, activation = "relu", 
              input_shape = nrow(train_data)) %>% 
  layer_dense(units = 64, activation = "relu") %>% 
  layer_dense(units = 64, activation = "relu") %>% 
  layer_dense(units = ncol(train_labels_onehot), activation = 'softmax')

nn_model %>% compile(
  optimizer = "rmsprop", 
  loss = "categorical_crossentropy", 
  metrics = c("accuracy")
)

fit <- nn_model %>% fit(
  t(train_data), 
  train_labels_onehot,
  validation_data = list(t(test_data), test_labels_onehot),
  epochs = 25, 
  batch_size = 7
)

# check overall accuracy:
nn_model %>% evaluate(t(test_data), test_labels_onehot)
# 0.1144942 0.9784173 

fit

png(file='accuracy nn_model.png', bg= "transparent")
plot(fit)
dev.off()

# save model
save_model_hdf5(object = nn_model, filepath = "./output/model_nn_955.hdf5")


### train random forest ----------------------------------------------------------


# initiate object to record classifier performance
rf_accuracy_histories <- NULL

# run k-fold crossvalidation
for (i in 1:k) {
  cat("Processing fold #", i , "\n")
  
  val_indices <- indices[which(folds == i, arr.ind = TRUE)]
  val_data <- t(train_data[, val_indices])
  val_labels <- train_labels[val_indices] %>% as.factor
  
  partial_train_data <- t(train_data[, -val_indices])
  partial_train_labels <- train_labels[-val_indices] %>% as.factor
  
  model <- randomForest(x = partial_train_data, y = as.factor(partial_train_labels), ntree =2000)
  
  predicted_train <- model$predicted
  predicted_val <- predict(model, newdata = val_data)
  
  train_accuracy <- sum(partial_train_labels == predicted_train) / length(partial_train_labels)
  val_accuracy <- sum(as.vector(val_labels) == as.vector(predicted_val)) / length(val_labels)
  
  rf_accuracy_histories <- rbind(rf_accuracy_histories, 
                                 c(train_accuracy, val_accuracy))
}

colnames(rf_accuracy_histories) <- c("train", "val")

apply(rf_accuracy_histories, 2, mean)


# train final model
rf_model <- randomForest(x = t(train_data), y = as.factor(train_labels), ntree=2000)
rf_model

# save model
save_model_hdf5(object = rf_model, filepath = "./output/model_rf_947.hdf5")


### train gradient boosting machines (XGBoost) -----------------------------------

xgb_params <- list("objective" = "multi:softprob",
                   "eval_metric" = "merror",
                   "num_class" = max(train_labels_0based) + 1)

xgb_cv <- xgb.cv(params = xgb_params, data = t(train_data), label = train_labels_0based, 
                 nrounds = 100, nfold = 5, 
                 showsd = TRUE, stratified = TRUE, print_every_n = 10, 
                 early_stop_round = 20, maximize = FALSE, prediction = TRUE)

gb_model <- xgboost(data = t(train_data), 
                    label = train_labels_0based, 
                    nrounds = 100,
                    subsample = 0.7, 
                    colsample_bytree = 0.7, 
                    params = xgb_params)




#####################################
### Model performance calculation ###
#####################################

test_anno %>% 
  group_by(tumorType) %>% 
  summarise(n = n())

### neural network -------------------------------------------------------------
# scores
pred_nn_scores <- predict(object = nn_model, t(test_data))
colnames(pred_nn_scores) <- colnames(test_labels_onehot)

# classes
pred_nn_class <- apply(pred_nn_scores, 1, function(x){
  colnames(pred_nn_scores)[which.max(x)]
})

# merge levels of training and test dataset together and print the confusionMatrix
table(prediction = as.factor(pred_nn_class), actual = as.factor(test_anno$tumorType)) %>% 
  caret::confusionMatrix()


### random forest  -------------------------------------------------------------

# scores
pred_rf_scores <- predict(object = rf_model, t(test_data), type = "prob")

# classes
pred_rf_class <- predict(object = rf_model, t(test_data))

# merge levels of training and test dataset together and print the confusionMatrix
table(prediction = pred_rf_class, actual = test_anno$tumorType) %>% 
  caret::confusionMatrix()

### xg boosting machines -------------------------------------------------------

# scores
pred_xgb_scores <- predict(object = gb_model, newdata = t(test_data))
pred_xgb_scores <- pred_xgb_scores %>% 
  matrix(nrow = 7)
pred_xgb_scores <- t(pred_xgb_scores)
colnames(pred_xgb_scores) <- c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN")

# classes
pred_xgb_class <- apply(pred_xgb_scores, 1, function(x){
  colnames(pred_xgb_scores)[which.max(x)]
})

# merge levels of training and test dataset together and print the confusionMatrix
table(prediction = pred_xgb_class, actual = test_anno$tumorType) %>% 
  caret::confusionMatrix()



###########################################
### algorithm performance visualization ###
###########################################

### add performance to annotation ----------------------------------------------
test_anno <- test_anno %>% 
  mutate(pred_nn = pred_nn_class, 
         pred_rf = pred_rf_class, 
         pred_xgb = pred_xgb_class)



# add performance to annotation
test_anno <- test_anno %>% 
  mutate(pred_nn = pred_nn_class, 
         pred_rf = pred_rf_class, 
         pred_xgb = pred_xgb_class, 
         pred_scores_nn = apply(pred_nn_scores, 1, max),
         pred_scores_rf = apply(pred_rf_scores, 1, max), 
         pred_scores_xgb = apply(pred_xgb_scores, 1, max))

test_anno %>% 
  mutate(nn = as.integer(tumorType == pred_nn), 
         rf = as.integer(tumorType == pred_rf), 
         xgb = as.integer(tumorType == pred_xgb)) %>% 
  select(nn, rf, xgb) %>% 
  pivot_longer(cols = everything()) %>% 
  group_by(name) %>% 
  summarise(accuracy = sum(value)/n()) %>% 
  ggplot(aes(name, accuracy, fill = name)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = branded_colors) +
  theme_bw(base_size = 18) +
  labs(x = NULL, y = "Accuracy (test cohort)")
ggsave("accuracy_models.png", path= "./output/")


x <- test_anno %>% 
  mutate(nn_corr = as.integer(tumorType == pred_nn),
         rf_corr = as.integer(tumorType == pred_rf),
         xgb_corr = as.integer(tumorType == pred_xgb)) %>% 
  group_by(tumorType) %>% 
  summarise(sum(nn_corr) / n(), sum(rf_corr) / n(), sum(xgb_corr) / n())



### dimensionality reduction and UMAP/tSNE for whole cohort ----------------------
#Why is this script included in the model building script?

dr_input <- 1 - cor(betas[top_var_probes, ])
dr_input <- betas[top_var_probes, ] %>% t

# run t-SNE
tsne <- Rtsne(dr_input, perplexity = 15)

# plot t-SNE
anno %>% 
  ggplot(aes(tsne$Y[, 1], tsne$Y[, 2], col = tumorType)) + 
  geom_point(size = 4, alpha = 0.7) +
  #geom_text(data = anno %>% 
  #filter(tumorType == "PanNET"), aes(tsne$Y[, 1], tsne$Y[, 2], label = sampleName)) +
  geom_label(aes(tsne$Y[, 1], tsne$Y[, 2], label = sampleName, alpha = 0.5)) +
  scale_colour_manual(values = branded_colors) +
  labs(x = "tSNE 1", y = "tSNE 2") +
  theme_bw(base_size = 20)
ggsave("tSNE-whole cohort_labeled.png", path= "./output/")

# run UMAP
umap <- umap(dr_input, min_dist = 0.5, n_neighbors = 10)
plot(umap$layout)

# add info to annotation
anno <- anno %>% 
  mutate(tsne_x = tsne$Y[, 1], 
         tsne_y = tsne$Y[, 2], 
         umap_x = umap$layout[, 1], 
         umap_y = umap$layout[, 2])


# show which samples are in test and training cohort
anno <- anno %>% 
  mutate(split = ifelse(arrayId %in% train_anno$arrayId, "train", "test"))

anno <- anno %>% 
  mutate(avg_beta = apply(betas, 2, mean, na.rm = TRUE))


# tumor type (plot UMAP)
anno %>% 
  ggplot(aes(umap_x, umap_y, color = tumorType)) +
  geom_point(size = 4, alpha = 0.7) +
  #geom_label(aes(tsne_x, tsne_y, label = sampleName, alpha = 0.5)) +
  scale_colour_manual(values = branded_colors) +
  theme_bw(base_size = 18) +
  labs(x = "Umap 1", y = "Umap 2")
ggsave("UMAP-whole cohort.png", path= "./output/")

# avg methylation 
anno %>% 
  ggplot(aes(tumorType, avg_beta, fill = tumorType)) +
  geom_violin() +
  scale_fill_manual(values = branded_colors) +
  theme_bw(base_size = 18) +
  labs(x = NULL, y = "Average Methylation") +
  theme(legend.position = "none")
ggsave("AVG methylation-whole cohort.png", path= "./output/")


anno %>% 
  ggplot(aes(tsne_x, tsne_y, color = split)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_fill_manual(values = branded_colors) +
  theme_bw(base_size = 18) +
  labs(x = "t-SNE 1", y = "t-SNE 2")
ggsave("tSNE by split-whole cohort.png", path= "./output/")



# boxplot of MVPs per tumorType #werkt niet
train_data %>%
  ggplot(aes(train_data$tumorType, train_data[top_var_probes[3], ]), color = tumorType) +
  geom_boxplot() +
  scale_fill_manual(values = branded_colors)

### project results onto UMAP ------------

class_results <- test_anno %>% 
  select(arrayId, starts_with("pred"))

anno_ext <- left_join(anno, class_results)
anno_ext <- anno_ext %>% 
  mutate(dataset = ifelse(is.na(pred_nn), "train", "test"))

anno_ext %>% 
  mutate(correct = as.integer(label == pred_nn)) %>% 
  ggplot(aes(umap_x, umap_y, colour = dataset)) +
  geom_point(size = 3) +
  scale_colour_manual(values = branded_colors) +
  theme_classic(base_size = 30) +
  labs(x = NULL, y = "Accuracy (test cohort)")


anno_ext %>% 
  mutate(correct = as.factor(label == pred_nn)) %>% 
  ggplot(aes(umap_x, umap_y, colour = correct)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_classic(base_size = 30) +
  labs(x = NULL, y = "Accuracy (test cohort)")

