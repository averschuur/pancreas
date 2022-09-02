# Christoph Geisenberger
# github: @cgeisenberger
# last edited 01/09/2022



### Preparation ----------------------------------------------------------------

library(tidyverse)
library(ggplot2)

library(keras)
library(randomForest)
library(xgboost)

library(Rtsne)
library(umap)

source("./code-christoph/0_functions.R")
source("./code-christoph/0_branded_colors.R")

# load annotation and data
anno <- readRDS("./input/sample_annotation_extended.rds")
betas <- readRDS(file = "./input/betas_filtered.rds")

# keep only primaries (PDAC, ACC, SPN, PanNET) and normal pancreas
anno <- anno %>% 
  filter(!label %in% c("MACNEC", "Mixed", "NORM_OTHER")) %>% 
  filter(location %in% c("pancreas", "primary"))

betas <- betas[, anno$array_id]

anno %>% 
  group_by(label) %>% 
  summarise(n = n())

### split sample into training and test cohort -----------------------


# use n = 10 samples per tumor type and study for the training cohort
train_anno <- anno %>% 
  filter(source != "UMCU") %>% 
  group_by(source, label) %>% 
  slice_head(n = 10) %>% 
  ungroup()

# use rest of the data set as test cohort
test_anno <- anno %>% 
  anti_join(y = train_anno)

# training set stats
train_anno %>% 
  group_by(label, source) %>%
  summarise(n = n())

# select train data
train_data <- betas[, train_anno$array_id]

# select most variable probes 
probes_var <- apply(train_data, 1, var)
probes_topvar <- order(probes_var, decreasing = TRUE)[1:5000]

# subset training and test data
train_data <- betas[probes_topvar, train_anno$array_id]
test_data <- betas[probes_topvar, test_anno$array_id]

# define target values, i.e. labels
train_labels <- train_anno$label
test_labels <- test_anno$label

# convert to one-hot encoding for neural network
train_labels_onehot <- to_one_hot(train_labels)
test_labels_onehot <- to_one_hot(test_labels)

# convert to 0-based numeric vector for XGBoost
train_labels_0based <- as.numeric(as.factor(train_labels)) - 1
test_labels_0based <- as.numeric(as.factor(test_labels)) - 1



# determine groups for 5-fold cross-validation ---------------------------------

set.seed(12435)
k <- 3
indices <- sample(1:ncol(train_data))
folds <- cut(1:length(indices), breaks = k, labels = FALSE)


### train neural network ---------------------------------------------------------

# set parameters
num_epochs <- 20
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
    batch_size = 4, 
    verbose = 1
  )
  
  accuracy_history <- c(history$metrics$accuracy, history$metrics$val_accuracy)
  
  all_accuracy_histories <- cbind(all_accuracy_histories, accuracy_history)
}

# assemble performance metrics into dataframe
nn_stats <- tibble(
  epochs = rep(1:num_epochs, 2), 
  measure = rep(c("accuracy", "val_accuracy"), each = 20),
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
  scale_colour_manual(values = branded_colors2) +
  scale_fill_manual(values = branded_colors2)


# train final model
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
  epochs = 50, 
  batch_size = 4
)

# check overall accuracy:
nn_model %>% evaluate(t(test_data), test_labels_onehot)






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
  
  model <- randomForest(x = partial_train_data, y = as.factor(partial_train_labels))
  
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
rf_model <- randomForest(x = t(train_data), y = as.factor(train_labels))
rf_model



### train gradient boosting machines (XGBoost) -----------------------------------

xgb_params <- list("objective" = "multi:softprob",
                   "eval_metric" = "merror",
                   "num_class" = max(train_labels_0based) + 1)

xgb_cv <- xgb.cv(params = xgb_params, data = t(train_data), label = train_labels_0based, 
                 nrounds = 100, nfold = 3, 
                 showsd = TRUE, stratified = TRUE, print_every_n = 10, 
                 early_stop_round = 20, maximize = FALSE, prediction = TRUE)

xgb_model <- xgboost(data = t(train_data), 
                    label = train_labels_0based, 
                    nrounds = 100,
                    subsample = 0.7, 
                    colsample_bytree = 0.7, 
                    params = xgb_params)



### calculate model performance in test data -------------------------------------

test_anno %>% 
  group_by(source) %>% 
  summarise(n = n())


# NEURAL NETWORK 

# scores
pred_nn_scores <- predict(object = nn_model, t(test_data))
colnames(pred_nn_scores) <- colnames(test_labels_onehot)

# classes
pred_nn_classes <- apply(pred_nn_scores, 1, function(x){
  colnames(pred_nn_scores)[which.max(x)]
})

#merge levels of training and test dataset together and print the confusionMatrix
table(prediction = pred_nn_classes, actual = test_anno$label) %>% 
  caret::confusionMatrix()


# RANDOM FOREST

# scores
pred_rf_scores <- predict(object = rf_model, t(test_data), type = "prob")

# classes
pred_rf_classes <- predict(object = rf_model, t(test_data))

#merge levels of training and test dataset together and print the confusionMatrix
table(prediction = pred_rf_classes, actual = test_anno$label) %>% 
  caret::confusionMatrix()


# XGBOOST 

# scores
pred_xgb_scores <- predict(object = xgb_model, newdata = t(test_data))
pred_xgb_scores <- pred_xgb_scores %>% 
  matrix(nrow = 5)
pred_xgb_scores <- t(pred_xgb_scores)
colnames(pred_xgb_scores) <- c("ACC", "NORM", "PanNET", "PDAC", "SPN")

# classes
pred_xgb_class <- apply(pred_xgb_scores, 1, function(x){
  colnames(pred_xgb_scores)[which.max(x)]
})

# merge levels of training and test dataset together and print the confusionMatrix
table(prediction = pred_xgb_class, actual = test_anno$label) %>% 
  caret::confusionMatrix()


# add performance to annotation
test_anno <- test_anno %>% 
  mutate(pred_nn = pred_nn_classes, 
         pred_rf = pred_rf_classes, 
         pred_xgb = pred_xgb_class, 
         pred_scores_nn = apply(pred_nn_scores, 1, max),
         pred_scores_rf = apply(pred_rf_scores, 1, max), 
         pred_scores_xgb = apply(pred_xgb_scores, 1, max))

test_anno %>% 
  mutate(nn = as.integer(label == pred_nn), 
         rf = as.integer(label == pred_rf), 
         xgb = as.integer(label == pred_xgb)) %>% 
  select(nn, rf, xgb) %>% 
  pivot_longer(cols = everything()) %>% 
  group_by(name) %>% 
  summarise(accuracy = sum(value)/n()) %>% 
  ggplot(aes(name, accuracy, fill = name)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = branded_colors2) +
  theme_bw(base_size = 30) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Accuracy (test cohort)") +
  ylim(0, 1)

test_anno %>% 
  mutate(nn_corr = as.integer(label == pred_nn)) %>% 
  group_by(label) %>% 
  summarise(accuracy = sum(nn_corr) / n()) %>% 
  ggplot(aes(label, accuracy, fill = label)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = branded_colors) +
  theme_bw(base_size = 30) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Accuracy (test cohort)")



# project results onto UMAP ------------

class_results <- test_anno %>% 
  select(array_id, starts_with("pred"))

anno_ext <- left_join(anno, class_results)
anno_ext <- anno_ext %>% 
  mutate(dataset = ifelse(is.na(pred_nn), "train", "test"))

anno_ext %>% 
  mutate(correct = as.integer(label == pred_nn)) %>% 
  ggplot(aes(umap_x, umap_y, colour = dataset)) +
  geom_point(size = 3) +
  scale_colour_manual(values = branded_colors2) +
  theme_classic(base_size = 30) +
  labs(x = NULL, y = "Accuracy (test cohort)")


anno_ext %>% 
  mutate(correct = as.factor(label == pred_nn)) %>% 
  ggplot(aes(umap_x, umap_y, colour = correct)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_classic(base_size = 30) +
  labs(x = NULL, y = "Accuracy (test cohort)")

