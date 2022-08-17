# Christoph Geisenberger
# Department of Pathology
# LMU Munich
# github: @cgeisenberger
# last edit 02/08/2022



# libraries
library(tidyverse)

library(keras)
library(randomForest)
library(xgboost)

library(Rtsne)
library(umap)

source("./00_christoph_code/functions.R")
source("./00_christoph_code/0_branded_colors.R")


# load annotation and data
anno <- readRDS(file = "./00_christoph_annotation/sample_annotation.rds")

# remove MACNECs and Mixed tumors
anno <- anno %>% 
  filter(!tumorType %in% c("MACNEC", "Mixed")) %>% 
  filter(!location %in% c("acinar cells", "alpha cells", "beta cells", 
                          "ductal cells", "MACNEC normal"))

# remove
data <- readRDS(file = "./00_christoph_data/methylation_data.rds")
data <- data[, anno$arrayId]



# split sample into training, validation and test cohort -----------------------

# use n = 10 samples per tumor type and study
train_anno <- anno %>% 
  filter(source != "UMCU") %>% 
  filter(tumorType %in% c("ACC", "normal", "PanNET", "SPN")) %>% 
  filter(location %in% c("primary", "pancreas", "acc normal")) %>% 
  group_by(source, tumorType, location) %>% 
  slice_head(n = 10) %>% 
  ungroup()

nrow(train_anno) # n = 76 total

# use rest of the data set as validation data
test_anno <- anno %>% 
  filter(!arrayId %in% train_anno$arrayId)

# look at statistics
train_anno %>% 
  group_by(source, tumorType, location) %>%
  summarise(n = n())

# select train data
train_data <- data[, train_anno$arrayId]

# select most variable probes 
top_var_probes <- apply(train_data, 1, var) %>% 
  order(decreasing = TRUE)
top_var_probes <- rownames(train_data)[top_var_probes[1:5000]]

# subset training and test data
train_data <- data[top_var_probes, train_anno$arrayId]
test_data <- data[top_var_probes, test_anno$arrayId]

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




# train neural network ---------------------------------------------------------

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
    layer_dense(units = 4, activation = 'softmax')
  
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

# train final model
nn_model <- keras_model_sequential() %>% 
  layer_dense(units = 64, activation = "relu", 
              input_shape = nrow(train_data)) %>% 
  layer_dense(units = 64, activation = "relu") %>% 
  layer_dense(units = 64, activation = "relu") %>% 
  layer_dense(units = 4, activation = 'softmax')

nn_model %>% compile(
  optimizer = "rmsprop", 
  loss = "categorical_crossentropy", 
  metrics = c("accuracy")
)

nn_model %>% fit(
  t(train_data), 
  train_labels_onehot,
  validation_data = list(t(test_data), test_labels_onehot),
  epochs = 25, 
  batch_size = 4
)



# train random forest ----------------------------------------------------------


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



# train gradient boosting machines (XGBoost) -----------------------------------

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



# calculate model performance in test data -------------------------------------


test_anno %>% 
  group_by(source) %>% 
  summarise(n = n())


# NEURAL NETWORK 

# scores
pred_nn_scores <- predict(object = nn_model, t(test_data))
colnames(pred_nn_scores) <- colnames(test_labels_onehot)
pred_nn_scores_max <- apply(pred_nn_scores, 1, max)

# classes
pred_nn_classes <- apply(pred_nn_scores, 1, function(x){
  colnames(pred_nn_scores)[which.max(x)]
})


# RANDOM FOREST

# scores
pred_rf_scores <- predict(object = rf_model, t(test_data), type = "prob")
pred_rf_scores_max <- apply(pred_rf_scores, 1, max)

# classes
pred_rf_classes <- predict(object = rf_model, t(test_data))


# XGBOOST 

# scores
pred_xgb_scores <- predict(object = gb_model, newdata = t(test_data))
pred_xgb_scores <- pred_xgb_scores %>% 
  matrix(nrow = 4)
pred_xgb_scores <- t(pred_xgb_scores)
colnames(pred_xgb_scores) <- c("ACC", "normal", "PanNET", "SPN")
pred_xgb_scores_max <- apply(pred_xgb_scores, 1, max)

# classes
pred_xgb_classes <- apply(pred_xgb_scores, 1, function(x){
  colnames(pred_xgb_scores)[which.max(x)]
})


# add performance to annotation
test_anno <- test_anno %>% 
  mutate(pred_nn = pred_nn_classes, 
         pred_rf = pred_rf_classes, 
         pred_xgb = pred_xgb_classes, 
         score_nn = pred_nn_scores_max, 
         score_rf = pred_rf_scores_max, 
         score_xgb = pred_xgb_scores_max)

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

test_anno %>% 
  mutate(nn_corr = as.integer(tumorType == pred_nn)) %>% 
  group_by(source) %>% 
  summarise(sum(nn_corr) / n())



# dimensionality reduction and UMAP/tSNE for whole cohort ----------------------

dr_input <- 1 - cor(data[top_var_probes, ])
dr_input <- data[top_var_probes, ] %>% t

# run t-SNE
tsne <- Rtsne(dr_input, perplexity = 15)

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
  mutate(avg_beta = apply(data, 2, mean, na.rm = TRUE))


# tumor type
anno %>% 
  ggplot(aes(umap_x, umap_y, color = tumorType)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_colour_manual(values = branded_colors) +
  theme_bw(base_size = 18) +
  labs(x = "Umap 1", y = "Umap 2") +
  theme(legend.position = "none")

# avg methylation 
anno %>% 
  ggplot(aes(tumorType, avg_beta, fill = tumorType)) +
  geom_violin() +
  scale_fill_manual(values = branded_colors) +
  theme_bw(base_size = 18) +
  labs(x = NULL, y = "Average Methylation (beta)") +
  theme(legend.position = "none")


anno %>% 
  ggplot(aes(tsne_x, tsne_y, color = split)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_color_manual(values = c("skyblue2", "salmon1")) +
  theme_bw(base_size = 18) +
  labs(x = "t-SNE 1", y = "t-SNE 2")



# heatmap of correlation matrix
pheatmap::pheatmap(dr_input, 
                   labels_row = anno$tumorType)
boxplot(train_data[top_var[3], ] ~ train_anno$tumorType)



### assess implementation of cutoff for classifier scores ----------------------

# get true classes
cutoff_true_classes <- test_anno$tumorType

# get highest score for each method
cutoff_scores_nn <- apply(pred_nn_scores, 1, max)
cutoff_scores_rf <- apply(pred_rf_scores, 1, max)
cutoff_scores_xgb <- apply(pred_xgb_scores, 1, max)

# determine cutoffs for scores
cutoffs <- seq(from = 0.3, to = 0.95, length.out = 14)


cutoff_nn <- slide_along_cutoff(label_real = cutoff_true_classes, 
                                label_pred = pred_nn_classes, 
                                scores = cutoff_scores_nn, 
                                cutoffs = cutoffs)

cutoff_rf <- slide_along_cutoff(label_real = cutoff_true_classes, 
                                label_pred = pred_rf_classes, 
                                scores = cutoff_scores_rf, 
                                cutoffs = cutoffs)

cutoff_xgb <- slide_along_cutoff(label_real = cutoff_true_classes, 
                                label_pred = pred_xgb_classes, 
                                scores = cutoff_scores_xgb, 
                                cutoffs = cutoffs)

# add name of method
cutoff_nn <- cutoff_nn %>% 
  add_column(method = "neuralNet", .before = TRUE)

cutoff_rf <- cutoff_rf %>% 
  add_column(method = "randomForest", .before = TRUE)

cutoff_xgb <- cutoff_xgb %>% 
  add_column(method = "xgBoost", .before = TRUE)

performance_cutoff <- rbind(cutoff_nn, cutoff_rf, cutoff_xgb)
rm(cutoff_nn, cutoff_rf, cutoff_xgb)



# plot cutoff vs. predictable / accuracy

performance_cutoff %>% 
  #filter(method == "neuralNet") %>% 
  pivot_longer(cols = accuracy:predictable, 
              names_to = "statistic") %>% 
  ggplot(aes(cutoff, value, color = statistic)) +
  geom_point(size = 3) +
  geom_line() +
  ylim(c(0.5, 1)) +
  scale_color_manual(values = branded_colors) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none") +
  facet_grid(cols = vars(method)) +
  labs(x = "Cutoff", y = "Accuracy/Predictable (%)")

test_anno %>% 
  select(source, arrayType, tumorType, location, score_nn, score_rf, score_xgb) %>% 
  mutate(location = ifelse(location == "primary", "primary", "other")) %>% 
  pivot_longer(cols = starts_with("score"), names_to = "method", values_to = "score") %>% 
  mutate(method = str_extract(string = method, pattern = "[^_]*$")) %>% 
  ggplot(aes(score, method, col = method)) +
  geom_jitter(size = 3, width = 0.1) +
  scale_color_manual(values = branded_colors) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none") +
  labs(x = "Score", y = NULL)

test_anno %>% 
  select(source, arrayType, tumorType, location, score_nn, score_rf, score_xgb) %>% 
  mutate(location = ifelse(location == "primary", "primary", "other")) %>% 
  pivot_longer(cols = starts_with("score"), names_to = "method", values_to = "score") %>% 
  mutate(method = str_extract(string = method, pattern = "[^_]*$")) %>% 
  ggplot(aes(source, score, col = tumorType)) +
  geom_jitter(size = 3, width = 0.1) +
  scale_color_manual(values = branded_colors) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top") +
  labs(x = "Score", y = NULL)

  
  
