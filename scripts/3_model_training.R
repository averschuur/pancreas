# Christoph Geisenberger
# github: @cgeisenberger
# last edit 26/02/2023 (CG)


# load libraries
library(tidyverse)
library(ggplot2)

reticulate::use_condaenv("r-tensorflow")
library(keras)
library(randomForest)
library(xgboost)
library(caret)

library(doParallel)

source(file = "./scripts/0_helpers.R")



# import annotation and data ---------------------------------------------------

anno <- readRDS("./output/sample_annotation_umap_purity_24102023.rds")
betas <- readRDS(file = "./input/betas_pancreas_everything.rds")

# pick 5,000 most variable probes
top_var_probes <- readRDS(file = "./output/pancreas_top_variable_probes_training_set_24102023.rds")
top_var_probes <- top_var_probes[1:5000]

# filter
betas <- betas[top_var_probes, anno$arrayId]

# split data
train_set <- list()
test_set <- list()

train_set$x <- t(betas[, anno$cohort == "train"])
train_set$y <- as.factor(anno$tumorType[anno$cohort == "train"])

test_set$x <- t(betas[, anno$cohort == "test"])
test_set$y <- as.factor(anno$tumorType[anno$cohort == "test"])

# upsample the training cohort
set.seed(234234)
train_set_upsampled <- upSample(train_set$x, train_set$y, list = TRUE)

# add one-hot converted y for neural networks
train_set_upsampled$onehot <- to_one_hot(train_set_upsampled$y)
test_set$onehot <- to_one_hot(test_set$y)

#saveRDS(object = train_set_upsampled, file = "./output/train_set_upsampled.rds")
#saveRDS(object = test_set, file = "./output/test_set.rds")

#train_set_upsampled <- readRDS("./output/train_set_upsampled_24102023.rds")
#test_set <- readRDS("./output/test_set_24102023.rds")



# default random forest, no tuning ---------------------------------------------

rf_tunegrid <- expand.grid(mtry = floor(sqrt(ncol(train_set$x))))

rf_model <- caret::train(x = train_set_upsampled$x,
                         y = train_set_upsampled$y, 
                         method = "rf",
                         trControl = control,
                         tuneGrid = rf_tunegrid,
                         metric = "Accuracy")

rf_pred <- predict(rf_model, newdata = test_set$x)
rf_cfmatrix <- confusionMatrix(rf_pred, test_set$y)
rf_cfmatrix

#saveRDS(object = rf_model, file = "./output/rf_model_default_24102023.rds")
#saveRDS(object = rf_cfmatrix, file = "./output/rf_confusion_matrix_24102023.rds")



# default xgb, no tuning -------------------------------------------------------

xgb_tunegrid <- expand.grid(nrounds = 100, 
                            max_depth = 6,
                            eta = 0.3,
                            gamma = 1,
                            colsample_bytree = 1,
                            min_child_weight = 1,
                            subsample = 1)

xgb_model <- caret::train(x = train_set_upsampled$x,
                         y = train_set_upsampled$y, 
                         method = "xgbTree",
                         trControl = control,
                         tuneGrid = xgb_tunegrid,
                         metric = "Accuracy")

xgb_pred <- predict(xgb_model, newdata = test_set$x)
xgb_cfmatrix <- confusionMatrix(xgb_pred, test_set$y)
xgb_cfmatrix

#saveRDS(object = xgb_model, file = "./output/xgb_model_default_24102023.rds")
#saveRDS(object = xgb_cfmatrix, file = "./output/xgb_confusion_matrix_24102023.rds")



# neural network ---------------------------------------------------------

# cross-validation parameters
cv_reps <- 3
cv_folds_n <- 10

# set parameters
n_epochs <- 100

# helper function for model
build_model <- function() {
  model <- keras_model_sequential() %>% 
    layer_dense(units = 64, activation = "relu", 
                input_shape = ncol(train_set_upsampled$x)) %>% 
    layer_dense(units = 64, activation = "relu") %>% 
    layer_dense(units = 64, activation = "relu") %>% 
    layer_dense(units = ncol(train_set_upsampled$onehot), activation = 'softmax')
  
  model %>% compile(
    optimizer = "rmsprop", 
    loss = "categorical_crossentropy", 
    metrics = c("accuracy")
  )
}


## run cross-validation -------------------------

perf_hist <- NULL

for(j in 1:cv_reps){
  cat("Processing rep #", j , "\n")
  cv_folds <- cut(seq(1, nrow(train_set_upsampled$x)), 
                  breaks = cv_folds_n, labels = FALSE) %>% sample()
  for (i in 1:cv_folds_n) {
    cat("Processing fold #", i , "\n")
  
    val_indices <- which(cv_folds == i, arr.ind = TRUE)
    val_data <- train_set_upsampled$x[val_indices, ]
    val_labels <- train_set_upsampled$onehot[val_indices, ]
  
    partial_train_data <- train_set_upsampled$x[-val_indices, ]
    partial_train_labels <- train_set_upsampled$onehot[-val_indices, ]
  
    model <- build_model()
  
    history <- model %>% fit(
      as.matrix(partial_train_data), 
      partial_train_labels, 
      validation_data = list(as.matrix(val_data), val_labels),
      epochs = n_epochs, 
      batch_size = 16, 
      verbose = 1
    )
  
    perf <- cbind(j, i, 1:n_epochs, history$metrics$accuracy, history$metrics$val_accuracy)
    perf_hist <- rbind(perf_hist, perf)
  }
}

colnames(perf_hist) <- c("rep", "fold", "epoch", "acc_train", "acc_test")
perf_hist <- as_tibble(perf_hist)

perf_hist %>%
  group_by(epoch) %>%
  summarise(train = mean(acc_train), 
            val = mean(acc_test),
            delta = train-val) %>% 
  ggplot(aes(epoch, delta)) +
  geom_line(lwd = 1) +
  geom_smooth() +
  theme_classic(base_size = 20)

perf_hist %>%
  group_by(epoch) %>%
  summarise(train = sd(acc_train), 
            val = sd(acc_test),
            delta = train-val) %>% 
  pivot_longer(cols = c(train, val)) %>% 
  ggplot(aes(epoch, value, col = name)) +
  geom_smooth() +
  theme_classic(base_size = 20)



## train final model -------------------------

nn_model <- keras_model_sequential() %>% 
  layer_dense(units = 64, activation = "relu", 
              input_shape = ncol(train_set_upsampled$x)) %>%
  layer_dense(units = 64, activation = "relu") %>% 
  layer_dense(units = 64, activation = "relu") %>% 
  layer_dense(units = ncol(train_set_upsampled$onehot), activation = 'softmax')

nn_model %>% compile(
  optimizer = "rmsprop", 
  loss = "categorical_crossentropy", 
  metrics = c("accuracy")
)

nn_model %>% fit(
  as.matrix(train_set_upsampled$x), 
  train_set_upsampled$onehot,
  validation_data = list(as.matrix(test_set$x), test_set$onehot),
  epochs = 50, 
  batch_size = 16
)

# check overall accuracy:
nn_model %>% evaluate(as.matrix(test_set$x), test_set$onehot)

# save model and performance history
save_model_hdf5(object = nn_model, filepath = "./output/nn_model_24102023.hdf5")
saveRDS(object = perf_hist, file = "./output/nn_model_performance_history_24102023.rds")



# compare model performance ------------------------------------------

rf_model <- readRDS(file = "./output/rf_model_default_24102023.rds")
xgb_model <- readRDS(file = "./output/xgb_model_default_24102023.rds")
nn_model <- load_model_hdf5(file = "./output/nn_model_24102023.hdf5")

# rf predictions
rf_pred_class <- predict(rf_model, newdata = t(betas))
rf_pred_scores <- predict(rf_model, newdata = t(betas), type = "prob")
rf_pred_scores_max <- apply(rf_pred_scores, 1, max)
rf_pred_class <- apply(rf_pred_scores, 1, function(x) colnames(rf_pred_scores)[which.max(x)]) %>% as.factor

# xgb predictions
xgb_pred_class <- predict(xgb_model, newdata = t(betas))
xgb_pred_scores <- predict(xgb_model, newdata = t(betas), type = "prob")
xgb_pred_scores_max <- apply(xgb_pred_scores, 1, max)

# nn predictions
nn_pred_scores <- predict(object = nn_model, x = t(betas))
rownames(nn_pred_scores) <- rownames(t(betas))
colnames(nn_pred_scores) <- colnames(rf_pred_scores)
nn_pred_scores_max <- apply(nn_pred_scores , 1, max)
nn_pred_class <- apply(nn_pred_scores, 1, function(x) colnames(nn_pred_scores)[which.max(x)]) %>% as.factor

# add performance to annotation
anno <- anno %>% 
  mutate(pred_nn = nn_pred_class, 
         pred_rf = rf_pred_class, 
         pred_xgb = xgb_pred_class, 
         pred_scores_nn = apply(nn_pred_scores, 1, max),
         pred_scores_rf = apply(rf_pred_scores, 1, max), 
         pred_scores_xgb = apply(xgb_pred_scores, 1, max))

#saveRDS(object = anno, file = "./output/sample_annotation_classifier_performance.rds")
#anno <- readRDS("./output/sample_annotation_classifier_performance.rds")


# confusion matrices
test_indices <- which(anno$cohort == "test", arr.ind = TRUE)
conf_mat <- list(rf_pred_class, xgb_pred_class, nn_pred_class)
conf_mat <- lapply(conf_mat, function(x) x[test_indices])
conf_mat <- lapply(conf_mat, function(x) confusionMatrix(x, test_set$y))
#saveRDS(object = conf_mat, file = "./output/confusion_matrices.rds")
#conf_mat <- readRDS("./output/confusion_matrices.rds")

# plot accuracy per algorithm
perf_acc <- sapply(conf_mat, function(x) x[["overall"]][c(1, 3, 4)])
colnames(perf_acc) <- c("rf", "xgb", "nn")
perf_acc <- as_tibble(t(perf_acc), rownames = "method")

perf_acc  %>% 
  ggplot(aes(x = method, y = Accuracy, fill = method)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = branded_colors2) +
  geom_errorbar(aes(ymin = AccuracyLower, ymax = AccuracyUpper), width = 0.3) +
  theme_bw(base_size = 30) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Accuracy (test cohort)") +
  ylim(0, 1)

# statistics
ANOVA1 <- perf_acc[,c(1:2)]
ANOVA2 <- perf_acc[,c(1,3)]
colnames(ANOVA2) <- c("method", "Accuracy")
ANOVA3 <- perf_acc[,c(1,4)]
colnames(ANOVA3) <- c("method", "Accuracy")
ANOVA <- rbind(ANOVA1, ANOVA2)
ANOVA <- rbind(ANOVA, ANOVA3)

attach(ANOVA)
boxplot(Accuracy~method)


mfit <- lm(Accuracy~factor(method))
anova(mfit)



# plot accuracy across different classes RF
perf_per_class <- conf_mat[[1]]$byClass %>% 
  as_tibble %>% 
  mutate(class = colnames(rf_pred_scores))
colnames(perf_per_class) <- make.names(colnames(perf_per_class))

perf_per_class %>% 
  ggplot(aes(class, Balanced.Accuracy, fill = class)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = branded_colors1) +
  theme_bw(base_size = 30) +
  theme(legend.position = "none")



# project results onto UMAP --------------------------------------------------

anno %>% 
  select(tumorType, absolute, estimate, avg_beta_unfiltered, avg_beta_filtered, conversion, 18:20) %>% 
  GGally::ggpairs()

anno %>% 
  mutate(correct = ifelse(tumorType == pred_rf, "correct", "incorrect")) %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType, shape = correct)) +
  geom_point(size = 3) +
  scale_colour_manual(values = branded_colors1) +
  theme_classic(base_size = 30) +
  labs(x = "Umap 1", y = "Umap 2")


### look into misclassified cases

# RF
incorrect_rf <- anno %>% 
  mutate(correct = ifelse(tumorType == pred_rf, "correct", "incorrect")) %>%
  subset(correct == "incorrect")

correct_rf <- anno %>% 
  mutate(correct = ifelse(tumorType == pred_rf, "correct", "incorrect")) %>%
  subset(correct == "correct")

# NN
incorrect_nn <- anno %>% 
  mutate(correct = ifelse(tumorType == pred_nn, "correct", "incorrect")) %>%
  subset(correct == "incorrect")

# XGB
incorrect_xgb <- anno %>% 
  mutate(correct = ifelse(tumorType == pred_xgb, "correct", "incorrect")) %>%
  subset(correct == "incorrect")



# look at accuracies at different cutoffs for scores ---------------------------

# determine cutoffs for scores
cutoffs <- seq(from = 0.1, to = 0.95, length.out = 18)
real_class <- anno$tumorType[test_indices]

cutoff_nn <- slide_along_cutoff(label_real = real_class, 
                                label_pred = nn_pred_class[test_indices], 
                                scores = nn_pred_scores_max[test_indices], 
                                cutoffs = cutoffs)

cutoff_rf <- slide_along_cutoff(label_real = real_class, 
                                label_pred = rf_pred_class[test_indices], 
                                scores = rf_pred_scores_max[test_indices], 
                                cutoffs = cutoffs)

cutoff_xgb <- slide_along_cutoff(label_real = real_class, 
                                 label_pred = xgb_pred_class[test_indices], 
                                 scores = xgb_pred_scores_max[test_indices], 
                                 cutoffs = cutoffs)

# add name of method
cutoff_nn <- cutoff_nn %>% 
  add_column(method = "NN", .before = TRUE)

cutoff_rf <- cutoff_rf %>% 
  add_column(method = "RF", .before = TRUE)

cutoff_xgb <- cutoff_xgb %>% 
  add_column(method = "XGB", .before = TRUE)

performance_cutoff <- rbind(cutoff_nn, cutoff_rf, cutoff_xgb)
#saveRDS(object = performance_cutoff, file = "./output/performance_cutoff.rds")
#performance_cutoff <- readRDS("./output/performance_cutoff.rds")
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
  scale_color_manual(values = branded_colors1) +
  theme_bw(base_size = 30) +
  theme(legend.title = element_blank(), 
        legend.position = "top") +
  facet_grid(cols = vars(method)) +
  labs(x = "Cutoff", y = "Accuracy/Predictable (%)")
