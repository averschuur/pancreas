### Model performance remocing outliers from oancreatic dataset

## 1: identify number of incorrectly classified tumors classified as outlier ---------------
# load data
anno <- readRDS("./output/sample_annotation_classifier_performance_01112023.rds")
rf_data <- readRDS(file = "./output/outlier_det_results_01112023.rds")

# select data
anno_outliers <- rf_data %>%
  filter(class_char == "pancreas")
anno_outliers <- cbind(anno, anno_outliers$od_class)
anno_outliers <- anno_outliers %>% 
  mutate(correct = ifelse(tumorType == pred_rf, "correct", "incorrect"))
colnames(anno_outliers) <-  c("arrayId","source","sampleName","arrayType","tumorType","location","cohort","absolute","estimate",
                              "avg_beta_unfiltered","conversion","avg_beta_filtered","umap_x","umap_y","pred_nn","pred_rf","pred_xgb",
                              "pred_scores_nn","pred_scores_rf","pred_scores_xgb","od_class","correct")

# identify number of incorrectly classified tumors classified as outlier
anno_outliers %>%
  group_by(od_class, correct) %>%
  summarise(n = n())



## 2: calculate classifier performances excluding outliers --------------------------------
# load data
anno <- readRDS("./output/sample_annotation_umap_purity_01112023.rds")
betas <- readRDS(file = "./input/betas_pancreas_everything.rds")
rf_data <- readRDS(file = "./output/outlier_det_results_01112023.rds")

# select data
anno_outliers <- rf_data %>%
  filter(class_char == "pancreas")
anno_outliers <- cbind(anno, anno_outliers$od_class)
colnames(anno_outliers) <-  c("arrayId","source","sampleName","arrayType","tumorType","location","cohort","absolute","estimate",
                              "avg_beta_unfiltered","conversion","avg_beta_filtered","umap_x","umap_y","od_class")

# remove outliers
anno_outliers <- anno_outliers %>%
  filter(od_class == "pancreas")

# pick 5,000 most variable probes
top_var_probes <- readRDS(file = "./output/pancreas_top_variable_probes_training_set_01112023.rds")
top_var_probes <- top_var_probes[1:5000]


# filter betas
betas2 <- betas[top_var_probes, anno_outliers$arrayId]


# load models
rf_model <- readRDS(file = "./output/rf_model_default_01112023.rds")
xgb_model <- readRDS(file = "./output/xgb_model_default_01112023.rds")
nn_model <- load_model_hdf5(file = "./output/nn_model_01112023.hdf5")

# get model predictions
rf_pred_class2 <- predict(rf_model, newdata = t(betas2))
rf_pred_scores2 <- predict(rf_model, newdata = t(betas2), type = "prob")
rf_pred_scores_max2 <- apply(rf_pred_scores2, 1, max)
rf_pred_class2 <- apply(rf_pred_scores2, 1, function(x) colnames(rf_pred_scores2)[which.max(x)]) %>% as.factor

# xgb predictions
xgb_pred_class2 <- predict(xgb_model, newdata = t(betas2))
xgb_pred_scores2 <- predict(xgb_model, newdata = t(betas2), type = "prob")
xgb_pred_scores_max2 <- apply(xgb_pred_scores2, 1, max)

# nn predictions
nn_pred_scores2 <- predict(object = nn_model, x = t(betas2))
rownames(nn_pred_scores2) <- rownames(t(betas2))
colnames(nn_pred_scores2) <- colnames(rf_pred_scores)
nn_pred_scores_max2 <- apply(nn_pred_scores2 , 1, max)
nn_pred_class2 <- apply(nn_pred_scores2, 1, function(x) colnames(nn_pred_scores2)[which.max(x)]) %>% as.factor

# add performance to annotation
anno_outliers <- anno_outliers %>% 
  mutate(pred_nn = nn_pred_class2, 
         pred_rf = rf_pred_class2, 
         pred_xgb = xgb_pred_class2, 
         pred_scores_nn = apply(nn_pred_scores2, 1, max),
         pred_scores_rf = apply(rf_pred_scores2, 1, max), 
         pred_scores_xgb = apply(xgb_pred_scores2, 1, max))



# split data
train_set <- list()
test_set <- list()

train_set$x <- t(betas[, anno_outliers$cohort == "train"])
train_set$y <- as.factor(anno_outliers$tumorType[anno_outliers$cohort == "train"])

test_set$x <- t(betas[, anno_outliers$cohort == "test"])
test_set$y <- as.factor(anno_outliers$tumorType[anno_outliers$cohort == "test"])

# confusion matrices
test_indices <- which(anno_outliers$cohort == "test", arr.ind = TRUE)
conf_mat <- list(rf_pred_class2, xgb_pred_class2, nn_pred_class2)
conf_mat <- lapply(conf_mat, function(x) x[test_indices])
conf_mat <- lapply(conf_mat, function(x) confusionMatrix(x, test_set$y))


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


#saveRDS(object = perf_acc, file = "./output/perf_acc_ex_outliers_01112023.rds")
