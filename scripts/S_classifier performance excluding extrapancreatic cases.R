

outliers <- rf_data %>%
       filter(class_char == "pancreas")

outliers <- left_join(anno, outliers)

outlies <- outliers[c(1:277,279:326),] #exclude ACC2_M

anno_ex_outliers <- outliers %>%
  filter(od_class == "pancreas")

betas2 <- betas[top_var_probes, anno_ex_outliers$arrayId]


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
anno_ex_outliers <- anno_ex_outliers %>% 
  mutate(pred_nn = nn_pred_class2, 
         pred_rf = rf_pred_class2, 
         pred_xgb = xgb_pred_class2, 
         pred_scores_nn = apply(nn_pred_scores2, 1, max),
         pred_scores_rf = apply(rf_pred_scores2, 1, max), 
         pred_scores_xgb = apply(xgb_pred_scores2, 1, max))



# split data
train_set <- list()
test_set <- list()

train_set$x <- t(betas[, anno_ex_outliers$cohort == "train"])
train_set$y <- as.factor(anno_ex_outliers$tumorType[anno_ex_outliers$cohort == "train"])

test_set$x <- t(betas[, anno_ex_outliers$cohort == "test"])
test_set$y <- as.factor(anno_ex_outliers$tumorType[anno_ex_outliers$cohort == "test"])

# confusion matrices
test_indices <- which(anno_ex_outliers$cohort == "test", arr.ind = TRUE)
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

#saveRDS(object = perf_acc, file = "./output/perf_acc_ex_outliers.rds")
