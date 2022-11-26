# Christoph Geisenberger
# github: @cgeisenberger
# last edited 01/09/2022


### assess implementation of cutoff for classifier scores ----------------------

# get true classes
cutoff_true_classes <- test_anno$label

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
                                 label_pred = pred_xgb_class, 
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
  scale_color_manual(values = branded_colors) +
  theme_bw(base_size = 30) +
  theme(legend.position = "none") +
  facet_grid(cols = vars(method)) +
  labs(x = "Cutoff", y = "Accuracy/Predictable (%)")

test_anno %>% 
  select(source, platform, label, location, starts_with("pred")) %>% 
  mutate(location = ifelse(location == "primary", "primary", "other")) %>% 
  pivot_longer(cols = starts_with("pred_scores"), names_to = "method", values_to = "score") %>% 
  mutate(method = str_extract(string = method, pattern = "[^_]*$")) %>% 
  ggplot(aes(method, score, col = method)) +
  geom_jitter(size = 3, width = 0.1) +
  scale_color_manual(values = branded_colors2) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none") +
  labs(x = "Score", y = NULL)

test_anno %>% 
  select(source, platform, label, location, starts_with("pred")) %>% 
  mutate(location = ifelse(location == "primary", "primary", "other")) %>% 
  pivot_longer(cols = starts_with("pred_scores"), names_to = "method", values_to = "score") %>% 
  mutate(method = str_extract(string = method, pattern = "[^_]*$")) %>% 
  filter(method == "nn") %>% 
  ggplot(aes(source, score, col = label)) +
  geom_jitter(size = 3, width = 0.1) +
  scale_color_manual(values = branded_colors2) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top") +
  labs(x = "Score", y = NULL)
