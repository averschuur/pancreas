

accuracy_nn_99 <- test_anno %>%
  filter(pred_scores_nn >= 0.99)

156/201

accuracy_rf_99 <- test_anno %>%
  filter(pred_scores_rf >= 0.99)

0/201

accuracy_xgb_99 <- test_anno %>%
  filter(pred_scores_xgb >= 0.99)

3/201

accuracy_nn_95 <- test_anno %>%
  filter(pred_scores_nn >= 0.95)

175/201

accuracy_rf_95 <- test_anno %>%
  filter(pred_scores_rf >= 0.95)

3/201

accuracy_xgb_95 <- test_anno %>%
  filter(pred_scores_xgb >= 0.95)

65/201


accuracy_nn_99 %>% 
  mutate(nn_corr = as.integer(tumorType == pred_nn)) %>% 
  group_by(tumorType) %>% 
  summarise(sum(nn_corr) / n())

accuracy_rf_99 %>% 
  mutate(rf_corr = as.integer(tumorType == pred_rf)) %>% 
  group_by(tumorType) %>% 
  summarise(sum(rf_corr) / n())

accuracy_xgb_99 %>% 
  mutate(xgb_corr = as.integer(tumorType == pred_xgb)) %>% 
  group_by(tumorType) %>% 
  summarise(sum(xgb_corr) / n())

