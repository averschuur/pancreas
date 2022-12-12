##############################
### FIGURE 4: Case Reports ###
##############################


test <- pancreas_scores_rf %>% 
  as_tibble(rownames = "sample_id") %>% 
  pivot_longer(cols = -sample_id, names_to = "tumor_type", values_to = "rf_score") %>% 
  mutate(rf_score = as.numeric(rf_score))

i = 4
test %>% 
  slice_head(n = ncol(pancreas_scores_rf) * i) %>% 
  ggplot(aes(tumor_type, rf_score, fill = tumor_type)) +
  geom_col() +
  theme_bw(base_size = 18) +
  theme(panel.border = element_blank(), 
        legend.position = "none") +
  labs(x = NULL, y = NULL) +
  coord_polar() + 
  facet_wrap(~ sample_id)