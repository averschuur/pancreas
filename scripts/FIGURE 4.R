###################################
### FIGURE 4: Outlier detection ###
###################################

# number of samples per dataset
anno_tcga <- anno_tcga %>%
  add_column(add_column = "TCGA")
anno_pancreas <- anno[,c(1,5)]
anno_pancreas <- anno_pancreas %>%
  add_column(add_column = "pancreas")

anno_total <- full_join(anno_tcga, anno_pancreas)

anno_total %>% 
  ggplot(aes(add_column)) +
  geom_bar(aes(fill = add_column)) +
  labs(#title="Tumor Type By Source",
    x = "", y = "No. of cases", fill = "Location") +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, vjust=0, hjust=0.5),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.4, 'cm'),
    panel.grid = element_blank())
ggsave("Figure 3_count Dataset.pdf", path= "./plots/", dpi=500)

## Tumors in TCGA dataset
colors <- c("#24693DFF", "#2E7542FF", "#368046FF", "#3E8A4AFF", "#46944EFF", 
            "#529E54FF", "#5EA85BFF", "#69B262FF", "#75BC69FF", "#8BC482FF",
            "#A0CB9CFF", "#B3D2B5FF", "#C7D9CFFF", "#B9CFD2FF", "#A6C3D1FF", 
            "#92B8D0FF", "#7DACCFFF", "#70A1C8FF", "#6697C0FF", "#5C8DB7FF",
            "#5283AEFF", "#4D7EAAFF", "#4979A5FF", "#4474A0FF", "#406F9BFF", 
            "#3C6A96FF", "#376591FF", "#33608DFF", "#2F5C88FF", "#2A5783FF")

anno_tcga %>% 
  ggplot(aes(tumorType)) +
  geom_bar(aes(fill = tumorType)) +
  coord_flip() +
  scale_fill_manual(values = colors) +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, vjust=0, hjust=0.5),
    axis.ticks = element_blank(),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("Figure 3_count TCGA Dataset.pdf", path= "./plots/", dpi=500)

### plot confusion matrix
rf_scores_tcga <- rf_scores_tcga %>%
  mutate(winning_class = rf_class_tcga, 
         winning_score = rf_scores_tcga_max, 
         tumorType = anno_tcga$tumorType,
         class_int = rep(0, length(rf_class_tcga)),
         class_char = ifelse(class_int == 0, "outlier", "pancreas"))

rf_tcga_plot <- rf_scores_tcga %>% 
  select(tumorType, winning_class) %>% 
  group_by(tumorType, winning_class) %>%
  count() 

rf_tcga_plot %>%
  ggplot(aes(winning_class, tumorType)) +
  geom_tile(aes(fill=n)) +
  geom_text(aes(label=n), color = "black") +
  scale_fill_gradient(low = "#C7D9CFFF", high = "#3A8548FF") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, vjust=0, hjust=0.5),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.4, 'cm'),
    panel.grid = element_blank()) +
  labs(x = NULL)
ggsave("Figure3_tcga prediction by rf.pdf", path= "./plots/", dpi=500)

### RF score distributions for the three classes from B for outliers vs. pancreatic cancers
rf_data %>% 
  filter(winning_class %in% c("ACC", "PanNET", "PDAC")) %>% 
  ggplot(aes(class_char, winning_score, fill = winning_class)) +
  geom_boxplot(alpha = 0.5) +
  facet_wrap(facets = vars(winning_class)) +
  labs(x = NULL, y = "Random Forest Score") +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", labels = c("Non-pancreatic Tumor", "Pancreatic Tumor")) +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, vjust=0, hjust=0.5),
    axis.ticks = element_blank(),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("Figure 4_RF score distributions for the three classes.pdf", width = 8, height = 3, units = "in", path= "./plots/", dpi=500)

### Outlier probability
rf_data %>% 
  ggplot(aes(class_char, od_prob, fill = class_char)) +
  geom_boxplot(alpha = 0.5) +
  labs(#title="Tumor Type By Source",
    x = "", y = "Outlier probability", fill = "Location") +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", labels = c("Non-pancreatic tumor", "Pancreatic Tumor")) +
  scale_x_discrete(limits = c("outlier", "pancreas"),
                   breaks = c("outlier", "pancreas"),
                   labels = c("Non-pancreatic tumor", "Pancreatic tumor")) +
  #geom_hline(yintercept = 0.5, lty = 2, col = "steelblue") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, vjust=0, hjust=0.5),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.4, 'cm'),
    panel.grid = element_blank())
ggsave("Figure 3_score distribution.pdf", path= "./plots/", dpi=500)

### Plot dropout and pass
cutoff_data %>% 
  pivot_longer(cols = -cutoff) %>% 
  filter(cutoff > 0) %>% 
  ggplot(aes(cutoff, value, col = name)) + 
  geom_line(lwd = 2) +
  geom_vline(xintercept = 0.5, col = "steelblue", lty = 2) +
  facet_wrap(facets = vars(name), scales = "free", ncol = 1) +
  theme_bw(base_size = 18) +
  labs(x = "Cutoff (outlier prob)", y = "Samples (%)") +theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, vjust=0, hjust=0.5),
    axis.ticks = element_blank(),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("Figure 4_dropout and pass.pdf", path= "./plots/", dpi=500)


# plot ROC curve
od_prob_roc <- rf_data %>%
  pROC::roc(class_int, od_prob)
plot.roc(od_prob_roc)
od_prob_roc

ggroc(od_prob_roc, color = 'steelblue', size =1) +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, vjust=0, hjust=0.5),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.4, 'cm'),
    panel.grid = element_blank())
ggsave("Figure 3_ROC.pdf", path= "./plots/", dpi=500)



anno %>%
  mutate(online = case_when(source == "Benhamida" ~ "online",
                            source == "Chan" ~ "online",
                            source == "DiDomenico" ~ "online",
                            source == "Endo" ~ "online",
                            source == "Jakel" ~ "online",
                            source == "Selenica" ~ "online",
                            source == "yachida" ~ "online",
                            source == "UMCU" ~ "local")) %>%
  ggplot(aes(online)) +
  geom_bar(aes(fill = online)) +
  labs(#title="Tumor Type By Source",
    x = "", y = "No. of cases", fill = "Location") +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, vjust=0, hjust=0.5),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.4, 'cm'),
    panel.grid = element_blank())
ggsave("Figure 2_count online Dataset.pdf", path= "./plots/", dpi=500)

correct = ifelse(tumorType == pred_rf, "correct", "incorrect"
                 