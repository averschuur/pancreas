###################################
### FIGURE 4: Outlier detection ###
###################################


# Figure 4A: Tumors in TCGA dataset --------------------------------------------
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
ggsave("Figure 4A_count TCGA Dataset_06112023.pdf", path= "./plots/", dpi=500)


### Figure 4B: plot confusion matrix RF ------------------------------------------
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
ggsave("Figure 4B_TCGA prediction by RF_06112023.pdf", path= "./plots/", dpi=500)


### Figure 4C: RF score distribution TCGA samples -------------------------
rf_data %>%
  ggplot(aes(class_char, winning_score, fill = class_char)) +
  geom_boxplot(alpha = 0.5) +
  labs(#title="Tumor Type By Source",
    x = "", y = "Random Forest Score Distribution", fill = "Location") +
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
ggsave("Figure 4C_RF score distribution TCGA samples_06112023.pdf", path= "./plots/", dpi=500)

### Figure 4D: ROC -------------------------
ggroc(rf_data_roc) +
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
ggsave("Figure 4D_RF ROC_no logistic regression_06112023.pdf", width = 8, height = 3, units = "in", path= "./plots/", dpi=500)


# Figure 4E: Outlier probability after logistic regression --------------------------------------------------

rf_data %>% 
  filter(dataset == "test") %>%
  ggplot(aes(class_char, od_score, fill = class_char)) +
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
ggsave("Figure 4E_score distribution logistic regresiion_06112023.pdf", path= "./plots/", dpi=500)

### Figure 4F: plot ROC after logistic regression -------------------------
ggroc(od_prob_roc) +
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
ggsave("Figure 4F_RF ROC_after logistic regression_06112023.pdf", width = 8, height = 3, units = "in", path= "./plots/", dpi=500)



# Figure 4GH: Plot dropout and pass------------------------------------------------
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
ggsave("Figure 4GH_dropout and pass_lr_06112023.pdf", path= "./plots/", dpi=500)
