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
ggsave("Figure 4A_count TCGA Dataset.pdf", path= "./plots/", dpi=500)


### Figure 4B: plot confusion matrix ------------------------------------------
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
ggsave("Figure 4B_tcga prediction by rf.pdf", path= "./plots/", dpi=500)


### Figure 4C-E: RF score distributions for the three classes from B for outliers vs. pancreatic cancers-------------------------
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
ggsave("Figure 4CDE_RF score distributions for the three classes.pdf", width = 8, height = 3, units = "in", path= "./plots/", dpi=500)

### stats ACC
rf_data_ACC <- rf_data %>%
  filter(winning_class %in% "ACC")

rf_data_ACC %>%
  group_by(class_char) %>% 
  summarise(n = n())
#1 outlier       76
#2 pancreas      48
  
t.test(rf_data_ACC$winning_score ~ rf_data_ACC$class_char,  paired = F)
#Welch Two Sample t-test

#data:  rf_data_ACC$winning_score by rf_data_ACC$class_char
#t = -17.655, df = 49.826, p-value < 0.00000000000000022
#alternative hypothesis: true difference in means between group outlier and group pancreas is not equal to 0
#95 percent confidence interval:
#  -0.5749651 -0.4574954
#sample estimates:
#  mean in group outlier mean in group pancreas 
#0.2953947              0.8116250 

t.test(rf_data_ACC$winning_score ~ rf_data_ACC$class_char,  paired = F, var.equal=TRUE)
#Two Sample t-test

#data:  rf_data_ACC$winning_score by rf_data_ACC$class_char
#t = -21.791, df = 122, p-value < 0.00000000000000022
#alternative hypothesis: true difference in means between group outlier and group pancreas is not equal to 0
#95 percent confidence interval:
#  -0.5631273 -0.4693332
#sample estimates:
#  mean in group outlier mean in group pancreas 
#0.2953947              0.8116250 

### stats pNET
rf_data_pNET <-rf_data %>%
  filter(winning_class %in% "PanNET")

t.test(rf_data_pNET$winning_score ~ rf_data_pNET$class_char,  paired = F, var.equal=TRUE)
#Two Sample t-test

#data:  rf_data_pNET$winning_score by rf_data_pNET$class_char
#t = -29.068, df = 394, p-value < 0.00000000000000022
#alternative hypothesis: true difference in means between group outlier and group pancreas is not equal to 0
#95 percent confidence interval:
#  -0.4917732 -0.4294649
#sample estimates:
#  mean in group outlier mean in group pancreas 
#0.4131111              0.8737302

### stats PDAC
rf_data_PDAC <-rf_data %>%
  filter(winning_class %in% "PanNET")

t.test(rf_data_PDAC$winning_score ~ rf_data_PDAC$class_char,  paired = F, var.equal=TRUE)
#Two Sample t-test

#data:  rf_data_PDAC$winning_score by rf_data_PDAC$class_char
#t = -29.068, df = 394, p-value < 0.00000000000000022
#alternative hypothesis: true difference in means between group outlier and group pancreas is not equal to 0
#95 percent confidence interval:
#  -0.4917732 -0.4294649
#sample estimates:
#  mean in group outlier mean in group pancreas 
#0.4131111              0.8737302 

# Figure 4F: Outlier probability --------------------------------------------------
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
ggsave("Figure 4F_score distribution logistic regresiion.pdf", path= "./plots/", dpi=500)


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
ggsave("Figure 4GH_dropout and pass_lr.pdf", path= "./plots/", dpi=500)


######### EXTRA FIGURES #########

# plot ROC curve
od_prob_roc <- rf_data %>%
  pROC::roc(class_int, od_prob)
pROC::plot.roc(od_prob_roc)
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
ggsave("Figure 4_ROC.pdf", path= "./plots/", dpi=500)

                 