########################################
### FIGURE 3: Classifier development ###
########################################

### Figure A: count training, test and validation set by location --------------------------------------------

anno %>% 
  ggplot(aes(cohort)) +
  geom_bar(aes(fill = tumorType)) +
  labs(#title="Tumor Type By Tumor Type",
    x = "", y = "No. of cases", fill = "TumorType") +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_x_discrete(limits = c("train", "test"),
                   breaks = c("train", "test"),
                   labels = c("Training Set", "Test Set")) +
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
ggsave("Figure 3A_count Dataset by TuorType.pdf", path= "./plots/", dpi=500)



### Figure B: model performance by test set -------------------------------------------------------------------------------

perf_acc %>%  
  ggplot(aes(x = method, y = Accuracy, fill = method)) +
  geom_col(width = 0.7) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  geom_errorbar(aes(ymin = AccuracyLower, ymax = AccuracyUpper), width = 0.3) +
  scale_x_discrete(limits = c("rf", "xgb", "nn"),
                   breaks = c("rf", "xgb", "nn"),
                   labels = c("Random Forest", "XGBoost", "Neural Network")) +
  geom_text(aes(label=round(Accuracy*100, digits=1)), vjust=15, color="black", size=4) +
  theme_bw(base_size = 30) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "none",
        panel.grid = element_blank()) +
  labs(title="Algorithm Performance",x = NULL, y = "Accuracy (test cohort)") +
  ylim(0, 1)
ggsave("Figure3B_Algorithm Performance.pdf", path= "./plots/", dpi=500)



### Figure C: plot accuracy across different classes RF/NN/xgb -----------------------------------------------------------

perf_per_class_rf <- conf_mat[[1]]$byClass %>% 
  as_tibble %>% 
  mutate(class = colnames(rf_pred_scores))
colnames(perf_per_class_rf) <- make.names(colnames(perf_per_class_rf))

perf_per_class_rf <- perf_per_class_rf[,c(11:12)]
colnames(perf_per_class_rf) <- c("Accuracy_rf", "class")

perf_per_class_xgb <- conf_mat[[2]]$byClass %>% 
  as_tibble %>% 
  mutate(class = colnames(xgb_pred_scores))
colnames(perf_per_class_xgb) <- make.names(colnames(perf_per_class_xgb))

perf_per_class_xgb <- perf_per_class_xgb[,c(11:12)]
colnames(perf_per_class_xgb) <- c("Accuracy_xgb", "class")

perf_per_class_nn <- conf_mat[[3]]$byClass %>% 
  as_tibble %>% 
  mutate(class = colnames(nn_pred_scores))
colnames(perf_per_class_nn) <- make.names(colnames(perf_per_class_nn))

perf_per_class_nn <- perf_per_class_nn[,c(11:12)]
colnames(perf_per_class_nn) <- c("Accuracy_nn", "class")

perf_per_class2 <- merge(perf_per_class_rf, perf_per_class_xgb, by = "class")
perf_per_class2 <- merge(perf_per_class2, perf_per_class_nn, by = "class")

perf_per_class2$class <- factor(perf_per_class2$class, levels = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"))

perf_per_class2 %>%
  gather(key = algorithm, value = Value, factor(c("Accuracy_rf", "Accuracy_xgb","Accuracy_nn"), levels=c("Accuracy_rf", "Accuracy_xgb","Accuracy_nn"))) %>%
  ggplot(aes(class, Value, fill = algorithm)) + 
  geom_col(position = "dodge") +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", limits = c("Accuracy_rf", "Accuracy_xgb","Accuracy_nn")) +
  
  scale_x_discrete(limits = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"),
                   breaks = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"),
                   labels = c("PanNET", "PDAC", "ACC", "SPN", "NORMAL", "PB", "PanNEC")) +
  #scale_fill_color(limits = c("Accuracy_rf", "Accuracy_xgb","Accuracy_nn")) +
  labs(title="Algorithm Performance By Tumor Type",
       x = "", 
       y = "Accuracy (test cohort)") +
  theme_bw(base_size = 30) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank()) +
  ylim(0, 1)
ggsave("Figure3C_model_comparison_accuracies_all models.pdf", path= "./plots/", dpi=500)



### Figure D, E, F: cutoff vs. predictable / accuracy -----------------------------------------------------------------------------------------------------

performance_cutoff %>% 
  #filter(method == "RF") %>% 
  #filter(method == "XGB") %>% 
  #filter(method == "NN") %>% 
  pivot_longer(cols = accuracy:predictable, 
               names_to = "statistic") %>% 
  ggplot(aes(cutoff, value, color = statistic)) +
  geom_point(size = 5) +
  geom_line() +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe")+
  #scale_color_manual(values = branded_colors1) +
  theme_bw(base_size = 18) +
  theme(legend.title = element_blank(), 
        legend.position = "top") +
  facet_grid(cols = vars(method)) +
  labs(x = "Cutoff", y = "Accuracy/Predictable (%)")
ggsave("Figure 3DEF_cutoff_vs_predictable_accuracy.pdf", path= "./plots/")



### Figures G, H, I: Confusionmatrix table for class prediction ----------------------------------------------------------------------

# nn
nn_cm <- table(prediction = as.factor(nn_pred_class), actual = anno$tumorType) %>% 
  caret::confusionMatrix()

plt <- as.data.frame(nn_cm$table)
plt$prediction <- factor(plt$prediction, levels=rev(levels(plt$prediction)))
plt$actual <- factor(plt$actual, levels = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"))
plt$prediction <- factor(plt$prediction, levels = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"))

ggplot(plt, aes(prediction,actual, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq), color = "white") +
  scale_fill_gradient(low = "#132B43", high = "#88CCEEFF") +
  #paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  labs(x = "True Class ",y = "Predicted Class (NN)") +
  scale_x_discrete(limits = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"),
                   breaks = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"),
                   labels = c("PanNET", "PDAC", "ACC", "SPN", "NORMAL", "PB", "PanNEC")) +
  theme_bw(base_size = 16) +
  theme(axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks=element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank())
ggsave("Figure3G_ConfusionMatrix_nn.pdf", path= "./plots/")

## rf
rf_cm <- table(prediction = rf_pred_class, actual = anno$tumorType) %>% 
  caret::confusionMatrix()
plt <- as.data.frame(rf_cm$table)
plt$prediction <- factor(plt$prediction, levels=rev(levels(plt$prediction)))
plt$actual <- factor(plt$actual, levels = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"))
plt$prediction <- factor(plt$prediction, levels = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"))


ggplot(plt, aes(prediction,actual, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq), color = "white") +
  labs(x = "True Class ",y = "Predicted Class (RF)") +
  scale_x_discrete(limits = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"),
                   breaks = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"),
                   labels = c("PanNET", "PDAC", "ACC", "SPN", "NORMAL", "PB", "PanNEC")) +
  theme_bw(base_size = 16) +
  theme(axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks=element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank())
ggsave("Figure3H_ConfusionMatrix_rf.pdf", path= "./plots/")

## xgb
xgb_nn<-table(prediction = xgb_pred_class, actual = anno$tumorType) %>% 
  caret::confusionMatrix()
plt <- as.data.frame(xgb_nn$table)
plt$prediction <- factor(plt$prediction, levels=rev(levels(plt$prediction)))
plt$actual <- factor(plt$actual, levels = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"))
plt$prediction <- factor(plt$prediction, levels = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"))


ggplot(plt, aes(prediction,actual, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq), color = "white") +
  labs(x = "True Class",y = "Predicted Class (XGB)") +
  scale_x_discrete(limits = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"),
                   breaks = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"),
                   labels = c("PanNET", "PDAC", "ACC", "SPN", "NORMAL", "PB", "PanNEC")) +
  theme_bw(base_size = 16) +
  theme(axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks=element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank())
ggsave("Figure3I_ConfusionMatrix_xgb.pdf", path= "./plots/")



### Figure J: correct and incorrect classification -------------------------------------------------------------------------
anno %>% 
  mutate(correct = ifelse(tumorType == pred_rf, "correct", "incorrect")) %>% 
  filter(cohort == "test") %>%
  ggplot(aes(umap_x, umap_y, col = tumorType, shape = correct)) +
  geom_point(size = 3) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe",
                                     limits = c("ACC", "NORM", "PanNEC","PanNET", "PB", "PDAC", "SPN"),
                                     breaks = c("ACC", "NORM", "PanNEC","PanNET", "PB", "PDAC", "SPN"),
                                     labels = c("ACC", "NORM", "PanNEC","PanNET", "PB", "PDAC", "SPN")) +
  labs(#title="UMAP Clustering",
    x = "UMAP 1", 
    y = "UMAP 2",
    col = "Tumor Type") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 12),
    axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=0, size = 10),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
    axis.text.x = element_text(face = "bold", size = 9),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.ticks.length.y =unit(.05, "cm"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, 'cm'),
    panel.grid = element_blank()) +
  guides(col = guide_legend(nrow = 1))
ggsave("Figure3J_correct classification TEST COHORT.pdf", path= "./plots/", dpi = 500, height = 6, width = 5.5, units = "in")



##### EXTRA ######

### Figure G: UMAP entire cohort ------------------------------------------------------------------------------------------
anno %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType)) + 
  geom_point(shape=19, size = 3) +
  scale_color_manual(limits = c("ACC", "NORM", "PanNEC","PanNET", "PB", "PDAC", "SPN"),
                     breaks = c("ACC", "NORM", "PanNEC","PanNET", "PB", "PDAC", "SPN"),
                     values = branded_colors1, 
                     labels = c("ACC", "NORM", "PanNEC","PanNET", "PB", "PDAC", "SPN")) +
  labs(title="UMAP Clustering",
       x = "UMAP 1", 
       y = "UMAP 2",
       col = "Tumor Type") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 8),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = c(.89, .2),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank())
ggsave("Figure 2_UMAP-all07042023.pdf", path= "./plots/")

### UMAP by cohort
anno %>% 
  #mutate(correct = ifelse(tumorType == pred_rf, "correct", "incorrect")) %>% 
  ggplot(aes(umap_x, umap_y, color = tumorType, shape = cohort)) +
  geom_point(size = 3) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe") +
  labs(#title="UMAP Clustering",
    x = "UMAP 1", 
    y = "UMAP 2",
    col = "Tumor Type") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 12),
    axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=0, size = 10),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
    axis.text.x = element_text(face = "bold", size = 9),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.ticks.length.y =unit(.05, "cm"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, 'cm'),
    panel.grid = element_blank()) +
  guides(col = guide_legend(nrow = 1))
ggsave("Figure3_UMAP by cohort.pdf", path= "./plots/", dpi = 500, height = 6, width = 5.5, units = "in")


### Figure 3_ probability score distribution per tumortype --------------------------------------------------------------------------------

anno %>% 
  select(tumorType, starts_with("pred")) %>%
  pivot_longer(cols = starts_with("pred_scores"), names_to = "method", values_to = "score") %>% 
  mutate(method = str_extract(string = method, pattern = "[^_]*$")) %>% 
  ggplot(aes(score, tumorType)) +
  geom_jitter(aes(col= method), height = 0.25, size = 2) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 30) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank())
ggsave("Figure 3_extra_distribution probability score per tumortype.pdf", path= "./plots/", dpi = 500)

### Figure 3_ probability score distribution by source --------------------------------------------------------------------------------

anno %>% 
  select(source, starts_with("pred")) %>%
  pivot_longer(cols = starts_with("pred_scores"), names_to = "method", values_to = "score") %>% 
  mutate(method = str_extract(string = method, pattern = "[^_]*$")) %>% 
  ggplot(aes(score, source)) +
  geom_jitter(aes(col= method), height = 0.25, size = 2) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 30) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank())
ggsave("Figure 3_extra_distribution probability score by source.pdf", path= "./plots/", dpi = 500)
