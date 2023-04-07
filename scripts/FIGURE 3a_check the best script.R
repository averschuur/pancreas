########################################
### FIGURE 3: Classifier development ###
########################################

library(gridExtra)

### Figure2; count training, test and validation set by location --------------------------------------------

# plot
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
ggsave("Figure 2A_count Dataset by TuorType.pdf", path= "./plots/", dpi=500)


### Figure2; model performance by test set --------------------------------------------

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
ggsave("Figure3_Algorithm Performance.pdf", path= "./plots/", dpi=500)

# plot accuracy across different classes RF
anno$tumorType <- factor(anno$tumorType, levels = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"))
anno$source <- factor(anno$source, levels = c("Benhamida", "Chan", "DiDomenico", "Endo", "Jakel", "Selenica", "yachida", "UMCU"))

perf_per_class %>% 
  ggplot(aes(class, Balanced.Accuracy, fill = class)) +
  geom_col(width = 0.7) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_x_discrete(limits = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"),
                   breaks = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"),
                   labels = c("PanNET", "PDAC", "ACC", "SPN", "NORMAL", "PB", "PanNEC")) +
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
ggsave("Figure3_model_comparison_accuracies_randforest.pdf", path= "./plots/", dpi=500)



# plot accuracy across different classes RF/NN/xgb
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
ggsave("Figure3_model_comparison_accuracies_all models.pdf", path= "./plots/", dpi=500)



# Figure 2 umap
anno %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType)) + 
  geom_point(shape=19, size = 3, alpha = 0.7) +
  scale_color_manual(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                     breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                     values = branded_colors1, 
                     labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
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
ggsave("Figure 2_UMAP-all.png", path= "./output/")


# Figure 2_ correct and incorrect classification
anno %>% 
  mutate(correct = ifelse(tumorType == pred_rf, "correct", "incorrect")) %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType, shape = correct)) +
  geom_point(size = 3) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe",
                                     limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                                     breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                                     labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
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
ggsave("Figure3_correct classification.pdf", path= "./plots/", dpi = 500, height = 6, width = 5.5, units = "in")

# Figure 2_ scores per method per tumorType

anno %>% 
  select(tumorType, starts_with("pred")) %>%
  pivot_longer(cols = starts_with("pred_scores"), names_to = "method", values_to = "score") %>% 
  mutate(method = str_extract(string = method, pattern = "[^_]*$")) %>% 
  ggplot(aes(score, tumorType)) +
  geom_jitter(aes(col= method), height = 0.25, size = 1) +
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
ggsave("Figure 2_score per tumortype.pdf", path= "./plots/", dpi = 500)



# Figure 2_UMAP Test and training cohort
anno %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType)) + 
  geom_point(aes(shape = split), size = 3, alpha = 0.7) +
  scale_color_manual(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                     breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                     values = branded_colors, 
                     labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  scale_shape_manual(limits = c("test", "train"),
                     breaks = c("test", "train"),
                     values = c(17, 19),
                     labels = c("Test Cohort", "Training Cohort")) +
  labs(title="UMAP Clustering",
       x = "UMAP 1", 
       y = "UMAP 2",
       col = "Tumor Type",
       shape = "Dataset") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 8),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank())
ggsave("Figure2_UMAP Test and training cohort.png", path= "./output/")

#Figure 2 t-SNE Test and Training
anno %>% 
  ggplot(aes(tsne_x, tsne_y, col = tumorType)) + 
  geom_point(aes(shape = split), size = 3, alpha = 0.7) +
  scale_color_manual(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                     breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                     values = branded_colors, 
                     labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  scale_shape_manual(limits = c("test", "train"),
                     breaks = c("test", "train"),
                     values = c(17, 19),
                     labels = c("Test Cohort", "Training Cohort")) +
  labs(title="t-SNE Clustering",
       x = "t-SNE 1", 
       y = "t-SNE 2",
       col = "Tumor Type",
       shape = "Dataset") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 8),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank())
ggsave("Figure2_t-SNE Test and training cohort.png", path= "./output/")

# Figure 2 Accuracy Test Cohort
df <-data.frame((c("nn", "rf", "xgb")), (c("93.9", "94.7", "92.8")))

test_anno %>% 
  mutate(nn = as.integer(tumorType == pred_nn), 
         rf = as.integer(tumorType == pred_rf), 
         xgb = as.integer(tumorType == pred_xgb)) %>% 
  select(nn, rf, xgb) %>% 
  pivot_longer(cols = everything()) %>% 
  group_by(name) %>% 
  summarise(accuracy = sum(value)/n()) %>% 
  ggplot(aes(name, accuracy, fill = name)) +
  geom_col(width = 0.7) +
  geom_text(aes(label=df[,2]), vjust=1.6, color="white", size=3.5)+
  scale_fill_manual(values = branded_colors) +
  scale_x_discrete(limits = c("nn", "rf", "xgb"),
                   breaks = c("nn", "rf", "xgb"),
                   labels = c("Neural Network", "Random Forest", " XGBoost")) +
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
ggsave("Figure2_Algorithm Performance.png", path= "./output/")

# algorithm performance by tumorType
x <- test_anno %>% 
  mutate(nn_corr = as.integer(tumorType == pred_nn),
         rf_corr = as.integer(tumorType == pred_rf),
         xgb_corr = as.integer(tumorType == pred_xgb)) %>% 
  group_by(tumorType) %>% 
  summarise(sum(nn_corr) / n(), sum(rf_corr) / n(), sum(xgb_corr) / n())

colnames(x) <- c("tumorType", "nn", "rf", "xgb"))


x %>% 
  gather(key = algorithm, value = Value, c("nn", "rf", "xgb")) %>%
  ggplot(aes(tumorType, Value, fill = algorithm)) + 
  geom_col(position = "dodge") +
  geom_text(aes(label=format(round(Value*100, 1), nsmall = 1)), hjust=1.6, position = position_dodge(0.9), angle = 90, color="white", size=2.2)+
  scale_fill_manual(values = branded_colors2) +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  labs(title="Algorithm Performance By Tumor Type",
       x = "", 
       y = "Accuracy") +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, angle = 90, vjust=0.1, hjust=1),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank())
ggsave("Figure2_Algorithm PerformanceByTumorType.png", path= "./output/")


##################################################
### Confusionmatrix table for class prediction ###
##################################################
library(cvms)
library(broom)    
library(tibble)   
library(ggimage)   
library(rsvg) 

# nn
nn_cm <- table(prediction = as.factor(nn_pred_class), actual = anno$tumorType) %>% 
  caret::confusionMatrix()

cvms::plot_confusion_matrix(
  conf_mat[[1]],
  add_sums = TRUE,
  sums_settings = sum_tile_settings(
    palette = "Oranges",
    label = "Total",
    tc_tile_border_color = "black"
  )
)

# cvms does not work, alternative
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
ggsave("Figure3_ConfusionMatrix_nn.pdf", path= "./plots/")

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
ggsave("Figure2_ConfusionMatrix_rf.pdf", path= "./plots/")

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
ggsave("Figure3_ConfusionMatrix_xgb.pdf", path= "./plots/")




## create horizontal bar as in Leitheiser Figure 2

anno %>%
  ggplot(aes(fill=tumorType, y="count")) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = branded_colors1) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank())

anno %>%
  gather(key = algorithm, value = Value, c("pred_nn", "pred_rf", "pred_xgb")) %>%
  ggplot(aes(fill = Value, algorithm)) +
  geom_bar(position = "fill", width = 1) +
  scale_fill_manual(values = branded_colors1) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank())+
  coord_flip()


test_anno <- test_anno %>%
  mutate(Entity = case_when(location == "acc normal" ~ "normal",
                            location == "pancreas" ~ "normal",
                            location == "bone metastasis" ~ "metastasis",
                            location == "metastasis" ~ "metastasis",
                            location == "liver metastasis" ~ "metastasis",
                            location == "lymph node metastasis" ~ "metastasis",
                            location == "primary" ~ "primary")) 
test_anno %>%
  group_by(tumorType, Entity) %>%
  summarise(count = n()) %>%
  ggplot(aes(fill="Entity")) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = branded_colors) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank())

"pred_nn", "pred_rf", "pred_xgb"

gather(key = algorithm, value = Value, c("pred_nn", "pred_rf", "pred_xgb")) %>%
  ggplot(aes(tumorType, Value, fill = algorithm)) +
  
  
  
  
  
  # heatmap of correlation matrix
  heat <- pheatmap::pheatmap(dr_input, 
                             labels_row = anno$tumorType,
                             labels_col = "")



save_pheatmap_pdf <- function(x, filename, width=30, height=30) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(heat, "T:/pathologie/PRL/Groep-Brosens/2. Anna Vera/2. ACN-SPN-NET/Pancreas-ID/output/heatmap_whole cohort.pdf")
