########################################
### FIGURE 2: Classifier development ###
########################################

library(gridExtra)

# Figure 2 umap
anno %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType)) + 
  geom_point(shape=19, size = 3, alpha = 0.7) +
  scale_color_manual(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                     breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                     values = branded_colors, 
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
nn_cm <- table(prediction = as.factor(pred_nn_class), actual = test_anno$tumorType) %>% 
  caret::confusionMatrix()

cvms::plot_confusion_matrix(
  nn_cm$`Confusion Matrix`[[1]],
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

ggplot(plt, aes(prediction,actual, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq), color = "white") +
  scale_color_viridis(option = "D") +
  labs(x = "True Class ",y = "Predicted Class (NN)") +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL"),
                   position = "top") +
  scale_y_discrete(limits = c("normal", "PB", "SPN", "PanNEC", "PanNET", "ACC",  "PDAC"),
                   breaks = c("normal", "PB", "SPN", "PanNEC", "PanNET", "ACC",  "PDAC"),
                   labels = c("NORMAL", "PB", "SPN", "PanNEC", "PanNET", "ACC",  "PDAC")) +
  theme_bw(base_size = 16) +
  theme(axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks=element_blank(),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank())
ggsave("Figure2_ConfusionMatrix_nn.png", path= "./output/")

## rf
rf_cm <- table(prediction = pred_rf_class, actual = test_anno$tumorType) %>% 
  caret::confusionMatrix()
plt <- as.data.frame(rf_cm$table)
plt$prediction <- factor(plt$prediction, levels=rev(levels(plt$prediction)))

ggplot(plt, aes(prediction,actual, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq), color = "white") +
  scale_color_viridis(option = "D") +
  labs(x = "True Class ",y = "Predicted Class (RF)") +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL"),
                   position = "top") +
  scale_y_discrete(limits = c("normal", "PB", "SPN", "PanNEC", "PanNET", "ACC",  "PDAC"),
                   breaks = c("normal", "PB", "SPN", "PanNEC", "PanNET", "ACC",  "PDAC"),
                   labels = c("NORMAL", "PB", "SPN", "PanNEC", "PanNET", "ACC",  "PDAC")) +
  theme_bw(base_size = 16) +
  theme(axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks=element_blank(),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank())
ggsave("Figure2_ConfusionMatrix_rf.png", path= "./output/")

## xgb
xgb_nn<-table(prediction = pred_xgb_class, actual = test_anno$tumorType) %>% 
  caret::confusionMatrix()
plt <- as.data.frame(xgb_nn$table)
plt$prediction <- factor(plt$prediction, levels=rev(levels(plt$prediction)))

ggplot(plt, aes(prediction,actual, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq), color = "white") +
  scale_color_viridis(option = "D") +
  labs(x = "True Class",y = "Predicted Class (XGB)") +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL"),
                   position = "top") +
  scale_y_discrete(limits = c("normal", "PB", "SPN", "PanNEC", "PanNET", "ACC",  "PDAC"),
                   breaks = c("normal", "PB", "SPN", "PanNEC", "PanNET", "ACC",  "PDAC"),
                   labels = c("NORMAL", "PB", "SPN", "PanNEC", "PanNET", "ACC",  "PDAC")) +
  theme_bw(base_size = 16) +
  theme(axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks=element_blank(),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank())
ggsave("Figure2_ConfusionMatrix_xgb.png", path= "./output/")




## create horizontal bar as in Leitheiser Figure 2

test_anno %>%
  ggplot(aes(fill=tumorType, y="count")) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = branded_colors) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank())

test_anno %>%
  gather(key = algorithm, value = Value, c("pred_nn", "pred_rf", "pred_xgb")) %>%
  ggplot(aes(fill = Value, algorithm)) +
  geom_bar(position = "fill", width = 1) +
  scale_fill_manual(values = branded_colors) +
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
