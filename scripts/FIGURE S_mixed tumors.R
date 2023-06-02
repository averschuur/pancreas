##############################
### FIGURE S: Mixed tumors ###
##############################

### Plot UMAP
anno %>% 
  ggplot(aes(umap_x, umap_y)) + 
  geom_point(aes(color = tumorType), shape=19, size = 4) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe",
                                     limits = c("ACC", "Mixed", "NORM", "PanNEC", "PanNET", "PB", "PDAC", "SPN"),
                                     breaks = c("ACC", "Mixed", "NORM", "PanNEC", "PanNET", "PB", "PDAC", "SPN"),
                                     labels = c("ACC", "MIXED", "NORM", "PanNEC", "PanNET", "PB", "PDAC", "SPN")) +
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
    #legend.margin=margin(0,0,0,0),
    #legend.box = "horizontal",
    #legend.box.just = "left",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, 'cm'),
    panel.grid = element_blank()) +
  guides(col = guide_legend(nrow = 1))
ggsave("Figure S_Mixed_UMAP-all.pdf", path= "./Plots/", dpi=500)

### RF Score entrophy mixed vs. non-mixed tumors
anno %>% 
  ggplot(aes(mixed, rf_entropy, fill = mixed)) +
  geom_boxplot() +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank()) +
  labs(x = NULL, y = "RF score entropy")
ggsave("Figure S_RF score entrophy.pdf", path= "./plots/", dpi=500)

## estimate scores mixed vs. non-mixed tumors
anno %>% 
  ggplot(aes(mixed, estimate, fill = mixed)) +
  geom_boxplot() +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank()) +
  labs(x = NULL, y = "Tumor Purity")
ggsave("Figure S_estimate scores mixed vs. non-mixed tumors.pdf", path= "./plots/", dpi=500)

### mean probability scores per tumor type
anno %>% 
  select(tumorType, starts_with("scores")) %>% 
  pivot_longer(cols = starts_with("scores")) %>% 
  mutate(name = str_replace_all(name, "scores_", "")) %>% 
  mutate(name = str_replace_all(name, "Pan", "")) %>% 
  group_by(tumorType, name) %>% 
  summarise(mean_score = mean(value)) %>% 
  ungroup %>% 
  ggplot(aes(name, mean_score, fill = name)) +
  geom_col() +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", breaks = c("ACC", "Mixed", "NORM", "PanNEC", "PanNET", "PB", "PDAC", "SPN")) +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank()) +
  labs(y = "Random Forest Score (mean)", fill = "Tumor Type") +
  facet_wrap(facets = vars(tumorType))
ggsave("Figure S_Probability scores per tumor type.pdf", path= "./plots/", dpi=500)

## extra figure
anno %>% 
  filter(tumorType %in% c("PB", "ACC", "Mixed")) %>% 
  ggplot(aes(scores_ACC, scores_PB)) +
  geom_point(aes(size = mixed, col = tumorType)) +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank()) +
  labs(x = "Score ACC", y = "Score PB")
ggsave("Figure S_Probability scores Mixed tumors vs ACC.pdf", path= "./plots/", dpi=500)

### RF score entropy incl TCGA data
anno_all %>% 
  ggplot(aes(tumorType, rf_entropy, fill = tumorType)) +
  geom_boxplot() +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank()) +
  labs(x = NULL, y = "RF score entropy")
ggsave("Figure S_RF score entrophy incl TCGA.pdf", path= "./plots/", dpi=500)

### UMAP incl TCGA and mixed tumors
anno_all %>% 
  ggplot(aes(umap_x, umap_y, col = annotation)) +
  geom_point(aes(color = annotation), shape=19, size = 1) +
  #paletteer::scale_color_paletteer_d("rcartocolor::Safe",
                                   #limits = c("ACC", "Mixed", "NORM", "PanNEC", "PanNET", "PB", "PDAC", "SPN", "TCGA"),
                                   #breaks = c("ACC", "Mixed", "NORM", "PanNEC", "PanNET", "PB", "PDAC", "SPN", "TCGA"),
                                   #labels = c("ACC", "MIXED", "NORM", "PanNEC", "PanNET", "PB", "PDAC", "SPN", "TCGA")) +
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
    #legend.margin=margin(0,0,0,0),
    #legend.box = "horizontal",
    #legend.box.just = "left",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, 'cm'),
    panel.grid = element_blank()) +
  guides(col = guide_legend(nrow = 1))
ggsave("Figure S_Mixed_UMAP-all incl mixed and TCGA_extra.pdf", path= "./Plots/", dpi=500)
