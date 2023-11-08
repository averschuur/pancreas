#####################################
### S_extra unsupervised_analysis ###
#####################################

library(paletteer)

### Tumor purity (estimate) -------------------------------------------------------
# Figure S: estimate tumor purity - umap
anno %>% 
  ggplot(aes(umap_x, umap_y, col = estimate)) + 
  geom_point(shape=19, size = 4) +
  paletteer::scale_color_paletteer_c("ggthemes::Green-Blue Diverging") +
  labs(#title="UMAP Clustering",
    x = "UMAP 1", 
    y = "UMAP 2",
    col = "Absolute Tumor Purity") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 12),
    axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=0, size = 10),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
    axis.text.x = element_text(face = "bold", size = 9),
    axis.text.y = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    #legend.margin=margin(0,0,0,0),
    legend.title = element_text(face = "bold", size = 10, vjust = 0.5),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.6, 'cm'),
    panel.grid = element_blank()) +
  guides(colour=guide_colourbar(barwidth=15,label.position="bottom", title.vjust=0.9, label.vjust=2))
ggsave("Figure S_UMAP-estimae_07222023.pdf", path= "./plots/", dpi=500)

# Figure S: estimate tumor purity - boxplot
anno %>% 
  ggplot(aes(tumorType, estimate, fill= tumorType, col = tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =2,alpha=0.7)+
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  labs(x = "Tumor Type", y = "Purity (ESTIMATE)") +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    panel.grid = element_blank(),
    legend.position = "none")
ggsave("Figure S_tumortype estimate boxplot_07222023.pdf", path= "./plots/", dpi=500)



### Tumor purity (absolute) -------------------------------------------------------
# Figure S: absolute tumor purity - umap
anno %>% 
  ggplot(aes(umap_x, umap_y, col = absolute)) + 
  geom_point(shape=19, size = 4) +
  paletteer::scale_color_paletteer_c("ggthemes::Green-Blue Diverging") +
  labs(#title="UMAP Clustering",
    x = "UMAP 1", 
    y = "UMAP 2",
    col = "Absolute Tumor Purity") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 12),
    axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=0, size = 10),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
    axis.text.x = element_text(face = "bold", size = 9),
    axis.text.y = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    #legend.margin=margin(0,0,0,0),
    legend.title = element_text(face = "bold", size = 10, vjust = 0.5),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.6, 'cm'),
    panel.grid = element_blank()) +
  guides(colour=guide_colourbar(barwidth=15,label.position="bottom", title.vjust=0.9, label.vjust=2))
ggsave("Figure S_UMAP-absolute_07222023.pdf", path= "./plots/", dpi=500)

# Figure S: absolute tumor purity - boxplot
anno %>% 
  ggplot(aes(tumorType, absolute, fill= tumorType, col = tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =2,alpha=0.7)+
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  labs(x = "Tumor Type", y = "Purity (ABSOLUTE)") +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    panel.grid = element_blank(),
    legend.position = "none")
ggsave("Figure S_tumortype absolute boxplot_07222023.pdf", path= "./plots/", dpi = 500)



### Average methylation -------------------------------------------------------
# Figure S: avg beta - umap
anno %>% 
  ggplot(aes(umap_x, umap_y, col = avg_beta_unfiltered)) + 
  geom_point(shape=19, size = 4) +
  paletteer::scale_color_paletteer_c("ggthemes::Green-Blue Diverging") +
  labs(#title="UMAP Clustering",
    x = "UMAP 1", 
    y = "UMAP 2",
    col = "Average Beta Values") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 12),
    axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=0, size = 10),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
    axis.text.x = element_text(face = "bold", size = 9),
    axis.text.y = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    #legend.margin=margin(0,0,0,0),
    legend.title = element_text(face = "bold", size = 10, vjust = 0.5),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.6, 'cm'),
    panel.grid = element_blank()) +
  guides(colour=guide_colourbar(barwidth=15,label.position="bottom", title.vjust=0.9, label.vjust=2))
ggsave("Figure S_UMAP-avg_beta_07222023.pdf", path= "./plots/", dpi=500)

# average methylation - boxplot
anno %>% 
  ggplot(aes(tumorType, avg_beta_unfiltered, fill= tumorType, col = tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =2,alpha=0.7)+
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  labs(x = "Tumor Type", y = "Average methylation") +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    panel.grid = element_blank(),
    legend.position = "none")
ggsave("Figure S_tumortype avg methylation boxplot_07222023.pdf", path= "./plots/", dpi = 500)



### Conversion -------------------------------------------------------
# Figure S: conversion - umap
anno %>% 
  ggplot(aes(umap_x, umap_y, col = conversion)) + 
  geom_point(shape=19, size = 4) +
  paletteer::scale_color_paletteer_c("ggthemes::Green-Blue Diverging") +
  labs(#title="UMAP Clustering",
    x = "UMAP 1", 
    y = "UMAP 2",
    col = "Absolute Tumor Purity") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 12),
    axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=0, size = 10),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
    axis.text.x = element_text(face = "bold", size = 9),
    axis.text.y = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    #legend.margin=margin(0,0,0,0),
    legend.title = element_text(face = "bold", size = 10, vjust = 0.5),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.6, 'cm'),
    panel.grid = element_blank()) +
  scale_color_stepsn(colors = c("#9FCB9B", "#7DACCF", "#41709C", "#2A5783"), breaks = c(0.8, 1.0, 1.2)) +
  guides(colour=guide_colourbar(barwidth=15,label.position="bottom", title.vjust=0.9, label.vjust=2)) +
  scale_fill_steps(breaks = c(0.8, 1.0, 1.2), guide = guide_coloursteps(even.steps = FALSE))
ggsave("Figure S_UMAP-conversion_07222023.pdf", path= "./plots/", dpi=500)

# Figure S: conversion - boxplot
anno %>% 
  ggplot(aes(tumorType, conversion, fill= tumorType, col = tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =2,alpha=0.7)+
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  labs(x = "Tumor Type", y = "Conversion Score") +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    panel.grid = element_blank(),
    legend.position = "none")
ggsave("Figure S_tumortype conversion score boxplot_07222023.pdf", path= "./plots/", dpi = 500)



### Tissue effect - umap ---------------------------------------------------
# tissues
tissues <- read_csv(file = "./annotation/annotation_tissues.csv")
anno <- left_join(anno, tissues)


anno %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType, shape = array)) + 
  geom_point(aes(color = tumorType), size = 3) +
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
    #legend.margin=margin(0,0,0,0),
    #legend.box = "horizontal",
    #legend.box.just = "left",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, 'cm'),
    panel.grid = element_blank()) +
  guides(col = guide_legend(nrow = 1))
ggsave("Figure S_umap_tissues_07112023.pdf", path= "./plots/", dpi=500)

