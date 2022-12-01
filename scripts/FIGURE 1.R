################
### FIGURE 1 ###
################

library(gridExtra)


# Figure 1B: count tumorType per source
anno %>% 
  ggplot(aes(tumorType)) +
  geom_bar(aes(fill = source)) +
  labs(#title="Tumor Type By Source",
    x = "", y = "No. of cases", fill = "Source") +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", labels = c("Benhamida et al.", "Chan et al.", "DiDomenico et al.", "Endo et al.", "Jäkel et al.", "Selenica et al.", "Yachida et al.")) +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
        axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
        axis.ticks = element_blank(),
        legend.position = c(.83, .79),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(0.4, 'cm'),
        panel.grid = element_blank())
ggsave("Figure 1B_count tumorType per source_colorSafe.pdf", path= "./output/", dpi=500)


# Figure 1C: average methylation per TumorType
anno %>% 
  ggplot(aes(tumorType, avg_beta, fill=tumorType, col=tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =2) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN")) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe", breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN")) +
  labs(#title="Average Methylation By Tumor Type",
       x = "", 
       y = "Average Methylation") +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
        axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.position = "none",
        panel.grid = element_blank())
ggsave("Figure 1C_average methylation_ColorSafe.pdf", path= "./output/", dpi=500)


# Figure 1D: tumor purity per tumorType_estimate
anno %>% 
  ggplot(aes(tumorType, estimate, fill=tumorType, col=tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =2,alpha=0.7) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN")) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe", breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN")) +
  labs(title="Estmate Tumor Purity By Tumor Type",
       x = "", 
       y = "Estimate Tumor Purity") +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 14),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
        axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.position = "none",
        panel.grid = element_blank())
ggsave("Figure 1D_tumor purity per tumorType_estimate_colorSafe.pdf", path= "./output/", dpi=500)


# Figure 1E: tumor purity per tumorType_absolute
anno %>% 
  ggplot(aes(tumorType, absolute, fill=tumorType, col=tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =2) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN")) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe", breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN")) +
  labs(#title="Absolute Tumor Purity By Tumor Type",
       x = "", 
       y = "Absolute Tumor Purity") +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
        axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.position = "none",
        panel.grid = element_blank())
ggsave("Figure 1E_tumor purity per tumorType_absolute_ColorSafe.pdf", path= "./output/", dpi=500)

# Figure 1F umap
anno %>% 
  ggplot(aes(umap_x, umap_y)) + 
  geom_point(aes(color = tumorType), shape=19, size = 4) +
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
        legend.margin=margin(0,0,0,0),
        legend.box = "horizontal",
        legend.box.just = "left",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank()) +
  guides(col = guide_legend(nrow = 1))
ggsave("Figure 1F_UMAP-all_colorSafe.pdf", path= "./output/", dpi=500)


# Figure 1F umap - avg beta
anno %>% 
  ggplot(aes(umap_x, umap_y, col = avg_beta)) + 
  geom_point(shape=19, size = 5) +
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
        legend.margin=margin(0,0,0,0),
        legend.title = element_text(face = "bold", size = 10, vjust = 0.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.6, 'cm'),
        panel.grid = element_blank()) +
  guides(colour=guide_colourbar(barwidth=15,label.position="bottom", title.vjust=0.9, label.vjust=2))
ggsave("Figure 1G_UMAP-avg_beta_colorSafe.pdf", path= "./output/", dpi=500)

# Figure 1F umap - absolute
anno %>% 
  ggplot(aes(umap_x, umap_y, col = absolute)) + 
  geom_point(shape=19, size = 5) +
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
        legend.margin=margin(0,0,0,0),
        legend.title = element_text(face = "bold", size = 10, vjust = 0.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.6, 'cm'),
        panel.grid = element_blank()) +
  guides(colour=guide_colourbar(barwidth=15,label.position="bottom", title.vjust=0.9, label.vjust=2))
ggsave("Figure 1H_UMAP-absolute_colorSafe.pdf", path= "./output/", dpi=500)



