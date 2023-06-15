####################################################
### FIGURE 2: Methods & baseline characteristics ###
####################################################

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



