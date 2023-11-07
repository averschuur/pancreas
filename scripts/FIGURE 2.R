####################################################
### FIGURE 2: Methods & baseline characteristics ###
####################################################

library(gridExtra)

# Figure 2: UMAP
anno %>% 
  ggplot(aes(umap_x, umap_y)) + 
  geom_point(aes(color = tumorType), shape=19, size = 3) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe",
                                     limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORM"),
                                     breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORM"),
                                     labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORM")) +
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
        legend.box = "horizontal",
        legend.box.just = "left",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank()) +
  guides(col = guide_legend(nrow = 1))
ggsave("Figure 2B_UMAP-all_06112023.pdf", path= "./plots/", dpi=500)





