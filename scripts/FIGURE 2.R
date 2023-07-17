####################################################
### FIGURE 2: Methods & baseline characteristics ###
####################################################

library(gridExtra)


# Figure 2A: count tumorType per source
anno %>% 
  ggplot(aes(tumorType, fill = tumorType)) +
  geom_bar() +
  labs(#title="Tumor Type By Source",
    x = "", y = "No. of cases", fill = "Source") +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal")) +
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
ggsave("Figure 2A_count tumorType.pdf", path= "./plots/", dpi=500)



# Figure 2B: UMAP
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
        #legend.margin=margin(0,0,0,0),
        legend.box = "horizontal",
        legend.box.just = "left",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank()) +
  guides(col = guide_legend(nrow = 1))
ggsave("Figure 2B_UMAP-all.pdf", path= "./plots/", dpi=500)


# Figure 2C umap - avg beta
anno %>% 
  ggplot(aes(umap_x, umap_y, col = avg_beta_unfiltered)) + 
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
        #legend.margin=margin(0,0,0,0),
        legend.title = element_text(face = "bold", size = 10, vjust = 0.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.6, 'cm'),
        panel.grid = element_blank()) +
  guides(colour=guide_colourbar(barwidth=15,label.position="bottom", title.vjust=0.9, label.vjust=2))
ggsave("Figure 2C_UMAP-avg_beta.pdf", path= "./plots/", dpi=500)

# Figure 2D umap - conversion
anno %>% 
  ggplot(aes(umap_x, umap_y, col = conversion)) + 
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
    #legend.margin=margin(0,0,0,0),
    legend.title = element_text(face = "bold", size = 10, vjust = 0.5),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.6, 'cm'),
    panel.grid = element_blank()) +
  scale_color_stepsn(colors = c("#9FCB9B", "#7DACCF", "#41709C", "#2A5783"), breaks = c(0.8, 1.0, 1.2)) +
  guides(colour=guide_colourbar(barwidth=15,label.position="bottom", title.vjust=0.9, label.vjust=2))
  scale_fill_steps(breaks = c(0.8, 1.0, 1.2), guide = guide_coloursteps(even.steps = FALSE)
  )
ggsave("Figure 1D_UMAP-conversion.pdf", path= "./plots/", dpi=500)

# Figure 2E umap - estimate
anno %>% 
  ggplot(aes(umap_x, umap_y, col = estimate)) + 
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
    #legend.margin=margin(0,0,0,0),
    legend.title = element_text(face = "bold", size = 10, vjust = 0.5),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.6, 'cm'),
    panel.grid = element_blank()) +
  guides(colour=guide_colourbar(barwidth=15,label.position="bottom", title.vjust=0.9, label.vjust=2))
ggsave("Figure 1E_UMAP-estimae.pdf", path= "./plots/", dpi=500)

# Figure 2F umap - absolute
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
        #legend.margin=margin(0,0,0,0),
        legend.title = element_text(face = "bold", size = 10, vjust = 0.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.6, 'cm'),
        panel.grid = element_blank()) +
  guides(colour=guide_colourbar(barwidth=15,label.position="bottom", title.vjust=0.9, label.vjust=2))
ggsave("Figure 1F_UMAP-absolute.pdf", path= "./plots/", dpi=500)



