
### normal tissue -----------------------------------------------------------------------------------
betas_normal_avg %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(name, value, fill = name)) +
  geom_violin(alpha = 0.8) +
  theme(legend.position = "none")+
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "none",
        panel.grid = element_blank()) +
  labs(title=NULL,x = "Beta", y = "Density") +
  theme(legend.position = "none") +
  ylim(0, 1)
ggsave("FigureS_batch effect_normal tissue_density.pdf", path= "./plots/", dpi=500)


betas_normal_avg %>% 
  sample_n(size = 30000) %>% 
  ggplot(aes(DiDomenico, Jakel)) +
  geom_bin2d(bins = 100) +
  paletteer::scale_fill_paletteer_c("scico::tokyo") +
  #paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "none",
        panel.grid = element_blank()) +
  geom_abline(slope = 1, intercept = 0, col = "grey", lty = 2) +
  theme(legend.position = "none")
ggsave("FigureS_batch effect_normal tissue_avg betas.pdf", path= "./plots/", dpi=500)


anno_normal %>% 
  ggplot(aes(umap_x, umap_y)) +
  geom_point(aes(color = source), shape=19, size = 4) +
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
ggsave("Figure S_Batch effect_normal_UMAP.pdf", path= "./Plots/", dpi=500)

### PNETs -----------------------------------------------------------------------

# plot sample numbers per study
anno_pnet %>% 
  group_by(source) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(source, n, fill = source)) +
  geom_col(width = 0.4) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", 
                                    limits = c("Chan", "DiDomenico", "Jakel", "UMCU", "yachida"),
                                    breaks = c("Chan", "DiDomenico", "Jakel", "UMCU", "yachida"),
                                    labels = c("Chan", "DiDomenico", "Jakel", "UMCU", "Yachida")) +
  labs(x = NULL, y = "Normal tissues (n)") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("Figure S_Batch effect_PNET_count.pdf", path= "./Plots/", dpi=500)


