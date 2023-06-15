
### normal tissue -----------------------------------------------------------------------------------
# plot sample numbers per study
anno_normal %>% 
  group_by(source) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(source, n, fill = source)) +
  geom_col(width = 0.4) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", 
                                    limits = c("DiDomenico", "Jakel"),
                                    breaks = c("DiDomenico", "Jakel"),
                                    labels = c("DiDomenico", "Jakel")) +
  labs(x = NULL, y = "Normal tissues (n)") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("Figure S_Batch effect_normal_count.pdf", path= "./Plots/", dpi=500)

# density plot
betas_normal_avg %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(value, fill = name)) +
  geom_density(alpha = 0.8) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "none",
        #panel.grid = element_blank(),
        strip.text = element_text(size = 10)) +
  labs(title=NULL,x = "Beta", y = "Density") +
  facet_grid(rows = vars(name))
ggsave("FigureS_batch effect_normal tissue_density.pdf", path= "./plots/", dpi=500)


# mean beta values
anno_normal %>% 
  ggplot(aes(source, mean_beta, fill = source)) +
  geom_boxplot(alpha = 0.8) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", 
                                    limits = c("DiDomenico", "Jakel"),
                                    breaks = c("DiDomenico", "Jakel"),
                                    labels = c("DiDomenico", "Jakel")) +
  labs(x = "", y = "Mean beta values") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("FigureS_batch effect_normal_mean betas.pdf", path= "./plots/", dpi=500)


# plot control score
anno_normal %>% 
  ggplot(aes(source, conversion, fill = source)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(size = 2, width = 0.1) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", 
                                    limits = c("DiDomenico", "Jakel"),
                                    breaks = c("DiDomenico", "Jakel"),
                                    labels = c("DiDomenico", "Jakel")) +
  labs(x = NULL, y = "Normal tissues (n)") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("Figure S_Batch effect_normal_controlscore.pdf", path= "./Plots/", dpi=500)

# umap
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

# plot control score
anno_pnet %>% 
  ggplot(aes(source, conversion, fill = source)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(size = 2, width = 0.1) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", 
                                    limits = c("Chan", "DiDomenico", "Jakel", "UMCU", "yachida"),
                                    breaks = c("Chan", "DiDomenico", "Jakel", "UMCU", "yachida"),
                                    labels = c("Chan", "DiDomenico", "Jakel", "UMCU", "yachida")) +
  labs(x = NULL, y = "Normal tissues (n)") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("Figure S_Batch effect_PNET_controlscore.pdf", path= "./Plots/", dpi=500)

betas_pnet_avg %>% 
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
ggsave("FigureS_batch effect_pnet_violinplot.pdf", path= "./plots/", dpi=500)

# mean beta values
anno_pnet %>% 
  ggplot(aes(source, mean_beta, fill = source)) +
  geom_boxplot(alpha = 0.8) +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 16),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=1, size = 16),
        axis.text = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = "none") +
  labs(x = "", y = "Mean beta values") +
  scale_fill_manual(values = branded_colors1)
ggsave("FigureS_batch effect_pnet_mean betas.pdf", path= "./plots/", dpi=500)


# density plot of beta value distribution
betas_pnet_avg %>%
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(value, fill = name)) +
  geom_density(alpha = 0.8) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "none",
        #panel.grid = element_blank(),
        strip.text = element_text(size = 10)) +
  labs(title=NULL,x = "Beta", y = "Density") +
  facet_grid(rows = vars(name))
ggsave("FigureS_batch effect_pnet_density.pdf", path= "./plots/", dpi=500)

# umap
anno_pnet %>% 
  ggplot(aes(umap_x, umap_y)) +
  geom_point(aes(color = source, shape = source), shape=19, size = 4) +
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
ggsave("Figure S_Batch effect_pnet_UMAP.pdf", path= "./Plots/", dpi=500)

### ACC tissue -----------------------------------------------------------------------------------
# plot sample numbers per study
anno_acc %>% 
  group_by(source) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(source, n, fill = source)) +
  geom_col(width = 0.4) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", 
                                    limits = c("Benhamida", "Jakel", "UMCU"),
                                    breaks = c("Benhamida", "Jakel", "UMCU"),
                                    labels = c("Benhamida", "Jakel", "UMCU")) +
  labs(x = NULL, y = "ACC tissues (n)") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("Figure S_Batch effect_acc_count.pdf", path= "./Plots/", dpi=500)

# plot control score
anno_acc %>% 
  ggplot(aes(source, conversion, fill = source)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(size = 2, width = 0.1) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", 
                                    limits = c("Benhamida", "Jakel", "UMCU"),
                                    breaks = c("Benhamida", "Jakel", "UMCU"),
                                    labels = c("Benhamida", "Jakel", "UMCU")) +
  labs(x = NULL, y = "ACC tissues (n)") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("Figure S_Batch effect_acc_controlscore.pdf", path= "./Plots/", dpi=500)


# mean beta values
anno_acc %>% 
  ggplot(aes(source, mean_beta, fill = source)) +
  geom_boxplot(alpha = 0.8) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", 
                                    limits = c("Benhamida", "Jakel", "UMCU"),
                                    breaks = c("Benhamida", "Jakel", "UMCU"),
                                    labels = c("Benhamida", "Jakel", "UMCU")) +
  labs(x = "", y = "Mean beta values") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("FigureS_batch effect_ACC_mean betas.pdf", path= "./plots/", dpi=500)


# density plot of beta value distribution
betas_acc_avg %>%
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(value, fill = name)) +
  geom_density(alpha = 0.8) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "none",
        #panel.grid = element_blank(),
        strip.text = element_text(size = 10)) +
  labs(title=NULL,x = "Beta", y = "Density") +
  facet_grid(rows = vars(name))
ggsave("FigureS_batch effect_acc_density.pdf", path= "./plots/", dpi=500)

# umap
anno_acc %>% 
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
ggsave("Figure S_Batch effect_acc_UMAP.pdf", path= "./Plots/", dpi=500)

### SPN tissue -----------------------------------------------------------------------------------
# plot sample numbers per study
anno_spn %>% 
  group_by(source) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(source, n, fill = source)) +
  geom_col(width = 0.4) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", 
                                    limits = c("Selenica", "UMCU"),
                                    breaks = c("Selenica", "UMCU"),
                                    labels = c("Selenica", "UMCU")) +
  labs(x = NULL, y = "ACC tissues (n)") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("Figure S_Batch effect_spn_count.pdf", path= "./Plots/", dpi=500)

# plot control score
anno_spn %>% 
  ggplot(aes(source, conversion, fill = source)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(size = 2, width = 0.1) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", 
                                    limits = c("Selenica", "UMCU"),
                                    breaks = c("Selenica", "UMCU"),
                                    labels = c("Selenica", "UMCU")) +
  labs(x = NULL, y = "ACC tissues (n)") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("Figure S_Batch effect_spn_controlscore.pdf", path= "./Plots/", dpi=500)


# mean beta values
anno_spn %>% 
  ggplot(aes(source, mean_beta, fill = source)) +
  geom_boxplot(alpha = 0.8) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", 
                                    limits = c("Selenica", "UMCU"),
                                    breaks = c("Selenica", "UMCU"),
                                    labels = c("Selenica", "UMCU")) +
  labs(x = "", y = "Mean beta values") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("FigureS_batch effect_spn_mean betas.pdf", path= "./plots/", dpi=500)


# density plot of beta value distribution
betas_spn_avg %>%
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(value, fill = name)) +
  geom_density(alpha = 0.8) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "none",
        #panel.grid = element_blank(),
        strip.text = element_text(size = 10)) +
  labs(title=NULL,x = "Beta", y = "Density") +
  facet_grid(rows = vars(name))
ggsave("FigureS_batch effect_spn_density.pdf", path= "./plots/", dpi=500)

# umap
anno_spn %>% 
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
ggsave("Figure S_Batch effect_spn_UMAP.pdf", path= "./Plots/", dpi=500)

