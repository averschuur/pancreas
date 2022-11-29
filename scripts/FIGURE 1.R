################
### FIGURE 1 ###
################

library(gridExtra)


# Figure 1B: count tumorType per source
anno %>% 
  ggplot(aes(tumorType)) +
  geom_bar(aes(fill = source)) +
  labs(title="Tumor Type By Source",x = "", y = "No. of cases", fill = "Source") +
  scale_fill_manual(values = branded_colors, labels = c("Benhamida et al.", "Chan et al.", "DiDomenico et al.", "Endo et al.", "Jäkel et al.", "Selenica et al.", "Yachida et al.")) +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, angle = 90, vjust=0.1, hjust=1),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = c(.85, .78),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank())
ggsave("Figure 1B_count tumorType per source.png", path= "./output/")


# Figure 1C: average methylation per TumorType
anno %>% 
  ggplot(aes(tumorType, avg_beta, fill=tumorType, col=tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =2,alpha=0.7) +
  scale_fill_manual(breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN"), 
                    values=branded_colors) +
  scale_color_manual(breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN"), 
                     values=branded_colors) +
  labs(title="Average Methylation By Tumor Type",
       x = "", 
       y = "Estimate Tumor Purity") +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, angle = 90, vjust=0.1, hjust=1),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = "none",
        panel.grid = element_blank())
ggsave("Figure 1C_average methylation.png", path= "./output/")


# Figure 1D: tumor purity per tumorType_estimate
anno %>% 
  ggplot(aes(tumorType, estimate, fill=tumorType, col=tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =2,alpha=0.7) +
  scale_fill_manual(breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN"), 
                    values=branded_colors) +
  scale_color_manual(breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN"), 
                    values=branded_colors) +
  labs(title="Estmate Tumor Purity By Tumor Type",
       x = "", 
       y = "Estimate Tumor Purity") +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, angle = 90, vjust=0.1, hjust=1),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = "none",
        panel.grid = element_blank())
ggsave("Figure 1D_tumor purity per tumorType_estimate.png", path= "./output/")


# Figure 1E: tumor purity per tumorType_absolute
anno %>% 
  ggplot(aes(tumorType, absolute, fill=tumorType, col=tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =2,alpha=0.7) +
  scale_fill_manual(breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN"), 
                    values=branded_colors) +
  scale_color_manual(breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN"), 
                     values=branded_colors) +
  labs(title="Absolute Tumor Purity By Tumor Type",
       x = "", 
       y = "Absolute Tumor Purity") +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, angle = 90, vjust=0.1, hjust=1),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = "none",
        panel.grid = element_blank())
ggsave("Figure 1E_tumor purity per tumorType_absolute.png", path= "./output/")

# Figure 1F umap
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
ggsave("UMAP-all.png", path= "./output/")

# Figure 1E t-SNE
anno %>% 
  ggplot(aes(tsne_x, tsne_y, col = tumorType)) + 
  geom_point(shape=19, size = 3, alpha = 0.7) +
  scale_color_manual(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                     breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                     values = branded_colors, 
                     labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  labs(title="t-SNE Clustering",
       x = "t-SNE 1", 
       y = "t-SNE 2",
       col = "Tumor Type") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 8),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = c(.89, .81),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank())
ggsave("tSNE-all.png", path= "./output/")


## Arrange Plots
lay <- rbind(c(1,5,5),
             c(2,6,6),
             c(3,6,6),
             c(4,6,6))

grid.arrange(arrangeGrob(G1, G2, G3, G4, nrow = 4), arrangeGrob(G5, G6, nrow = 2), ncol =2, layout_matrix = lay)

gs <- lapply(c("G1", "G2", "G3", "G4", "G5", "G6"), function(ii)
  grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))

grid.arrange(grobs = gs, layout_matrix = lay)

plot_grid(ggarrange(G1, G2, G3, G4,
                    labels = c("B", "C", "D", "E"),
                    ncol = 1, nrow = 4),
          ggarrange(G5, G6, 
                    labels = c("F", "G"),
                    ncol=1, nrow = 2))
### ------------------------------

# Figure 1A: count tumorType per source
G1 <- anno %>% 
  ggplot(aes(tumorType)) +
  geom_bar(aes(fill = source)) +
  labs(title="Tumor Type By Source",x = "", y = "No. of cases", fill = "Source") +
  scale_fill_manual(values = branded_colors, labels = c("Benhamida et al.", "Chan et al.", "DiDomenico et al.", "Endo et al.", "Jäkel et al.", "Selenica et al.", "Yachida et al.")) +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", hjust = 0, vjust = -1, size = 7),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=-0.5, size = 6),
        axis.text.x = element_text(face = "bold", size = 4, angle = 90, vjust=0.5, hjust=1),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y =unit(.05, "cm"),
        axis.text.y = element_text(face = "bold", size = 4),
        legend.position = c(.74, .57),
        legend.title = element_text(face = "bold", size = 6),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.1, 'cm'),
        panel.grid = element_blank())


# Figure 1C: average methylation per TumorType
G2 <- anno %>% 
  ggplot(aes(tumorType, avg_beta, fill=tumorType, col=tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =1,alpha=0.7) +
  scale_fill_manual(breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN"), 
                    values=branded_colors) +
  scale_color_manual(breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN"), 
                     values=branded_colors) +
  labs(title="Average Methylation By Tumor Type",
       x = "", 
       y = "Average Methylation") +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", hjust = 0, vjust = -1, size = 7),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=-0.5, size = 6),
        axis.text.x = element_text(face = "bold", size = 4, angle = 90, vjust=0.5, hjust=1),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y =unit(.05, "cm"),
        axis.text.y = element_text(face = "bold", size = 4),
        legend.position = "none",
        panel.grid = element_blank())


# Figure 1D: tumor purity per tumorType_estimate
G3 <- anno %>% 
  ggplot(aes(tumorType, estimate, fill=tumorType, col=tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =1,alpha=0.7) +
  scale_fill_manual(breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN"), 
                    values=branded_colors) +
  scale_color_manual(breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN"), 
                     values=branded_colors) +
  labs(title="Estimate Tumor Purity By Tumor Type",
       x = "", 
       y = "Estimate Tumor Purity") +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", hjust = 0, vjust = -1, size = 7),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=-0.5, size = 6),
        axis.text.x = element_text(face = "bold", size = 4, angle = 90, vjust=0.5, hjust=1),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y =unit(.05, "cm"),
        axis.text.y = element_text(face = "bold", size = 4),
        legend.position = "none",
        panel.grid = element_blank())


# Figure 1D: tumor purity per tumorType_absolute
G4 <- anno %>% 
  ggplot(aes(tumorType, absolute, fill=tumorType, col=tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =1,alpha=0.7) +
  scale_fill_manual(breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN"), 
                    values=branded_colors) +
  scale_color_manual(breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN"), 
                     values=branded_colors) +
  labs(title="Absolute Tumor Purity By Tumor Type",
       x = "", 
       y = "Absolute Tumor Purity") +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", hjust = 0, vjust = -1, size = 7),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=-0.5, size = 6),
        axis.text.x = element_text(face = "bold", size = 4, angle = 90, vjust=0.5, hjust=1),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y =unit(.05, "cm"),
        axis.text.y = element_text(face = "bold", size = 4),
        legend.position = "none",
        panel.grid = element_blank())

# Figure .. umap
G6 <- anno %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType)) + 
  geom_point(shape=19, size = 2, alpha = 0.7) +
  scale_color_manual(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                     breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                     values = branded_colors, 
                     labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  labs(title="UMAP Clustering",
       x = "UMAP 1", 
       y = "UMAP 2",
       col = "Tumor Type") +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", hjust = 0, vjust = -1, size = 7),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 6),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=-0.5, size = 6),
        axis.text.x = element_text(face = "bold", size = 4, vjust=4, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 4),
        axis.ticks.length=unit(.05, "cm"),
        legend.position = c(.88, .84),
        legend.title = element_text(face = "bold", size = 6),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.05, 'cm'),
        panel.grid = element_blank())

# Figure .. t-SNE
G5 <- anno %>% 
  ggplot(aes(tsne_x, tsne_y, col = tumorType)) + 
  geom_point(shape=19, size = 2, alpha = 0.7) +
  scale_color_manual(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                     breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                     values = branded_colors, 
                     labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  labs(title="t-SNE Clustering",
       x = "t-SNE 1", 
       y = "t-SNE 2",
       col = "Tumor Type") +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", hjust = 0, vjust = -1, size = 7),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 6),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=-0.5, size = 6),
        axis.text.x = element_text(face = "bold", size = 4, vjust=4, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 4),
        axis.ticks.length =unit(.05, "cm"),
        legend.position = c(.88, .84),
        legend.title = element_text(face = "bold", size = 6),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.05, 'cm'),
        panel.grid = element_blank())

Figure1 <- ggdraw() +
  draw_plot(G1, x = 0, y = .72, width = .4, height = .3) +
  draw_plot(G2, x = 0, y = .47, width = .4, height = .3) +
  draw_plot(G3, x = 0, y = .22, width = .4, height = .3) +
  draw_plot(G4, x = 0, y = -0.03, width = .4, height = .3) +
  draw_plot(G5, x = .39, y = .51, width = .6, height = .51) +
  draw_plot(G6, x = .39, y = 0.01, width = .6, height = .51) +
  draw_plot_label(label = c("B", "C", "D", "E", "F", "G"), size = 10,
                  x = c(.01, .01, .01, .01, .41, .41), y = c(1, .765, .505, .255, 1, 0.51))
ggsave("FIGURE1-TEST.png", Figure1, path= "./output/")


