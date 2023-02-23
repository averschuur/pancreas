

anno <- anno %>%
  mutate(conv_score = case_when(conversion >= 1 ~ 'good',
                                 conversion >= 0.8 ~ 'medium',
                                 conversion < 0.8 ~ 'bad'))


anno_goodconv <- anno %>%
  subset(conv_score == "good")


# plot UMAP
anno %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType)) +
  geom_point(size = 4) +
  #geom_text(aes(label = sampleName), size = 4) +
  theme_classic(base_size = 24) +
  scale_color_manual(values = branded_colors1) +
  theme(legend.direction = "horizontal", legend.position = "bottom") +
  labs(x = "UMAP 1", y = "UMAP 2")

anno_goodconv %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType)) +
  geom_point(size = 4) +
  #geom_text(aes(label = sampleName), size = 4) +
  theme_classic(base_size = 24) +
  scale_color_manual(values = branded_colors1) +
  theme(legend.direction = "horizontal", legend.position = "bottom") +
  labs(x = "UMAP 1", y = "UMAP 2")


anno %>% 
  ggplot(aes(umap_x, umap_y, col = conversion)) +
  geom_point(size = 3) +
  paletteer::scale_color_paletteer_c("grDevices::Blue-Red 2") +
  theme_classic(base_size = 24) +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks = element_blank(), 
        legend.direction = "horizontal", 
        legend.position = "bottom", 
        legend.key.width=unit(3,"cm")) +
  labs(x = NULL, y = NULL)

anno_goodconv %>% 
  ggplot(aes(umap_x, umap_y, col = conversion)) +
  geom_point(size = 3) +
  paletteer::scale_color_paletteer_c("grDevices::Blue-Red 2") +
  theme_classic(base_size = 24) +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks = element_blank(), 
        legend.direction = "horizontal", 
        legend.position = "bottom", 
        legend.key.width=unit(3,"cm")) +
  labs(x = NULL, y = NULL)