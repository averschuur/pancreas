###############################
#### FIGURE S_paired cases ####
###############################


### Matched cases biopsy--------------------------------------------------------
#
test %>% 
  as_tibble() %>% 
  slice_sample(n = 50000) %>% 
  ggplot(aes(sample1, sample2)) +
  #geom_bin2d(bins = 50) +
  geom_point(pch = ".", col = branded_colors2[2], alpha = 0.4) +
  scale_fill_distiller(palette = "PuRd", direction = 1) +
  theme_classic(base_size = 30) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, angle = 90, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank())
ggsave("Figure_S_matched cases_b_cor.pdf", path= "./plots/", dpi=500)

# correlation with other samples
test %>% 
  ggplot(aes(label_s, cor, col = as.factor(matched), label=label)) +
  geom_point() +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 30) +
  geom_text(#aes(label=ifelse(matched==1,as.character(label),'')), 
    size =2, hjust =1, vjust = -1) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, angle = 90, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank())
ggsave("Figure_S_matched cases_b_anno.pdf", path= "./plots/", dpi=500)

# correlation between matched cases and other samples
test %>% 
  ggplot(aes(as.factor(matched), cor, fill = matched)) +
  geom_boxplot() +
  #paletteer::scale_fill_paletteer_d("rcartocolor::Safe",
  #limits = c("Other samples", "Matched sample"),
  #breaks = c("Other samples", "Matched sample"),
  #labels = c("Other samples", "Matched sample")) +
  theme_bw(base_size = 18) +
  geom_text(y=0.8, label="P < .001") +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 0, vjust=1, hjust=0.5),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank()) +
  labs(x = NULL, y = "Correlation") +
  scale_x_discrete(labels=c("Other \n samples", "Matched \n sample"))
ggsave("Figure S_paired cases_b_correlation.pdf", path= "./plots/", dpi=500)

### Matched cases metastasis ---------------------------------------------------
#
test %>% 
  as_tibble() %>% 
  slice_sample(n = 50000) %>% 
  ggplot(aes(sample1, sample2)) +
  #geom_bin2d(bins = 50) +
  geom_point(pch = ".", col = branded_colors2[2], alpha = 0.4) +
  scale_fill_distiller(palette = "PuRd", direction = 1) +
  theme_classic(base_size = 30) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, angle = 90, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank())
ggsave("Figure_S_matched cases_m_cor.pdf", path= "./plots/", dpi=500)


# correlation with other samples
test %>% 
  ggplot(aes(label_s, cor, col = as.factor(matched), label=label)) +
  # geom_jitter(width = 0.1) +
  geom_jitter(position = position_jitter(seed = 1),
    #width = 0.1, height = 0.25, size = 1
              ) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 30) +
  geom_text(aes(label=ifelse(matched==1,as.character(label),'')), position = position_jitter(seed = 1), 
            size =2, hjust =1, vjust = -1) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, angle = 90, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank())
ggsave("Figure_S_matched cases_m.pdf", path= "./plots/", dpi=500)

test %>% 
  ggplot(aes(label_s, cor, col = as.factor(matched), label=label)) +
  geom_point() +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 30) +
  geom_text(#aes(label=ifelse(matched==1,as.character(label),'')), 
    size =2, hjust =1, vjust = -1) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, angle = 90, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank())
ggsave("Figure_S_matched cases_m_anno.pdf", path= "./plots/", dpi=500)

# correlation between matched cases and other samples
test %>% 
  ggplot(aes(as.factor(matched), cor, fill = matched)) +
  geom_boxplot() +
  #paletteer::scale_fill_paletteer_d("rcartocolor::Safe",
  #limits = c("Other samples", "Matched sample"),
  #breaks = c("Other samples", "Matched sample"),
  #labels = c("Other samples", "Matched sample")) +
  theme_bw(base_size = 18) +
  geom_text(y=0.8, label="P < 0.01") +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 0, vjust=1, hjust=0.5),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank()) +
  labs(x = NULL, y = "Correlation") +
  scale_x_discrete(labels=c("Other \n samples", "Matched \n sample"))
ggsave("Figure S_paired cases_M_cor.pdf", path= "./plots/", dpi=500)
