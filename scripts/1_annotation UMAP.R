#Annotation UMAP


annotation <- data.frame(
  x = c(0.15296519,-4.9704108, 0.27748376, 4.02438531, 7.36728043, 4.04303804, -3.4483243),
  y = c(-2.2111147,-1.48372349, -1.0113779, -0.8714847, -3.2585769, -0.4226539, -0.62975666),
  label = c("SPN9T", "13T", "14T", "ACC K11 T", "PNET 26T", "PNET39T", "PNET-NE040")
)


anno %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType)) + 
  geom_point(alpha = 0.7) +
  geom_text(data=annotation, aes(x=x, y=(y-0.15), label=label),                 , 
            color="black", 
            size=2 , fontface="bold" ) +
  annotate("text", x = 4.03022761, y = (-0.9404484-0.3), label = "ACC K15 T",
    color="black", 
    size=2 , 
    fontface="bold" ) +
  scale_colour_manual(values = branded_colors) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 12)
ggsave("UMAP-all_with annotation.png", path= "./output/")
