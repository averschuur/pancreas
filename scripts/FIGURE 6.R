##################################
### correlations matched cases ###
##################################

### biopsy cases

anno_paired_b <- anno_paired %>%
  filter(location_2 == "biopsy")

i = 4


test_b <- betas[, c(anno_paired_b$array_id_1[i], anno_paired_b$array_id_2[i])]
colnames(test_b) <- c("sample1", "sample2")
dim(test_b)

test_b %>% 
  as_tibble() %>% 
  slice_sample(n = 50000) %>% 
  ggplot(aes(sample1, sample2)) +
  #geom_bin2d(bins = 50) +
  geom_point(pch = ".", col = branded_colors2[2], alpha = 0.4) +
  scale_fill_distiller(palette = "PuRd", direction = 1) +
  theme_classic(base_size = 20)

test_b <- cor(betas[, anno_paired$array_id_1], betas[, anno_paired$array_id_2])
#test[lower.tri(test)] <- NA
test_b <- test_b %>% 
  as_tibble(rownames = "source") %>% 
  pivot_longer(cols = 2:10, names_to = "target", values_to = "cor")

test_b2 <- anno_paired %>% 
  select(array_id_1, array_id_2) %>% 
  rename(array_id_1 = "source", array_id_2 = "matched")

test_b <- test_b %>% 
  left_join(test_b2)

test_b <- test_b %>% 
  mutate(matched = ifelse(matched == target, 1, 0))

test_b %>% 
  ggplot(aes(source, cor, col = as.factor(matched))) +
  geom_jitter(width = 0.1) +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, vjust=0, hjust=0.5),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.4, 'cm'),
    panel.grid = element_blank())
ggsave("matched cases.pdf", path= "./plots/", dpi=500)

test_b %>% 
  ggplot(aes(as.factor(matched), cor)) +
  geom_boxplot() +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, vjust=0, hjust=0.5),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.4, 'cm'),
    panel.grid = element_blank())
ggsave("biopsy correlation.pdf", path= "./plots/", dpi=500)

# classification correct?
anno_biopsy <- anno %>%
  filter(source == "UMCU") %>%
  filter(sampleName %in% c("UMCU_panNET7_b", "UMCU_panNET7", "UMCU_panNET8_b", "UMCU_panNET8", "RB_ACC1_b", "RB_ACC1", "RB_SPN1_b", "RB_SPN1"))

anno_biopsy$pred_rf <- factor(anno_biopsy$pred_rf, levels = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"))
anno_biopsy$tumorType <- factor(anno_biopsy$tumorType, levels = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"))

plt <- caret::confusionMatrix(anno_biopsy$pred_rf,anno_biopsy$tumorType)

plt <- as.data.frame(plt$table)

ggplot(plt, aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq), color = "white") +
  scale_fill_gradient(low = "#132B43", high = "#88CCEEFF") +
  #paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  labs(x = "True Class ",y = "Predicted Class (NN)") +
  scale_x_discrete(limits = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"),
                   breaks = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"),
                   labels = c("PanNET", "PDAC", "ACC", "SPN", "NORMAL", "PB", "PanNEC")) +
  theme_bw(base_size = 16) +
  theme(axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks=element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank())
ggsave("confusion matrix biopsy.pdf", path= "./plots/", dpi=500)

### metastasized cases
anno_paired_m <- anno_paired %>%
  filter(location_2 == "metastasis") %>%
  filter(label != "Mixed")

i = 4


test_m <- betas[, c(anno_paired_b$array_id_1[i], anno_paired_b$array_id_2[i])]
colnames(test_m) <- c("sample1", "sample2")
dim(test_m)

test_m %>% 
  as_tibble() %>% 
  slice_sample(n = 50000) %>% 
  ggplot(aes(sample1, sample2)) +
  #geom_bin2d(bins = 50) +
  geom_point(pch = ".", col = branded_colors2[2], alpha = 0.4) +
  scale_fill_distiller(palette = "PuRd", direction = 1) +
  theme_classic(base_size = 20)

test_m <- cor(betas[, anno_paired$array_id_1], betas[, anno_paired$array_id_2])
#test[lower.tri(test)] <- NA
test_m <- test_m %>% 
  as_tibble(rownames = "source") %>% 
  pivot_longer(cols = 2:10, names_to = "target", values_to = "cor")

test_m2 <- anno_paired %>% 
  select(array_id_1, array_id_2) %>% 
  rename(array_id_1 = "source", array_id_2 = "matched")

test_m <- test_m %>% 
  left_join(test_m2)

test_m <- test_m %>% 
  mutate(matched = ifelse(matched == target, 1, 0))

test_m %>% 
  ggplot(aes(source, cor, col = as.factor(matched))) +
  geom_jitter(width = 0.1) +
  theme_classic(base_size = 20)
ggsave("matched cases.pdf", path= "./plots/", dpi=500)

test_m %>% 
  ggplot(aes(as.factor(matched), cor)) +
  geom_boxplot() +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, vjust=0, hjust=0.5),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.4, 'cm'),
    panel.grid = element_blank())
ggsave("met correlation.pdf", path= "./plots/", dpi=500)

# classification correct?
anno_met <- anno_UMC %>%
  filter(location != "primary") %>%
  filter(tumorType != "Mixed")

anno_met$pred_rf <- factor(anno_met$pred_rf, levels = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"))
anno_met$tumorType <- factor(anno_met$tumorType, levels = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"))

plt <- caret::confusionMatrix(anno_met$pred_rf,anno_met$tumorType)

plt <- as.data.frame(plt$table)

ggplot(plt, aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq), color = "white") +
  scale_fill_gradient(low = "#132B43", high = "#88CCEEFF") +
  #paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  labs(x = "True Class ",y = "Predicted Class (NN)") +
  scale_x_discrete(limits = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"),
                   breaks = c("PanNET", "PDAC", "ACC", "SPN", "NORM", "PB", "PanNEC"),
                   labels = c("PanNET", "PDAC", "ACC", "SPN", "NORMAL", "PB", "PanNEC")) +
  theme_bw(base_size = 16) +
  theme(axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
        axis.text.x = element_text(face = "bold", size = 9, vjust=0.1, hjust=0.5),
        axis.text.y = element_text(face = "bold", size = 8),
        axis.ticks=element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank())
ggsave("confusion matrix met.pdf", path= "./plots/", dpi=500)
