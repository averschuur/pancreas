# last edited 04/01/2023 by AV Verschuur

### Load required packages and sources  ----------------------------------------------

library(tidyverse)
library(minfi)
library(paletteer)
library(Rtsne)

source("./scripts/functions.R")


#### load data and annotation 
arrayRawEpic <- readRDS(file = "./data/preprocessing/arrayRawEpic_ChanJakelSelenicaDiDomenico.rds")
arrayRaw450k <- readRDS(file = "./data/preprocessing/arrayRaw450k_ChanJakelSelenicaDiDomenico.rds")
arrayEndo <- readRDS(file = "./data/preprocessing/arrayRawEndo.rds")
arrayYachida <- readRDS(file = "./data/preprocessing/arrayRawYachida.rds")
arrayBenhamida <- readRDS(file = "./data/preprocessing/arrayRawBenhamida.rds")

# investigate QC measures 1: bisulfite conversion efficiency 
conv_RawEpic <- getControlBeta(arrayRawEpic, controls = "BISULFITE CONVERSION II") %>% 
  as_tibble() %>% 
  group_by(arrayId) %>% 
  summarise(conversion = min(value[channel == "Red"]) / max(value[channel == "Green"]))

rm(arrayRawEpic)

conv_Raw450k <- getControlBeta(arrayRaw450k, controls = "BISULFITE CONVERSION II") %>% 
  as_tibble() %>% 
  group_by(arrayId) %>% 
  summarise(conversion = min(value[channel == "Red"]) / max(value[channel == "Green"]))

rm(arrayRaw450k)

conv_Benhamida <- getControlBeta(arrayBenhamida, controls = "BISULFITE CONVERSION II") %>% 
  as_tibble() %>% 
  group_by(arrayId) %>% 
  summarise(conversion = min(value[channel == "Red"]) / max(value[channel == "Green"]))

rm(arrayBenhamida)

conv_Endo <- getControlBeta(arrayEndo, controls = "BISULFITE CONVERSION II") %>% 
  as_tibble() %>% 
  group_by(arrayId) %>% 
  summarise(conversion = min(value[channel == "Red"]) / max(value[channel == "Green"]))

rm(arrayEndo)

conv_Yachida <- getControlBeta(arrayYachida, controls = "BISULFITE CONVERSION II") %>% 
  as_tibble() %>% 
  group_by(arrayId) %>% 
  summarise(conversion = min(value[channel == "Red"]) / max(value[channel == "Green"]))

rm(arrayYachida)

## merge conversions:
conversion <- full_join(conv_Benhamida, conv_Endo)
conversion <- full_join(conversion, conv_Raw450k)
conversion <- full_join(conversion, conv_RawEpic)
conversion <- full_join(conversion, conv_Yachida)

conversion[,1] <- as.matrix(conversion[,1]) %>%
  str_extract(pattern = "[0-9]*_R[0-9]{2}C[0-9]{2}$")

#load annotation
anno <- readRDS("./data/sample_annotation_extended.rds")


# merge conversion scores with annotation
anno <- anno %>% 
  left_join(conversion)
saveRDS(object = anno, file= "./data/sample_annotation_conversion scores.rds")

# save as csv
anno1 <- readRDS("./data/sample_annotation.rds")
anno1 <- anno1 %>% 
  left_join(conversion)

# add classification of conversion tot anno
anno <- anno %>%
  mutate(good_score1 = case_when(conversion >= 1 ~ 'good',
                                 conversion >= 0.8 ~ 'medium',
                                 conversion < 0.8 ~ 'bad'))

# stats
anno %>%
  group_by(good_score1) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)*100)

count <- anno %>%
  group_by(source, good_score1) %>%
  count(good_score1)

count %>%
  ggplot(aes(fill = good_score1, n, source)) +
  geom_bar(position="fill", stat="identity")
ggsave("Figure s1_Control Scores by source.pdf", path= "./output/", dpi=500)

# save anno
write.csv(anno, file= "./annotation/sample_annotation_conversion scores.csv")

# save anno ex bad and medium samples
anno_ebs <- anno %>%
  subset(good_score1 == "good")

saveRDS(anno_ebs, file ="./data/sample_annotation_extended_ex_bad_samples.rds")

# Prepare for data visualisation
anno %>%
  mutate(reply = anno$arrayId %in% colnames(betas)) %>%
  group_by(reply, sampleName) %>%
  summarize(reply_sum = sum(reply))

anno <- anno %>% 
  filter(!sampleName == "UMCU_ACC2")

anno1 <- anno %>% 
  filter(tumorType %in% c("PanNET", "ACC", "SPN", "PDAC", "normal", "acc normal", "PanNEC", "PB")) %>% 
  filter(location %in% c("primary", "pancreas")) %>% 
  filter(source != "UMCU")



### create UMAP:

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


# Figure 1S: UMAP conversion
anno1 %>% 
  ggplot(aes(umap_x, umap_y, col = conversion)) + 
  geom_point(shape=19, size = 5) +
  paletteer::scale_color_paletteer_c("ggthemes::Green-Blue-White Diverging") +
  labs(#title="UMAP Clustering",
    x = "UMAP 1", 
    y = "UMAP 2",
    col = "ControlScores") +
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
ggsave("Figure S1_UMAP-conversion.pdf", path= "./output/", dpi=500)

anno1 %>% 
  ggplot(aes(tumorType, conversion, fill=tumorType, col=tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =2) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN")) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe", breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN")) +
  labs(#title="Absolute Tumor Purity By Tumor Type",
    x = "", 
    y = "Control Scores") +
  geom_hline(yintercept=0.8, color = "grey", linetype = "dashed") +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("Figure s1_tumorType_Control Scores.pdf", path= "./output/", dpi=500)

anno1 %>% 
  ggplot(aes(source, conversion, fill=source, col=source)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=source), position=position_jitter(0.2),shape = 21,size =2) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", 
                                    #breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN")
                                    ) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe", 
                                     #breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN")
                                     ) +
  labs(#title="Absolute Tumor Purity By Tumor Type",
    x = "", 
    y = "Control Scores") +
  geom_hline(yintercept=0.8, color = "grey", linetype = "dashed") +
  #scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   #breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   #labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("Figure s1_Control Scores by study.pdf", path= "./output/", dpi=500)
