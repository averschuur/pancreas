### S_extra unsupervised_analysis

library(paletteer)




# estimate tumor type
anno %>% 
  ggplot(aes(tumorType, estimate, fill= tumorType, col = tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =2,alpha=0.7)+
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  labs(x = "Tumor Type", y = "Purity (ESTIMATE)") +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    panel.grid = element_blank())
ggsave("Figure S_tumortype estimate boxplot.pdf", path= "./plots/", dpi=500)

#stastiek
estimate <- anno[,c(5,8)]

boxplot(estimate~tumorType , data=estimate)

pg<-lm(estimate~tumorType , data=estimate)
anova(pg)

mfit2 <- aov(estimate~factor(tumorType), data=estimate)
TukeyHSD(mfit2)

#diff         lwr          upr     p adj

#PDAC-ACC      -0.2370668517 -0.27835689 -0.195776811 0.0000000
#SPN-ACC       -0.0717305889 -0.13117287 -0.012288309 0.0071839
#PDAC-NORM     -0.2402484392 -0.29690880 -0.183588078 0.0000000
#SPN-NORM      -0.0749121763 -0.14589748 -0.003926870 0.0310048
#PDAC-PanNEC   -0.2052563975 -0.27308055 -0.137432245 0.0000000
#PDAC-PanNET   -0.2044737219 -0.23670968 -0.172237759 0.0000000
#PDAC-PB       -0.1797991581 -0.24189282 -0.117705492 0.0000000
#SPN-PDAC       0.1653362628  0.10977102  0.220901506 0.0000000



# absolute tumor type
anno %>% 
  ggplot(aes(tumorType, absolute, fill= tumorType, col = tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =2,alpha=0.7)+
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  labs(x = "Tumor Type", y = "Purity (ABSOLUTE)") +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    panel.grid = element_blank())
ggsave("Figure S_tumortype absolute boxplot.pdf", path= "./plots/", dpi = 500)

#stastiek
absolute <- anno[,c(5,7)]

boxplot(absolute~tumorType , data=absolute)

pg<-lm(absolute~tumorType , data=absolute)
anova(pg)

mfit2 <- aov(absolute~factor(tumorType), data=absolute)
TukeyHSD(mfit2)
f <- TukeyHSD(mfit2, conf.level=0.95)
plot(f)
f %>%
  ggHSD() +
  theme_bw(base_size = 18) +
  theme(axis.text.y = element_text(face = "bold", size = 11, angle = 0, vjust=0.5, hjust=1))
ggsave("Figure S_tumortype absolute boxplot.pdf", path= "./plots/", dpi = 500)

#                      diff         lwr          upr     p adj
NORM-ACC      -0.060921297 -0.11962506 -0.002217536 0.0362014
PanNET-ACC    -0.100789534 -0.13820196 -0.063377112 0.0000000
PB-ACC        -0.065193097 -0.12886627 -0.001519928 0.0408406
PDAC-ACC      -0.311970256 -0.35205612 -0.271884395 0.0000000
SPN-ACC       -0.165772510 -0.22348122 -0.108063800 0.0000000
PDAC-NORM     -0.251048959 -0.30605688 -0.196041035 0.0000000
SPN-NORM      -0.104851213 -0.17376631 -0.035936116 0.0001800
PDAC-PanNEC   -0.245238441 -0.31108458 -0.179392306 0.0000000
SPN-PanNEC    -0.099040695 -0.17688112 -0.021200266 0.0035678
PDAC-PanNET   -0.211180722 -0.24247656 -0.179884886 0.0000000
SPN-PanNET    -0.064982976 -0.11697190 -0.012994052 0.0045509
PDAC-PB       -0.246777159 -0.30705993 -0.186494387 0.0000000
SPN-PB        -0.100579413 -0.17377389 -0.027384936 0.0011193
SPN-PDAC       0.146197746  0.09225300  0.200142489 0.0000000

library(MASS)

#make 4 sep groups
absolute$groups.7 <- with(absolute, paste(absolute,tumorType,sep = "."))

crab.aov <- aov(FL ~ groups.4, data = absolute)

crab.Tukey <- TukeyHSD(crab.aov)

plot(crab.Tukey)

# average methylation
anno %>% 
  ggplot(aes(tumorType, avg_beta_unfiltered, fill= tumorType, col = tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =2,alpha=0.7)+
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  labs(x = "Tumor Type", y = "Average methylation") +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    panel.grid = element_blank())
ggsave("Figure S_tumortype avg methylation boxplot.pdf", path= "./plots/", dpi = 500)

# conversion
anno %>% 
  ggplot(aes(tumorType, conversion, fill= tumorType, col = tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =2,alpha=0.7)+
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw(base_size = 18) +
  labs(x = "Tumor Type", y = "Conversion Score") +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    panel.grid = element_blank())
ggsave("Figure S_tumortype conversion score boxplot.pdf", path= "./plots/", dpi = 500)



### dimensionality reduction ---------------------------------------------------

# tissues
tissues <- read_csv(file = "./annotation/annotation_tissues.csv")
anno <- left_join(anno, tissues)

# load filtered beta values
betas <- readRDS("./input/betas_pancreas_filtered.rds")
betas <- betas[, anno$arrayId]
any(is.na(betas))

# add avg. methylation to annotation
anno <- anno %>% 
  add_column(avg_beta_filtered = apply(betas, 2, mean, na.rm = TRUE))

# get most variable probes across dataset and subset beta values
probes_topvar <- readRDS(file = "./output/pancreas_top_variable_probes_training_set.rds")

# pick betas for 5,000 top variable probes
betas_topvar <- betas[probes_topvar[1:5000], ]

# run UMAP
set.seed(45098)
umap_settings <- umap.defaults
umap_settings$n_neighbors = 15
umap_settings$min_dist = 0.2

umap <- umap(d = t(betas_topvar), config = umap_settings, ret_model = TRUE)

anno <- anno %>% 
  mutate(umap_x = umap$layout[, 1], 
         umap_y = umap$layout[, 2])


# plot UMAP
anno %>% 
  ggplot(aes(umap_x, umap_y, col = array)) +
  geom_point(size = 4) +
  #geom_text(aes(label = sampleName), size = 4) +
  theme_classic(base_size = 24) +
  scale_color_manual(values = branded_colors1) +
  theme(legend.direction = "horizontal", legend.position = "bottom") +
  labs(x = "UMAP 1", y = "UMAP 2")

anno %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType, shape = array)) + 
  geom_point(aes(color = tumorType), size = 4) +
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
    legend.margin=margin(0,0,0,0),
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, 'cm'),
    panel.grid = element_blank()) +
  guides(col = guide_legend(nrow = 1))
ggsave("Figure S_umap_tissues.pdf", path= "./plots/", dpi=500)

