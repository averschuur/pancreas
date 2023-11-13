##############################
### FIGURE 5: Case Reports ###
##############################

library(tidyverse)
library(minfi)
library(doParallel)

### prepare data
#idats <- list.files(path = "./data/UMCU/Idat Files_All/",
                    recursive = TRUE,
                    full.names = TRUE, 
                    pattern = "_Grn.idat")

#raw <- read.metharray(basenames = idats, force = TRUE)

# preprocess data
#preprocessed_epic <- raw %>%
#  preprocessNoob(dyeMethod= "single")

#betas <- getBeta(preprocessed_epic)

# save unfiltered betas
#saveRDS(object = betas, file = "./input/betas_UMCU_unfiltered.rds")
betas_UMC <- readRDS("./input/betas_UMCU_unfiltered.rds")

# create annotation
#anno_files <- list.files(path = "./annotation/",
#                         pattern = ".csv",
#                         full.names = TRUE)
#anno_files <- anno_files[grepl(x = anno_files, pattern = "*annotation_umcu")]
#anno_files <- anno_files[!grepl(x = anno_files, pattern = "*annotation_umcu_paired_samples")]
#anno_UMC <- lapply(as.list(anno_files), read_csv)
#anno_UMC <- Reduce(f = bind_rows, x = anno_UMC)


#saveRDS(object = anno_UMC, file = "./output/anno_UMCU.rds")

# open data
anno_UMC <- readRDS("./output/anno_UMCU.rds")
betas_UMC <- readRDS(file = "./input/betas_UMCU_unfiltered.rds")
betas_UMC <- betas_UMC[, anno_UMC$arrayId]

rf_model <- readRDS(file = "./output/rf_model_default_01112023.rds")
od_model <- readRDS("./output/outlier_det_model_01112023.rds")

# rf predictions
rf_pred_class <- predict(rf_model, newdata = t(betas_UMC))
rf_pred_scores <- predict(rf_model, newdata = t(betas_UMC), type = "prob")
rf_pred_scores_max <- apply(rf_pred_scores, 1, max)
rf_pred_class <- apply(rf_pred_scores, 1, function(x) colnames(rf_pred_scores)[which.max(x)]) %>% as.factor

# add performance to annotation
anno_UMC <- anno_UMC %>% 
  mutate(pred_rf = rf_pred_class, 
         pred_scores_rf = apply(rf_pred_scores, 1, max))

# calculate score entropy
rf_pred_scores <- rf_pred_scores %>% 
  mutate(entropy = apply(rf_pred_scores[, 1:7], 1, entropy))
anno_UMC <- cbind(anno_UMC, rf_pred_scores$entropy)
colnames(anno_UMC) <- c("arrayId","source","sampleName","arrayType","tumorType","location","winning_class","winning_score","entropy")

#outlier prediction
anno_UMC <- anno_UMC %>%
  mutate(od_score = predict(object = od_model, newdata = ., type = "response"), 
         od_class = ifelse(od_score > 0.5, "pancreas", "outlier"))


### extract classification for seperate samples
test <- rf_pred_scores %>% 
  as_tibble(rownames = "arrayId") %>% 
  pivot_longer(cols = -c(arrayId,entropy), names_to = "tumor_type", values_to = "rf_score") %>% 
  mutate(rf_score = as.numeric(rf_score))

### Extract data for UMCU_SPN1 = 205555390010_R06C01
UMCU_SPN1 <- test[test$arrayId == "205555390010_R06C01",]

UMCU_SPN1 %>%
  ggplot(aes(tumor_type, rf_score, fill = tumor_type)) +
  geom_col() +
  theme_bw(base_size = 18) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  theme(panel.border = element_blank(),
        panel.grid = element_line(color = "grey",
                                  size = 0.5),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_blank(),
        axis.line = element_blank()) +
  labs(x = NULL, y = NULL) +
  coord_polar()
ggsave("Figure 5_UMCU_SPN1_08112023.pdf", path= "./plots/", dpi=500)

### Extract data for UMCU_ACC1 = 205555390010_R07C01
UMCU_ACC1 <- test[test$arrayId == "205555390010_R07C01",]

UMCU_ACC1 %>%
  ggplot(aes(tumor_type, rf_score, fill = tumor_type)) +
  geom_col() +
  theme_bw(base_size = 18) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  theme(panel.border = element_blank(),
        panel.grid = element_line(color = "grey",
                                  size = 0.5),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_blank(),
        axis.line = element_blank()) +
  labs(x = NULL, y = NULL) +
  coord_polar()
ggsave("Figure 5_UMCU_ACC1_08112023.pdf", path= "./plots/", dpi=500)


### Extract data for UMCU_ACC1_m = 205555390025_R01C01
UMCU_ACC1_m <- test[test$arrayId == "205555390025_R01C01",]

UMCU_ACC1_m %>%
  ggplot(aes(tumor_type, rf_score, fill = tumor_type)) +
  geom_col() +
  theme_bw(base_size = 18) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  theme(panel.border = element_blank(),
        panel.grid = element_line(color = "grey",
                                  size = 0.5),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_blank(),
        axis.line = element_blank()) +
  labs(x = NULL, y = NULL) +
  coord_polar()
ggsave("Figure 5_UMCU_ACC1_m_08112023.pdf", path= "./plots/", dpi=500)

### Extract data for UMCU_ACC2 = 206601450125_R05C01
UMCU_ACC2 <- test[test$arrayId == "206601450125_R05C01",]

UMCU_ACC2 %>%
  ggplot(aes(tumor_type, rf_score, fill = tumor_type)) +
  geom_col() +
  theme_bw(base_size = 18) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe") +
  theme(panel.border = element_blank(),
        panel.grid = element_line(color = "grey",
                                  size = 0.5),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_blank(),
        axis.line = element_blank()) +
  labs(x = NULL, y = NULL) +
  coord_polar()
ggsave("Figure 5_UMCU_ACC2_08112023.pdf", path= "./plots/", dpi=500)
