#count training set and test set

train_anno$set <- "train"
test_anno$set <- "test"
annobind <- rbind(train_anno, test_anno)

annobind %>%
  ggplot(aes(set))+
  geom_bar(aes(fill = tumorType)) +
  scale_fill_manual(values = branded_colors) +
  theme_bw(base_size = 18)
ggsave("Figure 2_distribution over train and test set.png", path= "./output/")
