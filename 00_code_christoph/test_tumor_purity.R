# libraries
library(RFpurify)
library(minfi)
library(tidyverse)

# find idat files
test_data <- list.files(path = "../../epic-seq/analysis/data/idat/tech-val/",
                        pattern = "_Grn.idat", 
                        full.names = TRUE)

# laod data
test_data <- minfi::read.metharray(basenames = test_data) %>% 
  preprocessNoob %>% 
  getBeta

# get purity estimates
absolute <- RFpurify::predict_purity_betas(betas = test_data, method = "ABSOLUTE")
estimate <- RFpurify::predict_purity_betas(betas = test_data, method = "ESTIMATE")

# collect data into tibble (similiar to a dataframe)
data <- tibble(id = colnames(test_data),
               avg_beta = apply(test_data, 2, mean, na.rm = TRUE), 
               absolute = absolute, 
               estimate = estimate)

# plot 
data %>% 
  ggplot(aes(absolute, estimate)) + 
  geom_point(size = 3) + 
  theme_bw()