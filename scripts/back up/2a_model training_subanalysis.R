##############################################################
### Evaluate NN model performance in TEST set and UMCU set ###
##############################################################

anno <- readRDS("./data/sample_annotation.rds")
betas <- readRDS(file = "./data/methylation_data_filtered.rds")
nn_model <- load_model_hdf5(filepath = "./output/model_nn_955.hdf5")
rf_model <- load_model_hdf5(filepath = "./output/model_rf_947.hdf5")
  
### create datasets, obtain labels ---------------------------------------------
## split test set in non UMCU testset 
test_anno <- anno %>% 
  anti_join(y = train_anno) %>%
  filter(source != "UMCU")

# stats
test_anno %>% 
  group_by(tumorType, source, location) %>%
  summarise(n = n())

## creat an UMCU validation set
val_anno <- anno %>% 
  anti_join(y = train_anno) %>%
  anti_join(y=test_anno)

# stats
val_anno %>% 
  group_by(tumorType, source, location) %>%
  summarise(n = n())


## obtain data and create labels for TEST set
# get test_anno data
test_data <- betas[top_var_probes, test_anno$arrayId]

# define target values, i.e. labels
test_labels <- test_anno$tumorType

# convert to one-hot encoding for neural network
test_labels_onehot <- to_one_hot(test_labels)



### obtain data and create labels for UMC set
# get val_anno data
val_data <- betas[top_var_probes, val_anno$arrayId]

# define target values, i.e. labels
val_labels <- val_anno$tumorType

# convert to one-hot encoding for neural network
val_labels_onehot <- to_one_hot(val_labels)


### evaluate nn model on test data  ---------------------------------------------------------


# predict
test_scores_nn <- predict(object = nn_model, t(test_data))

# add labels to score matrixs
test_labels_nn <- c("ACC", "NORM", "PanNEC", "PanNET", "PB", "PDAC", "SPN")
colnames(test_scores_nn) <- test_labels_nn

# calculate max score and label for all cases
test_max_score <- apply(test_scores_nn, 1, max)
test_classes <- apply(test_scores_nn, 1, function(x) test_labels_nn[which.max(x)])

# add data to annotation
test_anno <- test_anno %>% 
  mutate(nn_class = test_classes, 
         nn_max = test_max_score,
         nn_score_acc = test_scores_nn[, 1], 
         nn_score_norm = test_scores_nn[, 2],
         nn_score_nec = test_scores_nn[, 3],
         nn_score_net = test_scores_nn[, 4], 
         nn_score_pb = test_scores_nn[, 5], 
         nn_score_pdac = test_scores_nn[, 6], 
         nn_score_spn = test_scores_nn[, 7])


# merge levels of training and test dataset together and print the confusionMatrix
prediction <- test_anno$nn_class
actual <- test_anno$tumorType

u <- union(prediction, actual)
t <- table(factor(prediction, u), factor(actual, u))
cM_test <- confusionMatrix(t)

# Accuracy : 0.8484  


### evaluate nn model on UMCU data  ---------------------------------------------------------

# predict
val_scores_nn <- predict(object = nn_model, t(val_data))

# add labels to score matrixs
val_labels_nn <- c("ACC", "NORM", "PanNEC", "PanNET", "PB", "PDAC", "SPN")
colnames(val_scores_nn) <- val_labels_nn

# calculate max score and label for all cases
val_max_score <- apply(val_scores_nn, 1, max)
val_classes <- apply(val_scores_nn, 1, function(x) val_labels_nn[which.max(x)])

# add data to annotation
val_anno <- val_anno %>% 
  mutate(nn_class = val_classes, 
         nn_max = val_max_score,
         nn_score_acc = val_scores_nn[, 1], 
         nn_score_norm = val_scores_nn[, 2],
         nn_score_nec = val_scores_nn[, 3],
         nn_score_net = val_scores_nn[, 4], 
         nn_score_pb = val_scores_nn[, 5], 
         nn_score_pdac = val_scores_nn[, 6], 
         nn_score_spn = val_scores_nn[, 7])


# merge levels of training and test dataset together and print the confusionMatrix
prediction <- val_anno$nn_class
actual <- val_anno$tumorType

u <- union(prediction, actual)
t <- table(factor(prediction, u), factor(actual, u))
cM_UMCU <- confusionMatrix(t)

# Accuracy : 1

### apply CM on a list

my.list <- list(test_anno, val_anno)

CM <- lapply(my.list, function(x) confusionMatrix(table(x[,5], x[,7])))



for (i in 1:ncol(a[[i]]))
{
  a[[i]] <- confusionMatrix(a[,5], a[,7])
}

cM <- cbind(cM_test, cM_UMCU)

knitr::kable(tm,"latex",longtable =T,booktabs =T,caption ="Longtable")%>%
  add_header_above(c(" ","p=50%"=2,"p=25%"=2,"p=75%"=2))%>%
  kable_styling(latex_options =c("repeat_header"))



#################################################################
### Evaluate NN model performance in primaries and metastases ###
#################################################################

### create dataset of test and UMC set with primaries and mets and obtain data and labels ------------------------------

# use rest of the data set as test cohort
test_anno <- anno %>% 
  anti_join(y = train_anno)

# split test set in primary group and met group
test_anno_primary <- test_anno %>%
  filter(source != "UMCU") %>%
  filter(location %in% c("primary", "pancreas", "acc normal"))

# stats
test_anno_primary %>% 
  group_by(tumorType, source, location) %>%
  summarise(n = n())

# obtain data
test_primary_data <- betas[top_var_probes, test_anno_primary$arrayId]

# define target values, i.e. labels
test_primary_labels <- test_anno_primary$tumorType

# convert to one-hot encoding for neural network
test_primary_labels_onehot <- to_one_hot(test_primary_labels)

# use rest of the data set as test cohort mets
test_anno_mets <- test_anno %>% 
  anti_join(y = test_anno_primary) %>%
  filter(source != "UMCU")

test_anno_mets %>% 
  group_by(tumorType, source, location) %>%
  summarise(n = n())

# obtain data
test_mets_data <- betas[top_var_probes, test_anno_mets$arrayId]

# define target values, i.e. labels
test_mets_labels <- test_anno_mets$tumorType

# convert to one-hot encoding for neural network
test_mets_labels_onehot <- to_one_hot(test_mets_labels)

### create UMCU validation primaries set
val_anno_primary <- test_anno %>%
  filter(source == "UMCU") %>%
  filter(location %in% c("primary", "pancreas"))

val_anno_primary %>% 
  group_by(tumorType, source, location) %>%
  summarise(n = n())

# obtain data
val_primary_data <- betas[top_var_probes, val_anno_primary$arrayId]

# define target values, i.e. labels
val_primary_labels <- val_anno_primary$tumorType

# convert to one-hot encoding for neural network
val_primary_labels_onehot <- to_one_hot(val_primary_labels)

### create UMCU validation mets set
val_anno_mets <- test_anno %>% 
  anti_join(y = val_anno_primary) %>%
  filter(source == "UMCU")

val_anno_mets %>% 
  group_by(tumorType, source, location) %>%
  summarise(n = n())

# obtain data
val_mets_data <- betas[top_var_probes, val_anno_mets$arrayId]

# define target values, i.e. labels
val_mets_labels <- val_anno_mets$tumorType

# convert to one-hot encoding for neural network
val_mets_labels_onehot <- to_one_hot(val_mets_labels)



### evaluate nn model on test data primaries ---------------------------------------------------------

# predict
test_primary_scores_nn <- predict(object = nn_model, t(test_primary_data))

# add labels to score matrixs
colnames(test_primary_scores_nn) <- colnames(test_primary_labels_onehot)

test_primary_labels_nn <- c("ACC", "NORM", "PanNEC", "PanNET", "PB", "PDAC", "SPN")
colnames(test_primary_scores_nn) <- test_primary_labels_nn

# calculate max score and label for all cases
test_primary_max_score <- apply(test_primary_scores_nn, 1, max)
test_primary_classes <- apply(test_primary_scores_nn, 1, function(x) test_primary_labels_nn[which.max(x)])

# add data to annotation
test_anno_primary <- test_anno_primary %>% 
  mutate(nn_class = test_primary_classes, 
         nn_max = test_primary_max_score,
         nn_score_acc = test_primary_scores_nn[, 1], 
         nn_score_norm = test_primary_scores_nn[, 2],
         nn_score_nec = test_primary_scores_nn[, 3],
         nn_score_net = test_primary_scores_nn[, 4], 
         nn_score_pb = test_primary_scores_nn[, 5], 
         nn_score_pdac = test_primary_scores_nn[, 6], 
         nn_score_spn = test_primary_scores_nn[, 7])


# merge levels of training and test dataset together and print the confusionMatrix
prediction <- test_anno_primary$nn_class
actual <- test_anno_primary$tumorType

u <- union(prediction, actual)
t <- table(factor(prediction, u), factor(actual, u))
caret::confusionMatrix(t)

#Accuracy : 0.8386  

### evaluate nn model on test data mets -----------------------------------------

# predict
test_mets_scores_nn <- predict(object = nn_model, t(test_mets_data))

# add labels to score matrixs
test_mets_labels_nn <- c("ACC", "NORM", "PanNEC", "PanNET", "PB", "PDAC", "SPN")
colnames(test_mets_scores_nn) <- test_mets_labels_nn

# calculate max score and label for all cases
test_mets_max_score <- apply(test_mets_scores_nn, 1, max)
test_mets_classes <- apply(test_mets_scores_nn, 1, function(x) test_mets_labels_nn[which.max(x)])

# add data to annotation
test_anno_mets <- test_anno_mets %>% 
  mutate(nn_class = test_mets_classes, 
         nn_max = test_mets_max_score,
         nn_score_acc = test_mets_scores_nn[, 1], 
         nn_score_norm = test_mets_scores_nn[, 2],
         nn_score_nec = test_mets_scores_nn[, 3],
         nn_score_net = test_mets_scores_nn[, 4], 
         nn_score_pb = test_mets_scores_nn[, 5], 
         nn_score_pdac = test_mets_scores_nn[, 6], 
         nn_score_spn = test_mets_scores_nn[, 7])


# merge levels of training and test dataset together and print the confusionMatrix
prediction <- test_anno_mets$nn_class
actual <- test_anno_mets$tumorType

u <- union(prediction, actual)
t <- table(factor(prediction, u), factor(actual, u))
caret::confusionMatrix(t)

# Accuracy : 0.9524 


### evaluate nn model on UMCU data primaries ---------------------------------------------------------

model_probes <- top_var_probes_names

# predict
val_primary_scores_nn <- predict(object = nn_model, t(val_primary_data))

# add labels to score matrixs
val_primary_labels_nn <- c("ACC", "NORM", "PanNEC", "PanNET", "PB", "PDAC", "SPN")
colnames(val_primary_scores_nn) <- val_primary_labels_nn

# calculate max score and label for all cases
val_primary_max_score <- apply(val_primary_scores_nn, 1, max)
val_primary_classes <- apply(val_primary_scores_nn, 1, function(x) val_primary_labels_nn[which.max(x)])

# add data to annotation
val_anno_primary <- val_anno_primary %>% 
  mutate(nn_class = val_primary_classes, 
         nn_max = val_primary_max_score,
         nn_score_acc = val_primary_scores_nn[, 1], 
         nn_score_norm = val_primary_scores_nn[, 2],
         nn_score_nec = val_primary_scores_nn[, 3],
         nn_score_net = val_primary_scores_nn[, 4], 
         nn_score_pb = val_primary_scores_nn[, 5], 
         nn_score_pdac = val_primary_scores_nn[, 6], 
         nn_score_spn = val_primary_scores_nn[, 7])


# merge levels of training and test dataset together and print the confusionMatrix
prediction <- val_anno_primary$nn_class
actual <- val_anno_primary$tumorType

u <- union(prediction, actual)
t <- table(factor(prediction, u), factor(actual, u))
caret::confusionMatrix(t)

# Accuracy : 1   


### evaluate nn model on UMC data mets -----------------------------------------

# predict
val_mets_scores_nn <- predict(object = nn_model, t(val_mets_data))

# add labels to score matrixs
val_mets_labels_nn <- c("ACC", "NORM", "PanNEC", "PanNET", "PB", "PDAC", "SPN")
colnames(val_mets_scores_nn) <- val_mets_labels_nn

# calculate max score and label for all cases
val_mets_max_score <- apply(val_mets_scores_nn, 1, max)
val_mets_classes <- apply(val_mets_scores_nn, 1, function(x) val_mets_labels_nn[which.max(x)])

# add data to annotation
val_anno_mets <- val_anno_mets %>% 
  mutate(nn_class = val_mets_classes, 
         nn_max = val_mets_max_score,
         nn_score_acc = val_mets_scores_nn[, 1], 
         nn_score_norm = val_mets_scores_nn[, 2],
         nn_score_nec = val_mets_scores_nn[, 3],
         nn_score_net = val_mets_scores_nn[, 4], 
         nn_score_pb = val_mets_scores_nn[, 5], 
         nn_score_pdac = val_mets_scores_nn[, 6], 
         nn_score_spn = val_mets_scores_nn[, 7])


# merge levels of training and test dataset together and print the confusionMatrix
prediction <- val_anno_mets$nn_class
actual <- val_anno_mets$tumorType

u <- union(prediction, actual)
t <- table(factor(prediction, u), factor(actual, u))
caret::confusionMatrix(t)

# Accuracy : 1  




##############################################################
### Evaluate RF model performance in TEST set and UMCU set ###
##############################################################

### evaluate random forest on test data -------------------------------------------------------------
# scores
test_scores_rf <- predict(object = rf_model, t(test_data), type = "prob")

# classes
test_classes_rf <- predict(object = rf_model, t(test_data))

# calculate max score and label for all cases
test_max_score <- apply(test_scores_rf, 1, max)

# add data to annotation
test_anno <- test_anno %>% 
  mutate(rf_class = test_classes, 
         rf_max = test_max_score,
         rf_score_acc = test_scores_rf[, 1], 
         rf_score_norm = test_scores_rf[, 2],
         rf_score_nec = test_scores_rf[, 3],
         rf_score_net = test_scores_rf[, 4], 
         rf_score_pb = test_scores_rf[, 5], 
         rf_score_pdac = test_scores_rf[, 6], 
         rf_score_spn = test_scores_rf[, 7])

# merge levels of training and test dataset together and print the confusionMatrix
table(prediction = test_classes_rf, actual = test_anno$tumorType) %>% 
  caret::confusionMatrix()

# Accuracy : 0.9467  


### evaluate random forest on UMCU data -------------------------------------------------------------
# scores
val_scores_rf <- predict(object = rf_model, t(val_data), type = "prob")

# classes
val_classes_rf <- predict(object = rf_model, t(val_data))

# calculate max score and label for all cases
val_max_score <- apply(val_scores_rf, 1, max)

# add data to annotation
val_anno <- val_anno %>% 
  mutate(rf_class = val_classes_rf, 
         rf_max = val_max_score,
         rf_score_acc = val_scores_rf[, 1], 
         rf_score_norm = val_scores_rf[, 2],
         rf_score_nec = val_scores_rf[, 3],
         rf_score_net = val_scores_rf[, 4], 
         rf_score_pb = val_scores_rf[, 5], 
         rf_score_pdac = val_scores_rf[, 6], 
         rf_score_spn = val_scores_rf[, 7])

# merge levels of training and test dataset together and print the confusionMatrix
prediction <- val_classes_rf
actual <- val_anno$tumorType

u <- union(prediction, actual)
t <- table(factor(prediction, u), factor(actual, u))
cM_UMCU <- caret::confusionMatrix(t)

#Accuracy : 0.95 



#################################################################
### Evaluate RF model performance in primaries and metastases ###
#################################################################

### evaluate rf model on test data primaries ---------------------------------------------------------

# predict
test_primary_scores_rf <- predict(object = rf_model, t(test_primary_data), type = "prob")

# classes
test_primary_classes_rf <- predict(object = rf_model, t(test_primary_data))

# merge levels of training and test dataset together and print the confusionMatrix
table(prediction = test_primary_classes_rf, actual = test_anno_primary$tumorType) %>% 
  caret::confusionMatrix()

#Accuracy : 0.9462  


### evaluate rf model on test data mets ---------------------------------------------------------

# predict
test_mets_scores_rf <- predict(object = rf_model, t(test_mets_data), type = "prob")

# classes
test_mets_classes_rf <- predict(object = rf_model, t(test_mets_data))

# merge levels of training and test dataset together and print the confusionMatrix
prediction <- test_mets_classes_rf
actual <- test_anno_mets$tumorType

u <- union(prediction, actual)
t <- table(factor(prediction, u), factor(actual, u))
caret::confusionMatrix(t)

# Accuracy : 0.9524

### evaluate rf model on UMCU data primaries ---------------------------------------------------------

# predict
val_primary_scores_rf <- predict(object = rf_model, t(val_primary_data), type = "prob")

# classes
val_primary_classes_rf <- predict(object = rf_model, t(val_primary_data))

# merge levels of training and test dataset together and print the confusionMatrix
prediction <- val_primary_classes_rf
actual <- val_anno_primary$tumorType

u <- union(prediction, actual)
t <- table(factor(prediction, u), factor(actual, u))
caret::confusionMatrix(t)

#Accuracy : 0.9375 


### evaluate rf model on test data mets ---------------------------------------------------------

# predict
val_mets_scores_rf <- predict(object = rf_model, t(val_mets_data), type = "prob")

# classes
val_mets_classes_rf <- predict(object = rf_model, t(val_mets_data))

# merge levels of training and test dataset together and print the confusionMatrix
prediction <- val_mets_classes_rf
actual <- val_anno_mets$tumorType

u <- union(prediction, actual)
t <- table(factor(prediction, u), factor(actual, u))
caret::confusionMatrix(t)

#Accuracy : 1 


#############################################
### Compare model performance in TEST set ###
#############################################

### add performance to annotation ----------------------------------------------
test_anno <-  test_anno %>% 
  mutate(pred_nn_test = nn_class, 
         pred_rf_test = rf_class, 
         pred_nn_test_scores = apply(test_scores_nn, 1, max),
         pred_rf_test_scores = apply(test_scores_rf, 1, max))
         


test_anno %>% 
  mutate(nn = as.integer(tumorType == pred_nn_test), 
         rf = as.integer(tumorType == pred_rf_test)) %>% 
  select(nn, rf) %>% 
  pivot_longer(cols = everything()) %>% 
  group_by(name) %>% 
  summarise(accuracy = sum(value)/n()) %>% 
  ggplot(aes(name, accuracy, fill = name)) +
  geom_col(width = 0.7) +
  geom_text(aes(label=accuracy), vjust=1.6, color="white", angle = 90, size=3.5)+
  scale_fill_manual(values = branded_colors) +
  theme_bw(base_size = 18) +
  labs(x = NULL, y = "Accuracy (test cohort)")

#############################################
### Compare model performance in TEST set ###
#############################################

### add performance to annotation ----------------------------------------------
val_anno <-  val_anno %>% 
  mutate(pred_nn_val = nn_class, 
         pred_rf_val = rf_class, 
         pred_nn_val_scores = apply(val_scores_nn, 1, max),
         pred_rf_val_scores = apply(val_scores_rf, 1, max))



val_anno %>% 
  mutate(nn = as.integer(tumorType == pred_nn_val), 
         rf = as.integer(tumorType == pred_rf_val)) %>% 
  select(nn, rf) %>% 
  pivot_longer(cols = everything()) %>% 
  group_by(name) %>% 
  summarise(accuracy = sum(value)/n()) %>% 
  ggplot(aes(name, accuracy, fill = name)) +
  geom_col(width = 0.7) +
  geom_text(aes(label=accuracy), hjust=2, color="white", angle = 90, size=3.5)+
  scale_fill_manual(values = branded_colors) +
  theme_bw(base_size = 18) +
  labs(x = NULL, y = "Accuracy (test cohort)")
