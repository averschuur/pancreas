library(keras)

model1 <- keras_model_sequential()
model1 %>%
  layer_dense(units = 64, activation = "tanh", input_shape = 2000) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 64, activation = "tanh") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 64, activation = "tanh") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 3, activation = 'softmax')

model2 <- keras_model_sequential()
model2 %>%
  layer_dense(units = 64, activation = "tanh", input_shape = 2000) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 64, activation = "tanh") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 64, activation = "tanh") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 3, activation = 'softmax')

model3 <- keras_model_sequential()
model3 %>%
  layer_dense(units = 64, activation = "tanh", input_shape = 2000) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 64, activation = "tanh") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 64, activation = "tanh") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 3, activation = 'softmax')

model4 <- keras_model_sequential()
model4 %>%
  layer_dense(units = 64, activation = "tanh", input_shape = 2000) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 64, activation = "tanh") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 64, activation = "tanh") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 3, activation = 'softmax')

model5 <- keras_model_sequential()
model5 %>%
  layer_dense(units = 64, activation = "tanh", input_shape = 2000) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 64, activation = "tanh") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 64, activation = "tanh") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 3, activation = 'softmax')

load_model_weights_hdf5(model1,'../epic-seq/analysis/input/jurmeister/models/nn_fold1_weights.hdf5')
load_model_weights_hdf5(model2,'../epic-seq/analysis/input/jurmeister/models/nn_fold2_weights.hdf5')
load_model_weights_hdf5(model3,'../epic-seq/analysis/input/jurmeister/models/nn_fold3_weights.hdf5')
load_model_weights_hdf5(model4,'../epic-seq/analysis/input/jurmeister/models/nn_fold4_weights.hdf5')
load_model_weights_hdf5(model5,'../epic-seq/analysis/input/jurmeister/models/nn_fold5_weights.hdf5')   


predictNN <- function(X, ...){
  p1 <- model1 %>% predict(X, ...)
  p2 <- model2 %>% predict(X, ...)
  p3 <- model3 %>% predict(X, ...)
  p4 <- model4 %>% predict(X, ...)
  p5 <- model5 %>% predict(X, ...)
  p = (p1+p2+p3+p4+p5) / 5
  rownames(p) <- rownames(X)
  colnames(p) <- c("HNSC","LUNGNORM","LUSC")
  p
}
