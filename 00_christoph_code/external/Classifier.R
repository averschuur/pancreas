### Code Supplement: Machine learning analysis of DNA methylation profiles distinguishes primary lung squamous cell carcinomas from head and neck metastases

###  Authors Philipp Jurmeister, Michael Bockmayr, David Capper and Frederick Klauschen

### Disclaimer: The version provided here is only intended for review purposes. The use of the results for clinical practice is in the sole responsibility of the treating physician.

################################################
####################### Set path to directory
################################################

### Please adapt path 
setwd("path/to/directory/stm_code")


################################################
#######################  Load R Packages
################################################

### All the package can be obtained from CRAN or Biocoductor
### For installation of the package keras please refer to the Readme.html file

library(caret)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(keras)
library(kernlab)
library(minfi)
library(randomForest)
library(rmarkdown)


#################################################
########################  Define cases
#################################################

### Please define cases
### As an example, four cases from the TCGA LUSC cohort are included in this dataset. Of note, due to the requierments of the package minfi, cases from 450k arrays and EPIC arrays cannot be preprocessed simultanously. Further, a minimum of two samples is required at each run.

idatpaths <- c("./stm_code/idat/6042308087_R03C02",
               "./stm_code/idat/6055424084_R05C02",
               "./stm_code/idat/6285633083_R01C02",
               "./stm_code/idat/9305216133_R03C02")


#################################################
########################  Load classifiers
#################################################

##### Load neural networks
# Due to incompatibilites between different versions of the R-package Keras, we initialize all trained models from the weights only.

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

load_model_weights_hdf5(model1,'data/nn_fold1_weights.hdf5')
load_model_weights_hdf5(model2,'data/nn_fold2_weights.hdf5')
load_model_weights_hdf5(model3,'data/nn_fold3_weights.hdf5')
load_model_weights_hdf5(model4,'data/nn_fold4_weights.hdf5')
load_model_weights_hdf5(model5,'data/nn_fold5_weights.hdf5')   


predictNN <- function(X){
    p1 <- model1 %>% predict_proba(X)
    p2 <- model2 %>% predict_proba(X)
    p3 <- model3 %>% predict_proba(X)
    p4 <- model4 %>% predict_proba(X)
    p5 <- model5 %>% predict_proba(X)
    p = (p1+p2+p3+p4+p5) / 5
    rownames(p) <- rownames(X)
    colnames(p) <- c("HNSC","LUNGNORM","LUSC")
    p
}

##### Load Support Vector Machine classifier
load("data/svm_model.RData")

predictSVM <- function(X){
    p <- predict(svm.model,X,type="prob")
    rownames(p) <- rownames(X)
    colnames(p) <- c("HNSC","LUNGNORM","LUSC")
    p
}


##### Load Random Forest classifier
load("data/rf_model.RData")

predictRF <- function(X){
    p <- predict(rf.model,X,type="prob")
    rownames(p) <- rownames(X)
    colnames(p) <- c("HNSC","LUNGNORM","LUSC")
    p
}


#### Read idat files
rgset <- read.metharray(idatpaths)

#### Preprocess data with normal-exponential out-of-band normalization
mset <- preprocessNoob(rgset)

#### Extract beta values
bval <- getBeta(mset)

#### Impute missing values (if any) from nearest CpG and select sites used for classification

preprocessSCC <- function(bval){
  cat("Selecting relevant sites and imputing missing values \n")
  cpg.mnp <- as.character(unlist(read.table("data/cpg_mnp_428799.txt",header=FALSE)))
  cpg.scc <- as.character(unlist(read.table("data/cpg_scc_2000.txt",header=FALSE)))
  X <- data.frame(bval)[cpg.scc,]
  Y <- data.frame(bval)[cpg.mnp,]
  rownames(X) <- cpg.scc
  rownames(Y) <- cpg.mnp
  index <- 0
  narows <- which(apply(X,1,function(x) sum(is.na(x)) > 0))
  if(length(narows) > 0) print(c("index","missing","replaced_by","dist"))
  for(i in narows){
    xxx <- X[i,]
    cases.na <- which(is.na(xxx))
    cpg.na <- rownames(X)[i]
    cpg.na.ind <- which(rownames(Y) == cpg.na)
    cpg.na.pos <- as.numeric(gsub("cg","",cpg.na) )
    neighbors <- max(0,cpg.na.ind-10):min(nrow(Y),cpg.na.ind+10)
    neighpos <- as.numeric(gsub("cg","",rownames(Y)[neighbors]))
    neighdis <- abs(neighpos- cpg.na.pos)
    neighord <- order(neighdis)[-1]
    for(case.na in cases.na){
      index <- index +1
      neigh.na <- is.na(Y[neighbors,case.na])
      neigh.sel <- neighbors[setdiff(neighord,which(neigh.na))[1]]
      print(c(index,cpg.na,
              rownames(Y)[neigh.sel],
              cpg.na.pos-as.numeric(gsub("cg","", rownames(Y)[neigh.sel]))
      ))
      X[cpg.na,case.na] <- Y[neigh.sel,case.na]
    }
  }
  return(t(X))
}

dset <- preprocessSCC(bval)

#### Apply classifiers
table_nn <- predictNN(dset)
table_svm <- predictSVM(dset)
table_rf <- predictRF(dset)

########################################################################
################################ Render output file using markdown
########################################################################

render("markdown/Classifier_output.Rmd",output_dir=getwd())
