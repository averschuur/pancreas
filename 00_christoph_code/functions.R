

getControlBeta <- function(rgSet, 
                           controls = "BISULFITE CONVERSION I",
                                        sampNames = NULL) {
  
  if (!controls %in% c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II")) {
    stop("controls must be BISULFITE CONVERSION I or BISULFITE CONVERSION II")
    }

  r <- getRed(rgSet)
  g <- getGreen(rgSet)
  

  ctrlAddress <- minfi::getControlAddress(rgSet, controlType = controls)
  
  # Red channel
  ctlWide <- as.matrix(log2(r[ctrlAddress, ,drop = FALSE]))
  if (!is.null(sampNames)) colnames(ctlWide) <- sampNames
  ctlR <- reshape::melt(ctlWide, varnames = c("address", "arrayId"))
  
  # Green channel
  ctlWide <- as.matrix(log2(g[ctrlAddress, ,drop = FALSE]))
  if (!is.null(sampNames)) colnames(ctlWide) <- sampNames
  ctlG <- reshape::melt(ctlWide, varnames = c("address", "arrayId"))
  
  # Plot
  ctl <- rbind(
    cbind(channel = "Red", ctlR),
    cbind(channel = "Green", ctlG))
  return(ctl)
}



# model training and evaluation ------------------------------------------------

to_one_hot <- function(labels){
  unique_labels <- sort(unique(labels))
  n_labels <- length(unique_labels)
  n_samples <- length(labels)
  
  onehot <- matrix(data = 0, ncol = n_labels, nrow = n_samples)
  
  for (i in 1:n_labels){
    onehot[labels == unique_labels[i], i] <- 1
  }
  
  colnames(onehot) <- unique_labels
  return(onehot)
}


slide_along_cutoff <- function(label_real, label_pred, scores, cutoffs){
  
  # instantiate empty variables
  predictable <- NULL
  accuracy <- NULL
  
  n <- length(cutoffs)
  n_samples <- length(label_real)
  
  for (i in 1:n){
    which_na <- which(scores < cutoffs[i])
    n_na <- length(which_na)
    label_pred[which_na] <- NA
    accuracy <- c(accuracy, sum(label_real == label_pred, na.rm = TRUE) / (n_samples - n_na))
    predictable <- c(predictable, (n_samples - n_na)/n_samples)
  }
  
  results <- tibble(cutoff = cutoffs, 
                    accuracy = accuracy, 
                    predictable = predictable) 
  return(results)
}
