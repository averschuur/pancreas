
# branded colors ---------------------------------------------------------------


branded_colors1 <- c("#e0607e", "#fea090", "#efdc60", "#9aa0a7", "#59c7eb", "#0a9086", "#3e5496", "#8e2043", "#007187")
branded_colors2 <- c("#f6d8ae", "#2e4057", "#da4167", "#3cdbd3", "#f4d35e", "#8d96a3")
branded_colors3 <- c("#ef476f", "#ffd166", "#06d6a0", "#118ab2", "#073b4c")




# detect idats
detect_idats <- function(dir){
  
  pat <- "_Grn.idat"
  
  # returns tibble with path, filename and array ID
  idats <- tibble(path = list.files(dir, pattern = pat, full.names = TRUE)) 
  idats <- idats %>% 
    mutate(filename = str_extract(string = path, 
                                  pattern = "[0-9]*_R[0-9]{2}C[0-9]{2}_Grn.idat")) %>% 
    mutate(array_id = str_replace(string = filename,
                                  pattern = pat, replacement = ""))
  return(idats)
}



# get control probes ------------------------------------------------
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



# save heatmap as pdf  ---------------------------------------------------------
save_pheatmap_pdf <- function(x, filename, width=30, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}




# miscellaneous ----------------------------------------------------------------

entropy <- function(x, na.rm = TRUE){
  if(na.rm){
    e <- -sum(x[!x == 0] * log2(x[!x == 0]))
  } else {
    e <- -sum(x * log2(x))
  }
}

to_cooccurrence <- function(vec){
  n <- length(vec)
  cmat <- matrix(data = NA, nrow = n, ncol = n)
  
  for(i in 1:n){
    cl <- vec[i]
    cmat[, i] <- ifelse(vec == cl, 1, 0) 
  }
  return(cmat)
}

