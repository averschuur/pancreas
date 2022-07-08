

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
