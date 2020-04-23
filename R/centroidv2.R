.getpeakpos <- function(mz, intensity){
  peakpos <- which(diff(sign(diff(intensity))) == -2) + 1
  peakstart <- which(diff(sign(diff(intensity))) == 2) + 1
  peakstart <- rep(peakstart,2)
  peakstartsimple <- which(diff(sign(diff(intensity))) == 1) + 1
  peakstarts <- sort(c(peakstart, peakstartsimple))
  peakstartsM <- matrix(peakstarts, ncol = 2, byrow = T)
  colnames(peakstartsM) <- c("start", "end")
  peakidx <- cbind(peakpos = peakpos , start = peakstartsM[,"start"], end  = peakstartsM[,"end"])
  return(peakidx)
}


centroidv2 <- function(mZ, intensity,
                       .funpos = .getpeakpos,
                       .funtrapez = protViz:::.trapez,
                       .funapex = weighted.mean){
  peakidx <- .funpos(mZ, intensity)




  getCentroid <- function(pidx){
    i <- pidx["start"]:pidx["end"]
    intensity.auc <- .funtrapez(mZ[i], intensity[i])
    apexgroup <- (pidx["peakpos"] - 1):(pidx["peakpos"] + 1)
    mZ.centroid <- .funapex( mZ[apexgroup], intensity[apexgroup])
    return(c(mZ = mZ.centroid, intensity = intensity.auc))
  }

  res <- data.frame(t(apply(peakidx,1,getCentroid)))
  return(res)
}



