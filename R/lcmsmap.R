#R
# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/lcmsmap.R $
# $Id: lcmsmap.R 6178 2014-02-27 09:33:30Z cpanse $
# $Date: 2014-02-27 10:33:30 +0100 (Thu, 27 Feb 2014) $


as.data.frame.psmSet <- function(x, ...){
  data <- x
  rtinseconds <- as.numeric(unlist(lapply(data, function(x){return (x$rtinseconds)})))
  pepmass <- as.numeric(unlist(lapply(data, function(x){return (x$pepmass)})))
  intensity <- as.numeric(unlist(lapply(data, function(x){return (sum(x$intensity))})))
  score <- as.numeric(unlist(lapply(data, function(x){return (x$mascotScore)})))
  charge <- as.integer(unlist((lapply(data, function(x){return (x$charge)}))))
  
  data.frame(RTINSECONDS = rtinseconds,
                  moverz = pepmass,
                  intensity = intensity,
                  score = score,
                  charge = charge, ...)
}

lcmsmap <- function(data, charges = 2:3, score.cutoff = 30, ...){
  
  if (is.mascot(data)){
    S <- as.data.frame.mascot(data)
  }else if(is.psmSet(data)){
    S <- as.data.frame.psmSet(data)
  }else{return(NULL)}
  
  cm <- rev(rainbow(max(charges), alpha=0.3))
  
  plot(S$RTINSECONDS, S$moverz,
       type = 'n',
       main = 'LC-MS2 overview',
       xlab = 'rt [seconds]',
       ylab = 'pepmass [m/Z]', ...)
  
  legend("topleft", paste(charges, "+", sep=''),
         title = paste('score.cutoff >= ', score.cutoff), 
         pch = 16,
         col = cm[charges]) 
  
  
  S.f <- S[!(S$charge %in% charges) & S$score < score.cutoff, ]
  points(S.f$RTINSECONDS, S.f$moverz, 
         pch = 16, 
         col = rgb(0.1,0.1,0.1, alpha=0.05), 
         cex = 0.75)
  
  S.f <- S[S$charge %in% charges & S$score >= score.cutoff, ] 
  points(S.f$RTINSECONDS, S.f$moverz, 
         pch = 16,
         col = cm[S.f$charge],
         cex = 1.0)
}
