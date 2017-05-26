#R
# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/lcmsmap.R $
# $Id: lcmsmap.R 6178 2014-02-27 09:33:30Z cpanse $
# $Date: 2014-02-27 10:33:30 +0100 (Thu, 27 Feb 2014) $


lcmsmap <- function(data, charges=2:3, score.cutoff = 30, ...){
  
  s.rtinseconds <- NULL
  s.pepmass <- NULL 
  s.intensity <- NULL
  s.score <- NULL
  s.charge <- NULL
  
    if (is.mascot(data)){
      s.rtinseconds <- as.numeric(unlist(lapply(data$queries, function(x){x$RTINSECONDS}))) 
      s.pepmass <- as.numeric(unlist(lapply(data$queries, function(x){x$query_moverz})))
      s.intensity <- as.numeric(lapply(data$queries, function(x){x$TotalIonsIntensity}))
      s.score <- .mascot.get.pep_score(data)
      s.charge <- as.integer(gsub("[+]", "", unlist(lapply(data$queries, function(x){(x$query_charge)}))))
    }else{
      s.rtinseconds <- lapply(data, function(x){return (x$rtinseconds)})
      s.pepmass <- lapply(data, function(x){return (x$pepmass)})
      s.intensity <- lapply(data, function(x){return (sum(x$intensity))})
      s.score <- lapply(data, function(x){return (x$mascotScore)})
      s.charge <- as.double(lapply(data, function(x){return (x$charge)}))
    }

  S <- data.frame(rt = s.rtinseconds,
                  pepmass = s.pepmass,
                  intensity = s.intensity,
                  score = s.score,
                  charge = s.charge)
  
    for (c in charges){
        plot(S$rt, S$pepmass,
            type='n',
            main=paste('LC-MS overview [', c, ' +]', sep=''),
            xlab='rt [seconds]',
            ylab='pepmass [m/Z]', ...)
      
        S.f <- S[S$charge == c & S$score < score.cutoff, ]
      
        points(S.f$rt, S.f$pepmass, 
            pch=16, col=rgb(0.1,0.1,0.1, alpha=0.1), cex=0.75)
        
        
        S.f <- S[S$charge == c & S$score >= score.cutoff, ] 
        points(S.f$rt, S.f$pepmass, 
            pch=16, col=rgb(0.1,0.1,0.8, alpha=0.4), cex=1)

        legend("topleft", c(paste('score.cutoff = ', score.cutoff), 'rest'), 
               lwd=c(4), 
               col=c(rgb(0.1,0.1,0.8, alpha=0.9), 
                     'grey')) 
    }
}
