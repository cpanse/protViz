#R
# Christian Panse <cp@fgcz.ethz.ch> 
# 2017-11-21

as.data.frame.pepxml <- function(x, ...){
  
  do.call('rbind', 
          lapply(x$msms_run_summary, function(y){
            if("search_result" %in% names(y)){
              expect <- NA
              pepmodseq <- NA
              deltacn <- NA
              xcorr <- NA
              
              z <- y$search_result
              
              idx.expect <- which(sapply(z$search_hit, function(yy){yy[1] == "expect"}))
              idx.xcorr <- which(sapply(z$search_hit, function(yy){yy[1] == "xcorr"}))
              idx.deltacn <- which(sapply(z$search_hit, function(yy){yy[1] == "deltacn"}))
              
              if (length(idx.expect) == 1){
                expect <- as.numeric(z$search_hit[idx.expect]$search_score[2])
              }
              
              if (length(idx.xcorr) == 1){
                xcorr <- as.numeric(z$search_hit[idx.xcorr]$search_score[2])
              }
              
              if (length(idx.deltacn) == 1){
                deltacn <- as.numeric(z$search_hit[idx.deltacn]$search_score[2])
              }
              
              if("modification_info" %in% names(z$search_hit)){
                pepmodseq <- z$search_hit$modification_info$.attrs[1]
              }
              
              rv <- as.data.frame(t(c(y$.attr, z$search_hit$.attrs)))
              
              
              rv$pepmodseq <- pepmodseq
              rv$deltacn <- deltacn
              rv$xcorr <- xcorr
              rv$expect <- expect
              
              return(rv)
            }
            return (NULL)
          }))
  
}