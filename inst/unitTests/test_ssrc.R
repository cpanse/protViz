#R
# $HeadURL: https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/specL/inst/unitTests/test_ssrc.R $
# $Id: test_ssrc.R 103296 2015-04-30 06:08:56Z c.panse $
# $Date: 2015-04-30 08:08:56 +0200 (Thu, 30 Apr 2015) $

test_ssrc <-
function(){
        
  # example of table iv [PMID:15238601]
	checkEqualsNumeric(sapply(c("SCHTAVGR", "SCHTGLGR", "EDLIAYLK"), ssrc),
		c(3.20805, 5.95145, 29.60045),
		tolerance = 0.01)

}


test_ssrc_irt <- function(){
  irtPeptide.ssrc <- c(7.33685 , 11.93745 , 29.22845 , 16.65045 , 20.58245 ,
                       27.03345 , 19.42845 , 22.32845 , 22.71845 , 36.67245 ,
                       25.31445 , 34.21845 , 38.108815 , 39.461215 , 40.516115 ,
                       40.488815 , 46.237915 , 41.876215 , 41.476515 , 46.025115 , 
                       7.33685 , 11.59745 , 16.65045 , 20.58245 , 19.42845 ,
                       22.32845 , 22.79845 , 25.57445 , 37.89545 , 40.488815 ,
                       46.111915)

  checkEqualsNumeric(sapply(as.character(iRTpeptides$peptide), ssrc),
                     irtPeptide.ssrc,
                     tolerance = 0.01)
                     
}

test_ssrc()
test_ssrc_irt()
