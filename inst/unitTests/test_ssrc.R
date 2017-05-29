#R
# $HeadURL: https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/specL/inst/unitTests/test_ssrc.R $
# $Id: test_ssrc.R 103296 2015-04-30 06:08:56Z c.panse $
# $Date: 2015-04-30 08:08:56 +0200 (Thu, 30 Apr 2015) $

test_ssrc <-
function(){
        

	checkEqualsNumeric(sapply(c("SCHTAVGR", "SCHTGLGR", "EDLIAYLK"), ssrc),
		c(3.20805, 5.95145, 29.60045),
		tolerance = 0.01)

}


test_ssrc()
