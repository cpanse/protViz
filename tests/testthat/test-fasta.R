#R

context("FASTA")

test_that("tryptic digest of P12763", {
  GT <- c('MK', 'SFVLLFCLAQLWGCHSIPLDPVAGYK', 'EPACDDPDTEQAALAAVDYINK',
            'HLPR', 'GYK', 'HTLNQIDSVK', 'VWPR',
            'RPTGEVYDIEIDTLETTCHVLDPTPLANCSVR', 'QQTQHAVEGDCDIHVLK',
            'QDGQFSVLFTK', 'CDSSPDSAEDVR', 'K', 'LCPDCPLLAPLNDSR',
            'VVHAVEVALATFNAESNGSYLQLVEISR', 'AQFVPLPVSVSVEFAVAATDCIAK',
            'EVVDPTK', 'CNLLAEK', 'QYGFCK', 'GSVIQK', 'ALGGEDVR',
            'VTCTLFQTQPVIPQPQPDGAEAEAPSAVPDAAGPTPSAAGPPVASVVVGPSVVAVPLPLHR',
            'AHYDLR', 'HTFSGVASVESSSGEAFHVGK', 'TPIVGQPSIPGGPVR', 'LCPGR', 'IR',
            'YFK', 'I')
  
  fname <- system.file("extdata", name='P12763.fasta', package = "protViz")
  TS <- Fasta$new(fname)

  TS$getTrypticPeptides()

  expect_equal(length(TS$getTrypticPeptides()), length(GT))
  expect_equal(sum(nchar(TS$getTrypticPeptides())), sum(nchar(GT)))
  expect_equal(sum(sapply(1:length(GT), 
                                function(i){
                                  GT[i] == TS$getTrypticPeptides()[i]
                                  })),
                         28)
})
