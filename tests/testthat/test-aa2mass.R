#R

context("amino acid to mass")
test_that("MELIVSK", {

    expect_equal(aa2mass("MELIVSK")[[1]], 
    c(131.04049, 129.04259, 113.08406, 113.08406, 99.06841, 87.03203, 128.09496),
    tolerance=0.001)
})


test_that("aa2mass verus parentIonMass", {

    # just a test
    peptides<-c('HTLNQIDSVK', 'ALGGEDVR', 'TPIVGQPSIPGGPVR')
    C_term <- 17.002740
    N_term <- 1.007825
    H_ <- 1.008

    expect_equal(parentIonMass(peptides),
        unlist(lapply(aa2mass(peptides), sum)) + C_term + N_term + H_, 
        tolerance=0.001)

})

