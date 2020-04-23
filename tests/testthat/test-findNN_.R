#R
context("findNN_") 


test_that("Test 1", {
    # the drawback of the findNN implementation. Use findNN_!
    expect_equal(findNN(3.5, 1:5), findNN(3.5, 1:6), tolerance=1.0)

    expect_equal(findNN_(3.5, 1:5), findNN_(3.5, 1:6), tolerance=0.0)
})

test_that("Test 2", {
    DB<-sort(rnorm(100, mean=100, sd=10))
    expect_equal(unique(DB[findNN(DB,DB)] - DB), 0, tolerance=0.0)
    expect_equal(unique(DB[findNN_(DB,DB)] - DB), 0, tolerance=0.0)
})

test_that("Test 3 -- testing lower and upper index", {
    DB <- seq(-1,1,length=101)
    query <- c(-1000,0,0.001,10,10000)
    result <- c(1, 51, 51, 101, 101)

    expect_equal(findNN(q=query,  vec=DB), result)
    expect_equal(findNN_(q=query,  vec=DB), result)
    expect_equal(findNN(q=query,  vec=DB), findNN_(q=query,  vec=DB))
})
