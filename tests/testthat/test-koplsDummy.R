test_that("koplsDummy", {
    expect_equal(ConsensusOPLS:::koplsDummy(c(rep(1, 7), rep(2, 7))), 
                 demo_3_Omics$Y, 
                 ignore_attr = TRUE)
})
test_that("koplsReDummy", {
    expect_equal(ConsensusOPLS:::koplsReDummy(matrix(demo_3_Omics$Y, byrow=F, ncol=2, dimnames=list(NULL, c(1,2)))),
                 c(rep('1', 7), rep('2', 7)))
})
