test_that("demo_3_Omics", {
    expect_equal(demo_3_Omics[[1]][1,2], c(1.243))
    expect_equal(demo_3_Omics[[2]][2,4], c(219.365316))
    expect_equal(demo_3_Omics[[3]][3,5], c(4.739317))
})
