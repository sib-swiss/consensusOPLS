test_that("koplsModel", {
    kernel <- koplsKernel(X1 = demo_3_Omics[[1]], X2 = NULL, Ktype = 'p', params = c(order=1))
    expect_equal(kernel[1,2], -44.6770)
    expect_equal(kernel[6,3], 0.8963)
    expect_equal(kernel[9,11], -13.2810)
})
