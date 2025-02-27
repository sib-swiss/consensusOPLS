test_that("koplsKernel", {
    kernel <- koplsKernel(X1 = demo_3_Omics[[1]], X2 = NULL, type = 'p', params = c(order=1))
    expect_equal(kernel[1,2], 175.6044, tolerance=1e-4)
    expect_equal(kernel[6,3], 228.6084, tolerance=1e-4)
    expect_equal(kernel[9,11], 145.23854, tolerance=1e-4)
    
    kernel <- koplsKernel(X1 = demo_3_Omics[[2]], X2 = NULL, type = 'p', params = c(order=0.5))
    expect_equal(kernel[1,2], 9085.007, tolerance=1e-4)
    expect_equal(kernel[6,3], 11383.501, tolerance=1e-4)
    expect_equal(kernel[9,11], 9438.85, tolerance=1e-4)
    
    kernel <- koplsKernel(X1 = demo_3_Omics[[3]], X2 = NULL, type = 'p', params = c(order=0.1))
    expect_equal(kernel[1,2], 8.560449, tolerance=1e-5)
    expect_equal(kernel[6,3], 9.715536, tolerance=1e-5)
    expect_equal(kernel[9,11], 9.892736, tolerance=1e-5)
})
