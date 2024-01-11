test_that("koplsKernel", {
    VIP <- MBVIP(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")],
                 Y=demo_3_Omics$Y)
    expect_equal(VIP[[1]][17], c(MT15724=0.2592), tolerance=1e-3)
    expect_equal(VIP[[2]][6], c(GC26860=0.0692), tolerance=1e-3)
    expect_equal(VIP[[3]][21], c(MT1146=2.5998), tolerance=1e-3)
    # do not match 0.000* in VIP[[3]]
})
