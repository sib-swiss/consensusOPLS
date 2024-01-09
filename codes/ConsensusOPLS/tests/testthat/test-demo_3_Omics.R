test_that("demo_3_Omics", {
    expect_equal(demo_3_Omics[[1]][1,2], c(MT15709=1.78947304144416))
    expect_equal(demo_3_Omics[[2]][2,4], c(GC26858=0.433230880710958))
    expect_equal(demo_3_Omics[[3]][3,6], c(MT1130=1.00069814258038))
})
