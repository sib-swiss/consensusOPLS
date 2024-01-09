test_that("demo_3_Omics", {
    expect_equal(demo_3_Omics[[1]][1,2], c(MT15709=1.243))
    expect_equal(demo_3_Omics[[2]][2,4], c(GC26858=219.365316))
    expect_equal(demo_3_Omics[[3]][3,5], c(MT113=4.739317))
})
