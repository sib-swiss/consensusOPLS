test_that("VIP", {
    VIP <- VIP(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")],
               Y=demo_3_Omics$Y)
    
    ## MetaboData
    expect_equal(VIP[[1]][11:13,1], 
                 unname(c(MT15718=0.764957797016524,MT15719=0.00856865416329684,MT15720=3.04645000211245)), tolerance=1e-6)
     
    ## MicroData   
    expect_equal(VIP[[2]][7:10,1], 
                 unname(c(GC26861=0.0214735865930633, GC26862=0.00142403129691162, 
                          GC26863=0.00345724527441899, GC26864=0.0253278644980716)), tolerance=1e-6)
    
    ## ProteoData
    expect_equal(VIP[[3]][53:58,1], 
                 unname(c(MT1303=0.00022579624142751, MT1305=8.24879827177417, MT1307=3.68609245397797, 
                          MT1308=2.41876426797679, MT1351=0.000103541922864548, MT1368=2.49498609312164e-05)), tolerance=1e-6)
})
