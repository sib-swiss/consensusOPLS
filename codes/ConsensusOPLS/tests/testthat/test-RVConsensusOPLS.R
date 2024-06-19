test_that("RVConsensusOPLS", {
    rvcopls <- RVConsensusOPLS(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")],
                               Y=demo_3_Omics$Y,
                               maxOcomp = 10,
                               nfold=14)
    rvcopls.fac <- RVConsensusOPLS(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")],
                               Y=factor(c(rep('M',7), rep('F',7))),
                               nfold=3)
    rvcopls.reg <- RVConsensusOPLS(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")],
                                   Y=c(rnorm(7, mean=1, sd=0.01), rnorm(7, mean=0, sd=0.01)),
                                   modelType="reg",
                                   nfold=3)
    copls.reg <- ConsensusOPLS(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")],
                               Y=matrix(c(rnorm(7, mean=1, sd=0.01), rnorm(7, mean=0, sd=0.01))),
                               maxPcomp=1,
                               maxOcomp=3,
                               nperm=100,
                               modelType="reg",
                               mc.cores=1,
                               kernelParams = list(type = "p", params = c(order = 2)),
                               nfold=3, 
                               verbose=T)
    copls.da <- ConsensusOPLS(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")],
                              Y=demo_3_Omics$Y,
                              maxPcomp=1,
                              maxOcomp=3,
                              nperm=100,
                              modelType="da",
                              mc.cores=1,
                              nfold=14,
                              verbose=T)
    ## RV
    expect_equal(rvcopls$RV, c(MetaboData=0.770310488895171,
                               MicroData=0.835375809822633,
                               ProteoData=0.74971449858494), tolerance=1e-6)
    
    ## lambda
    expect_equal(rvcopls$Model$lambda[,1], 
                 c(MetaboData=0.0455910270422634,
                   MicroData=0.00474677906005336,
                   ProteoData=0.00357785429350788), tolerance=1e-6)
    expect_equal(rvcopls$Model$lambda[,2], 
                 c(MetaboData=0.236514988519443,
                   MicroData=0.0231772225562188,
                   ProteoData=0.153126721989076), tolerance=1e-6)
    
    ## blockContribution
    expect_equal(rvcopls$Model$blockContribution[,1], 
                 c(MetaboData=0.84559897268353,
                   MicroData=0.0880408220024505,
                   ProteoData=0.0663602053140197), tolerance=1e-6)
    expect_equal(rvcopls$Model$blockContribution[,2], 
                 c(MetaboData=0.572926698791582,
                   MicroData=0.0561437974371785,
                   ProteoData=0.37092950377124), tolerance=1e-6)
    
    ## normKernels
    expect_equal(rvcopls$normKernels[[1]][1,2],  c(0.0328819715926973), tolerance=1e-6)
    expect_equal(rvcopls$normKernels[[2]][7,7],  c(0.104945253831102),  tolerance=1e-6)
    expect_equal(rvcopls$normKernels[[3]][9,11], c(0.121529241109287),  tolerance=1e-6)
    
    #### Optimal model
    ## Kpreproc
    expect_equal(rvcopls$Model$K[1,1][[1]][9,11], c(0.006208603), tolerance=1e-6)
    expect_equal(rvcopls$Model$K[1,1][[1]][3,13], c(-0.01054343), tolerance=1e-6)
    
    ## Cp
    expect_equal(unname(rvcopls$Model$Cp[,1]),
                 c(-0.707106781186547,
                    0.707106781186548), tolerance = 1e-6)
    
    ## Sp
    expect_equal(rvcopls$Model$Sp[1,1], c(0.9574551557441), tolerance=1e-6)
    
    ## Up
    expect_equal(unique(rvcopls$Model$Up[,1]),
                c(-0.707106781186547, 0.707106781186547), tolerance=1e-6)
    
    ## params
    expect_equal(rvcopls$Model$params$nOcomp, c(1))
    
    ## Tp
    expect_equal(rvcopls$Model$Tp[[1]][1:3],
                 c(0.0942610932887749, -0.154506780676698, -0.0608792282254192), tolerance=1e-6)
    expect_equal(rvcopls$Model$Tp[[2]][4:6],
                 c(-0.00916472674616893, -0.208638510110777, -0.115250463686483), tolerance=1e-6)
    
    ## scores
    expect_equal(unname(rvcopls$Model$scores[6:8,'p_1']),
                 c(-0.115250463686483,-0.134547425243546,0.169584736593916), tolerance=1e-6)
    expect_equal(unname(rvcopls$Model$scores[6:8,'o_1']),
                 c(0.225639604453228, -0.750090705302492, -0.246527124390688), tolerance=1e-6)
                 
    ## loadings
    expect_equal(rvcopls$Model$loadings$MetaboData[1:3,'p_1'],
                 c(MT15708=-3.7904893970217, MT15709=0.477951498442441, MT15710=0.597346765073993), tolerance=1e-6)
    expect_equal(rvcopls$Model$loadings$MetaboData[1:3,'o_1'],
                 c(MT15708=-0.0906993946312754, MT15709=0.391637521782826, MT15710=-0.102030862736768), tolerance=1e-6)
    
    expect_equal(rvcopls$Model$loadings$MicroData[1:3,'p_1'],
                 c(GC26855=-121.851575745773, GC26856=132.585031777677, GC26857=75.1435270079714), tolerance=1e-6)
    expect_equal(rvcopls$Model$loadings$MicroData[1:3,'o_1'],
                 c(GC26855=-28.2239137515635, GC26856=-31.9262881875121, GC26857=1.89187255025149), tolerance=1e-6)
    
    expect_equal(rvcopls$Model$loadings$ProteoData[1:3,'p_1'],
                 c(MT1126=-8838.55937473373, MT1127=1072.45583338495, MT1128=-232.657268097895), tolerance=1e-6)
    expect_equal(rvcopls$Model$loadings$ProteoData[1:3,'o_1'],
                 c(MT1126=162.821016719371, MT1127=-41.4788108974671, MT1128=-18.1358103752954), tolerance=1e-6)
    
    #### CV
    ## AllYhat
    expect_equal(unname(rvcopls$cv$AllYhat[1:3,1]),
                 c(0.1381898, 0.6109852, 0.5683056), tolerance=1e-6)
    expect_equal(unname(rvcopls$cv$AllYhat[4:6,2]),
                 c(0.60339737, 0.23570761, 0.48595049), tolerance=1e-6)
    expect_equal(unname(rvcopls$cv$AllYhat[1:3,3]),
                 c(0.136238770, 0.524588637, 0.740291299), tolerance=1e-6)
    expect_equal(unname(rvcopls$cv$AllYhat[7:9,4]),
                 c(-0.2020121, 0.9947210, 0.4679297), tolerance=1e-6)
    
    ## Q2Yhat
    expect_equal(unname(rvcopls$cv$Q2Yhat[1:5]),
                 c(0.24144459290868, 0.254116314089248, 0.0688378757417817, -0.058455510451356, 0.159444129092519), tolerance=1e-6)
    
    ## DQ2Yhat
    expect_equal(unname(rvcopls$cv$DQ2Yhat[1:5]),
                 c(0.120292310258979, 0.151671627027632, -0.0418685758525927, -0.171662986447783, 0.0788120118584003), tolerance=1e-6)
    
    # ## Tcv
    # expect_equal(rvcopls$cv$Tcv[6:8,1],
    #              c(-3.05298886630587,-10.7977788716645, 2.7915511228586), tolerance=1e-6)
    # 
    # ## Yhat (at maxOcomp)
    # expect_equal(rvcopls$cv$Yhat[6:8,1],
    #              c(0.80176823634328,1.56729072041118,0.224073196497945), tolerance=1e-6)
    # expect_equal(rvcopls$cv$Yhat[6:8,2],
    #              c(0.19823176365672,-0.567290720411185,0.775926803502055), tolerance=1e-6)
    
    # expect_equal(rvcopls$Model$R2X,
    #              c(0.18253923588794,0.257726827483885))
    # expect_equal(rvcopls$Model$R2XO,
    #              c(0,0.0811103596697501))
    # expect_equal(rvcopls$Model$R2XC,
    #              c(0.18253923588794,0.176616467814135))
})
