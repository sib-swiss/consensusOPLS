test_that("koplsModel", {
    rvcopls.dummy <- RVConsensusOPLS(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")],
                                     Y=demo_3_Omics$Y,
                                     nrcv=14)
    rvcopls <- RVConsensusOPLS(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")],
                               Y=factor(c(rep('M',7), rep('F',7))),
                               nrcv=3)
    rvcopls.reg <- RVConsensusOPLS(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")],
                                   Y=c(rnorm(7, mean=1, sd=0.01), rnorm(7, mean=0, sd=0.01)),
                                   modelType="reg",
                                   nrcv=3)
    rvcopls.perm <- RVConsensusOPLSPerm(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")],
                                   Y=matrix(c(rnorm(7, mean=1, sd=0.01), rnorm(7, mean=0, sd=0.01))),
                                   PredLVs=1,
                                   maxOrtholvs=3,
                                   nbruns=10,
                                   modelType="reg",
                                   mc.cores=1,
                                   nrcv=3)
    expect_equal(rvcopls$Model$Cp[,1],
                 c(F=-0.707106781186547,
                   M=0.707106781186548))
    expect_equal(rvcopls$RV, c(MetaboData=0.770310488895171,
                               MicroData=0.835375809822633,
                               ProteoData=0.74971449858494))
    expect_equal(rvcopls$AMat[[1]][1,2], c(0.0328819715926973))
    expect_equal(rvcopls$AMat[[2]][7,7], c(0.104945253831102))
    expect_equal(rvcopls$AMat[[3]][9,11], c(0.121529241109287))
    #expect_equal(rvcopls$Model$Sp,
    #             1.12130867444689, tolerance=1e-4)
    #expect_equal(unique(rvcopls$Model$Up[,1]),
    #             c(-0.76149961050859, 0.652713951864506), tolerance=1e-4)
    # expect_equal(rvcopls$Model$Tp,
    #              c(-0.299320882028158,-0.234737503178361,-0.239561232530216,
    #                -0.329485836706056,-0.247801779546492,-0.414865384358841,
    #                0.201072327452167,0.245312729047383,0.232741457001741,
    #                0.239571529778011,0.277280779035913,0.23723617065337,
    #                0.332557625379539,-0.293162869190637,-0.275152737076048,
    #                -0.295631263365155,-0.3005047448671,-0.302314707969838,
    #                -0.299006295879345,0.240391760356673,0.252064460741037,
    #                0.247904497175615,0.243334708084674,0.293420333759489,
    #                0.246726860097512,0.241929998133122))
    # expect_equal(rvcopls$Model$T,
    #              c(-0.293162869190637,-0.275152737076048,-0.295631263365155,
    #                -0.3005047448671,-0.302314707969838,-0.299006295879345,
    #                0.240391760356673,0.252064460741037,0.247904497175615,
    #                0.243334708084674,0.293420333759489,0.246726860097512,
    #                0.241929998133122))
    # expect_equal(rvcopls$Model$To,
    #              c(-0.0341849943797184,0.224357204198932,0.311261723470934,
    #                -0.160882818495652,0.302617776374132,-0.643168891168641,
    #                -0.218274080993194,-0.0374809075737232,-0.0841746285368153,
    #                -0.0208905425593303,-0.0895955565657284,-0.05268569161318,
    #                0.503101407841982))
    # expect_equal(rvcopls$Model$R2X,
    #              c(0.18253923588794,0.257726827483885))
    # expect_equal(rvcopls$Model$R2XO,
    #              c(0,0.0811103596697501))
    # expect_equal(rvcopls$Model$R2XC,
    #              c(0.18253923588794,0.176616467814135))
    ##TODO K, etc.
})
