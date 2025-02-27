# test da functionality
copls.da <- ConsensusOPLS(
    data=lapply(demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")], scale),
    Y=demo_3_Omics$Y,
    maxPcomp=1,
    maxOcomp=3,
    nperm=100,
    modelType="da",
    mc.cores=1,
    nfold=14,
    verbose=T)

# test predict
pred.da <- predict(copls.da)
# predicted Y
expect_equal(pred.da$Y[6:9,1],
             c(HCT15=1.014562172,
               KM12=0.996538043,
               NCIADRRES=0.035633998,
               OVCAR3=-0.003387164),
             tolerance=1e-7
             )
# predicted class
expect_equal(pred.da$class$class[6:9],
             c("Colon", "Colon", "Ovarian", "Ovarian"))
expect_equal(pred.da$class$margin[3:5],
             c(0.8632528, 1.0196975, 1.0395004),
             tolerance=1e-7)
expect_equal(pred.da$class$softmax.Colon[3:5],
             c(0.9999968, 1, 1),
             tolerance=1e-7)

# test regression functionality
copls.reg <- ConsensusOPLS(
    data=lapply(demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")], scale),
    Y=matrix(c(rnorm(7, mean=1, sd=0.01), rnorm(7, mean=0, sd=0.01))),
    maxPcomp=1,
    maxOcomp=3,
    nperm=100,
    modelType="reg",
    mc.cores=1,
    kernelParams=list(type = "p", params = c(order = 2)),
    nfold=3, 
    verbose=T)
