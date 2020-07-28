context("test-TADrfe")

test_that("Whether TADrfe gives us the same output", {

    data(arrowhead_gm12878_5kb)

    bounds.GR <- extractBoundaries(domains.mat=arrowhead_gm12878_5kb,
                                   preprocess=FALSE,
                                   CHR="CHR22",
                                   resolution=5000)

    data(tfbsList)

    set.seed(123)

    tadData <- createTADdata(bounds.GR=bounds.GR,
                             resolution=5000,
                             genomicElements.GR=tfbsList,
                             featureType="oc",
                             resampling="rus",
                             trainCHR="CHR22",
                             predictCHR=NULL)

    set.seed(123)

    RFEres <- TADrfe(trainData=tadData[[1]],
                     tuneParams=list(ntree=100,
                                     nodesize=1),
                     cvFolds=5,
                     cvMetric="Accuracy",
                     verbose=TRUE)

    expect_equal(round(max(RFEres[[1]]$Accuracy),3), 0.642)

    expect_equal(round(mean(RFEres[[2]]$Overall),3), 5.058)

})
