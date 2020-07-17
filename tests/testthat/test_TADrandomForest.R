context("test-TADrandomForest")

test_that("Whether TADrandomForest gives us the same output", {

    data(arrowhead_gm12878_5kb)

    bounds.GR <- extractBoundaries(domains.mat=arrowhead_gm12878_5kb,
                                   preprocess=FALSE,
                                   CHR=c("CHR21","CHR22"),
                                   resolution=5000)

    data(tfbsList)

    tadData <- createTADdata(bounds.GR=bounds.GR,
                             resolution=5000,
                             genomicElements.GR=tfbsList,
                             featureType="distance",
                             resampling="rus",
                             trainCHR="CHR21",
                             predictCHR="CHR22",
                             seed=123)

    tadModel <- TADrandomForest(trainData=tadData[[1]],
                                testData=tadData[[2]],
                                tuneParams=list(mtry=ceiling(seq.int(2,ncol(tadData[[1]])-1,length.out = 10)),
                                                ntree=500,
                                                nodesize=1),
                                cvFolds=3,
                                cvMetric="Accuracy",
                                verbose=TRUE,
                                seed=123,
                                model=TRUE,
                                importances=TRUE,
                                impMeasure="MDA",
                                performances=TRUE)

    expect_equal(tadModel[[1]]$bestTune$mtry, 8)

    expect_equal(round(tadModel[[2]]$Importance[which(tadModel[[2]]$Feature=="`Gm12878-Ctcf-Broad`")],3), 20.526)

    expect_equal(round(tadModel[[3]]$Performance[which(tadModel[[3]]$Metric=="BalancedAccuracy")],3), 0.762)
})
