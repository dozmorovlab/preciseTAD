context("test-juicer_func")

test_that("Whether juicer_func gives us the same output", {

    data(arrowhead_gm12878_5kb)

    bounds.GR <- extractBoundaries(domains.mat=arrowhead_gm12878_5kb,
                                   preprocess=FALSE,
                                   CHR=c("CHR21","CHR22"),
                                   resolution=5000)

    data(tfbsList)

    tfbsList_filt <- tfbsList[which(names(tfbsList) %in% c("Gm12878-Ctcf-Broad", "Gm12878-Rad21-Haib", "Gm12878-Smc3-Sydh", "Gm12878-Znf143-Sydh"))]

    set.seed(123)

    tadData <- createTADdata(bounds.GR=bounds.GR,
                             resolution=5000,
                             genomicElements.GR=tfbsList_filt,
                             featureType="distance",
                             resampling="rus",
                             trainCHR="CHR21",
                             predictCHR="CHR22")

    set.seed(123)

    tadModel <- TADrandomForest(trainData=tadData[[1]],
                                testData=tadData[[2]],
                                tuneParams=list(mtry=2,
                                                ntree=500,
                                                nodesize=1),
                                cvFolds=3,
                                cvMetric="Accuracy",
                                verbose=TRUE,
                                model=TRUE,
                                importances=TRUE,
                                impMeasure="MDA",
                                performances=TRUE)

    set.seed(123)

    pt <- preciseTAD(genomicElements.GR=tfbsList_filt,
                     featureType="distance",
                     CHR="CHR22",
                     chromCoords=list(17000000,19000000),
                     tadModel=tadModel[[1]],
                     threshold=1.0,
                     flank=NULL,
                     verbose=TRUE,
                     parallel=TRUE,
                     cores=2,
                     splits=2,
                     DBSCAN_params=list(5000,3),
                     samples=100)

    pt_juice <- juicer_func(pt[[2]])

    expect_equal(nrow(pt_juice), length(pt[[2]])-1)

    expect_equal(ncol(pt_juice), 6)
})
