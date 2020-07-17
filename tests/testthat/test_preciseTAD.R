context("test-preciseTAD")

test_that("Whether preciseTAD gives us the same output", {

    data(arrowhead_gm12878_5kb)

    bounds.GR <- extractBoundaries(domains.mat=arrowhead_gm12878_5kb,
                                   preprocess=FALSE,
                                   CHR=c("CHR21","CHR22"),
                                   resolution=5000)

    data(tfbsList)

    tfbsList_filt <- tfbsList[which(names(tfbsList) %in% c("Gm12878-Ctcf-Broad", "Gm12878-Rad21-Haib", "Gm12878-Smc3-Sydh", "Gm12878-Znf143-Sydh"))]

    tadData <- createTADdata(bounds.GR=bounds.GR,
                             resolution=5000,
                             genomicElements.GR=tfbsList_filt,
                             featureType="distance",
                             resampling="rus",
                             trainCHR="CHR21",
                             predictCHR="CHR22",
                             seed=123)

    tadModel <- TADrandomForest(trainData=tadData[[1]],
                                testData=tadData[[2]],
                                tuneParams=list(mtry=2,
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

    bounds.GR <- extractBoundaries(domains.mat=arrowhead_gm12878_5kb,
                                   preprocess=FALSE,
                                   CHR="CHR22",
                                   resolution=5000)

    pt <- preciseTAD(bounds.GR=bounds.GR,
                     genomicElements.GR=tfbsList_filt,
                     featureType="distance",
                     CHR="CHR22",
                     chromCoords=list(17000000,19000000),
                     tadModel=tadModel[[1]],
                     threshold=1.0,
                     flank=NULL,
                     verbose=TRUE,
                     seed=123,
                     parallel=TRUE,
                     cores=2,
                     splits=2,
                     DBSCAN=TRUE,
                     DBSCAN_params=list(5000,3),
                     method.Clust=NULL,
                     PTBR=TRUE,
                     CLARA=TRUE,
                     method.Dist="euclidean",
                     samples=100,
                     juicer=FALSE)

    expect_equal(length(pt[[1]]), 12)

    expect_equal(length(pt[[2]]), 7)

    expect_equal(IRanges::start(pt[[2]])[1], 17403327)
})
