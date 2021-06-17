context("test-preciseTAD")

test_that("Whether preciseTAD gives us the same output", {

    data(arrowhead_gm12878_5kb)

    bounds.GR <- extractBoundaries(domains.mat=arrowhead_gm12878_5kb,
                                   filter=FALSE,
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

    # genomicElements.GR=tfbsList_filt; featureType="distance"; CHR="CHR22"; chromCoords=list(17000000,19000000); tadModel=tadModel[[1]]; threshold=0.99; verbose=TRUE; parallel=NULL; DBSCAN_params=list(eps=c(5000), MinPts=c(1000)); flank=5000; genome="hg19"; BaseProbs=FALSE
    pt <- preciseTAD(genomicElements.GR  = tfbsList_filt,
                     featureType         = "distance",
                     CHR                 = "CHR22",
                     chromCoords         = list(17000000,19000000),
                     tadModel            = tadModel[[1]],
                     threshold           = 0.99,
                     verbose             = TRUE,
                     parallel            = NULL,
                     DBSCAN_params       = list(eps = c(5000), 
                                                MinPts = c(1000)),
                     flank               = 5000,
                     genome              = "hg19",
                     BaseProbs           = FALSE)

    expect_equal(width(pt$PTBR)[1], 13752)

    expect_equal(pt$Summaries$PTBRWidth$median, 11475)

    expect_equal(pt$Summaries$PTBRCoverage$median, 0.3979696)

    expect_equal(length(pt$PTBP), 14)

    expect_equal(IRanges::start(pt$PTBP)[1], 17398701)
})
