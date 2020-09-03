#' Precise TAD boundary prediction at base-level resolution using density-based
#' spatial clustering and partitioning techniques
#'
#' @param genomicElements.GR \code{GRangesList} object containing GRanges from
#' each ChIP-seq BED file that was used to train a predictive model (can be
#' obtained using the \code{\link{bedToGRangesList}}). Required.
#' @param featureType Controls how the feature space is constructed (one of
#' either "binary", "oc", "op", "signal, or "distance" (log2- transformed).
#' Default is "distance".
#' @param CHR Controls which chromosome to predict boundaries on at base-level
#' resolution. Required.
#' @param chromCoords List containing the starting bp coordinate and ending bp
#' coordinate that defines the region of the linear genome to make predictions
#' on. If chromCoords is not specified, then predictions will be made on the
#' entire chromosome. Default is NULL.
#' @param tadModel Model object used to obtain predicted probabilities at
#' base-level resolution (examples include \code{gbm}, \code{glmnet},
#' \code{svm}, \code{glm}, etc). For a random forest model, can be obtained
#' using \code{preciseTAD::randomForest}). Required.
#' @param threshold Bases with predicted probabilities that are greater
#' than or equal to this value are labeled as potential TAD boundaries. Values
#' in the range of .95-1.0 are suggested. Default is 1.
#' @param verbose Option to print progress. Default is TRUE.
#' @param parallel Option to parallelise the process for obtaining predicted
#' probabilities. Must be number to indicate the number of cores to use in
#' parallel. Default is NULL.
#' @param DBSCAN_params Parameters passed to \code{\link{dbscan}} in list form
#' containing 1) eps and 2) MinPts. Required.
#' @param flank Controls how much to flank the predicted TAD boundaries for
#' calculating normalized enrichment. Normalized enrichment is calculated as
#' the total number of peak regions that overlap with flanked predicted
#' boundaries divided by the number of predicted boundaries. Recommended value
#' is resolution. Required.
#'
#' @return A list containing 3 elements including:
#' 1) the genomic coordinates spanning each preciseTAD predicted region (PTBR),
#' 2) the genomic coordinates of preciseTAD predicted boundaries points (PTBP).
#' 3) a named list including summary statistics of the following:
#' PTBRWidth - PTBR width, PTBRCoverage - the proportion of bases within a PTBR
#' with probabilities that equal to or exceed the threshold (t=1 by default),
#' DistanceBetweenPTBR - the genomic distance between the end of the previous
#' PTBR and the start of the subsequent PTBR, NumSubRegions - the number of
#' the subregions in each PTBR cluster, SubRegionWidth - the width of
#' the subregion forming each PTBR, DistBetweenSubRegions -
#' the genomic distance between the end of the previous PTBR-specific subregion
#' and the start of the subsequent PTBR-specific subregion, and the normalized
#' enrichment of the genomic annotations used in the model around flanked PTBPs.
#' @export
#'
#' @importFrom pROC roc
#' @importFrom pROC coords
#' @importFrom stats predict
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @importFrom stats dist
#' @importFrom dbscan dbscan
#' @importFrom S4Vectors subjectHits
#' @import pbapply parallel doSNOW foreach cluster IRanges GenomicRanges
#'
#' @examples
#' # Read in ARROWHEAD-called TADs at 5kb
#' data(arrowhead_gm12878_5kb)
#'
#' # Extract unique boundaries
#' bounds.GR <- extractBoundaries(domains.mat = arrowhead_gm12878_5kb,
#'                                preprocess = FALSE,
#'                                CHR = c("CHR21", "CHR22"),
#'                                resolution = 5000)
#'
#' # Read in GRangesList of 26 TFBS and filter to include only CTCF, RAD21,
#' #SMC3, and ZNF143
#' data(tfbsList)
#'
#' tfbsList_filt <- tfbsList[which(names(tfbsList) %in%
#'                                                  c("Gm12878-Ctcf-Broad",
#'                                                    "Gm12878-Rad21-Haib",
#'                                                    "Gm12878-Smc3-Sydh",
#'                                                    "Gm12878-Znf143-Sydh"))]
#'
#' # Create the binned data matrix for CHR1 (training) and CHR22 (testing)
#' # using 5 kb binning, distance-type predictors from 4 TFBS from
#' # the GM12878 cell line, and random under-sampling
#' set.seed(123)
#' tadData <- createTADdata(bounds.GR = bounds.GR,
#'                          resolution = 5000,
#'                          genomicElements.GR = tfbsList_filt,
#'                          featureType = "distance",
#'                          resampling = "rus",
#'                          trainCHR = "CHR21",
#'                          predictCHR = "CHR22")
#'
#' # Perform random forest using TADrandomForest by tuning mtry over 10 values
#' # using 3-fold CV
#' set.seed(123)
#' tadModel <- TADrandomForest(trainData = tadData[[1]],
#'                             testData = tadData[[2]],
#'                             tuneParams = list(mtry = 2,
#'                                             ntree = 500,
#'                                             nodesize = 1),
#'                             cvFolds = 3,
#'                             cvMetric = "Accuracy",
#'                             verbose = TRUE,
#'                             model = TRUE,
#'                             importances = TRUE,
#'                             impMeasure = "MDA",
#'                             performances = TRUE)
#'
#' # Apply preciseTAD on a specific 2mb section of CHR22:17000000-19000000
#' set.seed(123)
#' pt <- preciseTAD(genomicElements.GR = tfbsList_filt,
#'                  featureType = "distance",
#'                  CHR = "CHR22",
#'                  chromCoords = list(17000000, 19000000),
#'                  tadModel = tadModel[[1]],
#'                  threshold = 1.0,
#'                  verbose = TRUE,
#'                  parallel = NULL,
#'                  DBSCAN_params = list(10000, 3),
#'                  flank = 5000)
preciseTAD = function(genomicElements.GR, featureType = "distance", CHR,
                      chromCoords = NULL, tadModel, threshold = 1,
                      verbose = TRUE, parallel = NULL, DBSCAN_params,
                      flank) {

    #ESTABLISHING CHROMOSOME-SPECIFIC SEQINFO#

    #LOADING CHROMOSOME LENGTHS#

    hg19 <- preciseTAD:::hg19

    if(is.list(chromCoords)) {
        seqDataTest <- c(chromCoords[[1]]:chromCoords[[2]])
        if("TRUE" %in% table(seqDataTest %in% c(hg19$centromerStart[hg19$chrom == CHR]:hg19$centromerEnd[hg19$chrom == CHR]))) {
            centromereTestStart <- hg19$centromerStart[hg19$chrom==CHR]
            centromereTestEnd <- hg19$centromerEnd[hg19$chrom==CHR]
            seqDataTest <- seqDataTest[-which(seqDataTest %in% c(centromereTestStart:centromereTestEnd))]
        }
    }else {
        seqLengthTest <- hg19$length[hg19$chrom == CHR]
        seqDataTest <- seq_along(0:seqLengthTest)
        centromereTestStart <- hg19$centromerStart[hg19$chrom == CHR]
        centromereTestEnd <- hg19$centromerEnd[hg19$chrom == CHR]
        seqDataTest <- seqDataTest[-which(seqDataTest %in% c(centromereTestStart:centromereTestEnd))]
    }

    #CREATING BP RESOLUTION TEST DATA#

    test_data <- matrix(nrow = length(seqDataTest),
                        ncol = length(genomicElements.GR),
                        dimnames = list(NULL,names(genomicElements.GR)))

    #g <- split(GRanges(seqnames=tolower(CHR),IRanges(start=seqDataTest,end=seqDataTest)),
    #           ceiling(seq_along(GRanges(seqnames=tolower(CHR),IRanges(start=seqDataTest,end=seqDataTest)))/10000000))

    if(verbose==TRUE){print(paste0("Establishing bp resolution test data using a ", featureType, " type feature space"))}
    for(j in seq_len(length(genomicElements.GR))){
        p <-list()
        for(i in seq_len(length(split(GRanges(seqnames=tolower(CHR),
                                        IRanges(start=seqDataTest,
                                                end=seqDataTest)),
                                ceiling(seq_along(GRanges(seqnames=tolower(CHR),
                                                          IRanges(start=seqDataTest,
                                                                  end=seqDataTest)))/10000000))))){
            if(featureType=="distance"){
                p[[i]] <- log(distance_func(split(GRanges(seqnames=tolower(CHR),
                                                          IRanges(start=seqDataTest,
                                                                  end=seqDataTest)),
                                                  ceiling(seq_along(GRanges(seqnames=tolower(CHR),
                                                                            IRanges(start=seqDataTest,
                                                                                    end=seqDataTest)))/10000000))[[i]],
                                            genomicElements.GR[[j]]) + 1,
                              base = 2)
            }else if(featureType=="binary"){
                p[[i]] <- binary_func(split(GRanges(seqnames=tolower(CHR),
                                                    IRanges(start=seqDataTest,
                                                            end=seqDataTest)),
                                            ceiling(seq_along(GRanges(seqnames=tolower(CHR),
                                                                      IRanges(start=seqDataTest,
                                                                              end=seqDataTest)))/10000000))[[i]],
                                      genomicElements.GR[[j]])
            }else if(featureType=="oc"){
                p[[i]] <- count_func(split(GRanges(seqnames=tolower(CHR),
                                                   IRanges(start=seqDataTest,
                                                           end=seqDataTest)),
                                           ceiling(seq_along(GRanges(seqnames=tolower(CHR),
                                                                     IRanges(start=seqDataTest,
                                                                             end=seqDataTest)))/10000000))[[i]],
                                     genomicElements.GR[[j]])
            }else{
                p[[i]] <- percent_func(split(GRanges(seqnames=tolower(CHR),
                                                     IRanges(start=seqDataTest,
                                                             end=seqDataTest)),
                                             ceiling(seq_along(GRanges(seqnames=tolower(CHR),
                                                                       IRanges(start=seqDataTest,
                                                                               end=seqDataTest)))/10000000))[[i]],
                                       genomicElements.GR[[j]])
            }
        }
        p <- unlist(p)
        test_data[,j] <- p
    }

    #rm("p", "g")

    #PREDICTING AT BP RESOLUTION #

    if (verbose == TRUE) {
        print("Establishing probability vector")
    }

    if(!is.null(parallel)){
        parallel_predictions<-function(fit, testing, c, n){
            cl<-makeCluster(c)
            registerDoSNOW(cl)
            split_testing<-sort(rank(seq_len(nrow(testing))) %% n)
            predictions<-foreach(i=unique(split_testing),
                                 .combine=c,
                                 .packages=c("caret")) %dopar% {
                                     as.numeric(predict(fit, newdata=testing[split_testing==i, ], type="prob")[, "Yes"])
                                 }
            stopCluster(cl)
            predictions
        }
        predictions <- parallel_predictions(fit=tadModel, testing=test_data, c=parallel, n=parallel)
    }else{
        if(!is.list(chromCoords)){
            array_split <- function(data, number_of_chunks){
                rowIdx <- seq_len(nrow(data))
                lapply(split(rowIdx, cut(rowIdx, pretty(rowIdx, number_of_chunks))), function(x) data[x, ])
            }
            split_ind <- as.integer(nrow(test_data)/10000000) + 1
            test_data <- array_split(test_data, split_ind)
            predictions <- list()
            for(i in seq_len(length(test_data))){
                predictions[[i]] <- predict(tadModel,newdata=test_data[[i]],type="prob")[,"Yes"]
            }
            predictions <- unlist(predictions)
        }else{
            predictions <- predict(tadModel,newdata=test_data,type="prob")[,"Yes"]
        }
    }

    if (verbose == TRUE) {
        print(paste0("preciseTAD identified a total of ", length(diff(seqDataTest[which(predictions >= threshold)])), " base pairs whose predictive probability was equal to or exceeded a threshold of ",
                     threshold))
    }

    if (verbose == TRUE) {
        print("Initializing DBSCAN")
    }
    res <- dbscan::dbscan(as.matrix(seqDataTest[which(predictions>=threshold)]), eps = DBSCAN_params[[1]], minPts = DBSCAN_params[[2]])
    if (0 %in% unique(res$cluster)) {
        k = length(unique(res$cluster)) - 1
    } else {
        k = length(unique(res$cluster))
    }

    if (verbose == TRUE) {
        print(paste0("preciseTAD identified ", k, " PTBRs"))
        print(paste0("Establishing PTBPs"))
    }

    medoids <- numeric()
    grlist <- GRangesList()
    numsubreg <- numeric()
    distbtwsubreg <- numeric()
    coverage <- numeric()
    widthsubreg <- numeric()
    for(i in seq_len(k)){
        if(verbose == TRUE){print(paste0("Cluster ", i, " out of ", k))}

        bp <- seqDataTest[which(predictions>=threshold)][which(res$cluster==i)]

        grlist[[i]] <- GRanges(seqnames = tolower(CHR),
                               IRanges(start=min(bp),
                                       end=max(bp)))

        if(0 %in% unique(res$cluster)){
            coverage[i] <- as.vector(table(res$cluster))[-1][i]/(max(bp)-min(bp)+1)
        }else{
            coverage[i] <- as.vector(table(res$cluster))[i]/(max(bp)-min(bp)+1)
        }

        if(identical(diff(bp[-which(diff(bp)==1)]), integer(0))){
            numsubreg[i] <- 1
            distbtwsubreg <- c(distbtwsubreg, 0)
        }else{
            numsubreg[i] <- length(diff(bp)[-which(diff(bp)==1)])+1
            distbtwsubreg <- c(distbtwsubreg, diff(bp)[-which(diff(bp)==1)])
        }

        retain <- c(1,cumsum(ifelse(diff(bp) != 1, 1, 0)) + 1)

        widthsubreg <- c(widthsubreg, as.vector(table(retain)))

        c <- pam(x = as.matrix(bp),
                 k = 1,
                 diss = FALSE,
                 metric = "euclidean",
                 medoids = NULL,
                 stand = FALSE,
                 cluster.only = FALSE,
                 do.swap = TRUE,
                 keep.diss = FALSE,
                 trace.lev = 0)
        medoids[i] <- c$medoids[1]

    }

    grlist <- unlist(grlist)

    predBound_gr <- GRanges(seqnames=tolower(CHR),
                            IRanges(start=medoids,
                                    end=medoids))

    if(0 %in% unique(res$cluster)){res$cluster <- res$cluster[-which(res$cluster==0)]}

    bp_results <- list(PTBR=grlist,
                       PTBP=predBound_gr,
                       Summaries=list(PTBRWidth = data.frame(min=min(width(grlist)),
                                                             max=max(width(grlist)),
                                                             median=median(width(grlist)),
                                                             iqr=IQR(width(grlist)),
                                                             mean=mean(width(grlist)),
                                                             sd=sd(width(grlist))),
                                      #TotalCoordperPTBR = data.frame(min=min(as.vector(table(res$cluster))),
                                      #                               max=max(as.vector(table(res$cluster))),
                                      #                               median=median(as.vector(table(res$cluster))),
                                      #                               iqr=IQR(as.vector(table(res$cluster))),
                                      #                               mean=mean(as.vector(table(res$cluster))),
                                      #                               sd=sd(as.vector(table(res$cluster)))),
                                      PTBRCoverage = data.frame(min=min(coverage),
                                                                max=max(coverage),
                                                                median=median(coverage),
                                                                iqr=IQR(coverage),
                                                                mean=mean(coverage),
                                                                sd=sd(coverage)),
                                      DistanceBetweenPTBR = data.frame(min=min(start(grlist)[-1]-end(grlist)[-length(end(grlist))]),
                                                                       max=max(start(grlist)[-1]-end(grlist)[-length(end(grlist))]),
                                                                       median=median(start(grlist)[-1]-end(grlist)[-length(end(grlist))]),
                                                                       iqr=IQR(start(grlist)[-1]-end(grlist)[-length(end(grlist))]),
                                                                       mean=mean(start(grlist)[-1]-end(grlist)[-length(end(grlist))]),
                                                                       sd=sd(start(grlist)[-1]-end(grlist)[-length(end(grlist))])),
                                      NumSubRegions = data.frame(min=min(numsubreg),
                                                                 max=max(numsubreg),
                                                                 median=median(numsubreg),
                                                                 iqr=IQR(numsubreg),
                                                                 mean=mean(numsubreg),
                                                                 sd=sd(numsubreg)),
                                      SubRegionWidth = data.frame(min=min(widthsubreg),
                                                                  max=max(widthsubreg),
                                                                  median=median(widthsubreg),
                                                                  iqr=IQR(widthsubreg),
                                                                  mean=mean(widthsubreg),
                                                                  sd=sd(widthsubreg)),
                                      DistBetweenSubRegions = data.frame(min=min(distbtwsubreg),
                                                                         max=max(distbtwsubreg),
                                                                         median=median(distbtwsubreg),
                                                                         iqr=IQR(distbtwsubreg),
                                                                         mean=mean(distbtwsubreg),
                                                                         sd=sd(distbtwsubreg)),
                                      NormilizedEnrichment = unlist(
                                          lapply(genomicElements.GR,
                                                 function(x){length(unique(subjectHits(findOverlaps(flank(predBound_gr,
                                                                                                          width = 5000,
                                                                                                          both = TRUE),
                                                                                                    x))))/length(predBound_gr)
                                                 }
                                          ))))

    return(bp_results)
}
