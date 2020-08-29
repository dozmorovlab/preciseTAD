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
#' @param flank Controls how much to flank the final predicted TAD boundaries
#' (necessary for evaluating overlaps, etc.). Default is NULL, i.e., no
#' flanking.
#' @param verbose Option to print progress. Default is TRUE.
#' @param parallel Option to parallelise the process for obtaining predicted
#' probabilities. Must be number to indicate the number of cores to use in
#' parallel. Default is NULL.
#' @param DBSCAN_params Parameters passed to \code{\link{dbscan}} in list form
#' containing 1) eps and 2) MinPts. Required.
#'
#' @return A list object containing 2 \code{GRanges} elements including:
#' 1) the genomic coordinates of preciseTAD predicted regions (PTBRs), and 2)
#' the genomic coordinates of preciseTAD predicted boundaries (PTBP)
#' @export
#'
#' @importFrom pROC roc
#' @importFrom pROC coords
#' @importFrom stats predict
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @importFrom stats dist
#' @importFrom dbscan dbscan
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
#' tfbsList_filt <- tfbsList[which(names(tfbsList) %in% c("Gm12878-Ctcf-Broad", "Gm12878-Rad21-Haib", "Gm12878-Smc3-Sydh", "Gm12878-Znf143-Sydh"))]
#'
#' # Create the binned data matrix for CHR1 (training) and CHR22 (testing)
#' # using 5 kb binning, distance-type predictors from 26 different TFBS from
#' # the GM12878 cell line, and random under-sampling
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
#' tadModel <- TADrandomForest(trainData = tadData[[1]],
#'                             testData = tadData[[2]],
#'                             tuneParams = list(mtry = c(2, 5, 8, 10, 13, 16, 18, 21, 24, 26),
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
#' pt <- preciseTAD(genomicElements.GR = tfbsList_filt,
#'                  featureType = "distance",
#'                  CHR = "CHR22",
#'                  chromCoords = list(17000000, 19000000),
#'                  tadModel = tadModel[[1]],
#'                  threshold = 1,
#'                  flank = NULL,
#'                  verbose = TRUE,
#'                  parallel = NULL,
#'                  DBSCAN_params = list(5000, 3))
preciseTAD = function(genomicElements.GR, featureType = "distance", CHR,
                      chromCoords = NULL, tadModel, threshold = 1, flank = NULL,
                      verbose = TRUE, parallel = NULL, DBSCAN_params) {

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

    g <- split(GRanges(seqnames=tolower(CHR),IRanges(start=seqDataTest,end=seqDataTest)),
               ceiling(seq_along(GRanges(seqnames=tolower(CHR),IRanges(start=seqDataTest,end=seqDataTest)))/10000000))

    if(verbose==TRUE){print(paste0("Establishing bp resolution test data using a ", featureType, " type feature space"))}
    for(j in seq_len(length(genomicElements.GR))){
        p <-list()
        for(i in seq_len(length(g))){
            if(featureType=="distance"){
                p[[i]] <- log(distance_func(g[[i]],genomicElements.GR[[j]]) + 1, base = 2)
            }else if(featureType=="binary"){
                p[[i]] <- binary_func(g[[i]],genomicElements.GR[[j]])
            }else if(featureType=="oc"){
                p[[i]] <- count_func(g[[i]],genomicElements.GR[[j]])
            }else{
                p[[i]] <- percent_func(g[[i]],genomicElements.GR[[j]])
            }
        }
        p <- unlist(p)
        test_data[,j] <- p
    }

    rm("p", "g")

    #PREDICTING AT BP RESOLUTION #

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

    retain <- c(1, cumsum(ifelse(diff(seqDataTest[which(predictions >= threshold)]) !=
                                     1, 1, 0)) + 1)

    mid <- tapply(seqDataTest[which(predictions >= threshold)], retain, function(x) {
        return(ceiling((x[1] + x[length(x)])/2))
    })

    mid <- as.vector(mid)

    #if (DBSCAN == FALSE) {
    #    dist_mat <- dist(mid, method = "euclidean")
    #    if (verbose == TRUE) {
    #        print("Initializing hierarchical clustering")
    #    }
    #    hc1 <- hclust(dist_mat, method = method.Clust)
    #    k <- pbsapply(2:(length(mid) - 1), function(i) {
    #        mean(silhouette(cutree(hc1, i), dist = dist_mat)[, "sil_width"])
    #    })
    #    k = which.max(k) + 1
    #    if (verbose == TRUE) {
    #        print(paste0("preciseTAD identified ", k, " PTBRs"))
    #    }
    #} else {
        if (verbose == TRUE) {
            print("Initializing DBSCAN")
        }
        res <- dbscan::dbscan(as.matrix(mid), eps = DBSCAN_params[[1]], minPts = DBSCAN_params[[2]])
        if (0 %in% unique(res$cluster)) {
            k = length(unique(res$cluster)) - 1
        } else {
            k = length(unique(res$cluster))
        }
        if (verbose == TRUE) {
            print(paste0("preciseTAD identified ", k, " PTBRs"))
        }
    #}

    #if (CLARA == TRUE) {
        if (verbose == TRUE) {
            print(paste0("Initializing CLARA with ", k, " clusters"))
            c <- clara(mid, k = k, samples = 100, metric = "euclidean", stand = FALSE,
                       trace = 2, medoids.x = TRUE, rngR = TRUE)
            medoids <- c$medoids
            clustering <- c$clustering
        } else {
            c <- clara(mid, k = k, samples = 100, metric = "euclidean", stand = FALSE,
                       trace = 0, medoids.x = TRUE, keep.data = FALSE, rngR = TRUE)
            medoids <- c$medoids
            clustering <- c$clustering
        }
    #} else {
    #    if (verbose == TRUE) {
    #        print(paste0("Initializing PAM with ", k, " clusters"))
    #        c <- pam(x = mid, k = k, diss = FALSE, metric = method.Dist, medoids = NULL,
    #                 stand = FALSE, cluster.only = FALSE, do.swap = TRUE, keep.diss = FALSE,
    #                 trace.lev = 2)
    #        medoids <- c$medoids
    #        clustering <- c$clustering
    #    } else {
    #        c <- pam(x = mid, k = k, diss = FALSE, metric = method.Dist, medoids = NULL,
    #                 stand = FALSE, cluster.only = FALSE, do.swap = TRUE, keep.diss = FALSE,
    #                 trace.lev = 0)
    #        medoids <- c$medoids
    #        clustering <- c$clustering
    #    }
    #}

    predBound_gr <- GRanges(seqnames = tolower(CHR), IRanges(start = seqDataTest[which(seqDataTest %in%
                                                                                           medoids)], end = seqDataTest[which(seqDataTest %in% medoids)]))

    if (!is.null(flank)) {
        predBound_gr <- flank(predBound_gr, width = flank, both = TRUE)
    }

        grlist <- GRangesList()
        for (i in seq_len(length(unique(clustering)))) {

            PTBAs = c(seq_len(length(clustering)))[which(clustering == i)][1]
            PTBAe = c(seq_len(length(clustering)))[which(clustering == i)][length(which(clustering ==
                                                                                     i))]

            PTBRs = mid[PTBAs]
            PTBRe = mid[PTBAe]

            grlist[[i]] <- GRanges(seqnames = tolower(CHR), IRanges(start = PTBRs,
                                                                    end = PTBRe))
        }
        gr <- unlist(grlist)

        #if (juicer == TRUE) {
        #    bp_results <- list(PTBR = gr,
        #                       PTBP = juicer_func(predBound_gr))
        #} else {
            bp_results <- list(PTBR = gr,
                               PTBP = predBound_gr)
        #}

    return(bp_results)
}
