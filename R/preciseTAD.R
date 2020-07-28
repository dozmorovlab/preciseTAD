#' Precise TAD boundary prediction at base-level resolution using density-based
#' spatial clustering and partitioning techniques
#'
#' @param bounds.GR \code{GRanges} object with chromosomal coordinates
#' of TAD boundaries used to identify positive cases in binary classification.
#' framework (can be obtained using \code{\link{extractBoundaries}})
#' @param genomicElements.GR \code{GRangesList} object containing GRanges from
#' each ChIP-seq BED file that was used to train a predictive model (can be
#' obtained using the \code{\link{bedToGRangesList}}).
#' @param featureType Controls how the feature space is constructed (one of
#' either "binary", "oc", "op", "signal, or "distance" (log2- transformed).
#' Default is "distance".
#' @param CHR Controls which chromosome to predict boundaries on at base-level
#' resolution.
#' @param chromCoords List containing the starting bp coordinate and ending bp
#' coordinate that defines the region of the linear genome to make predictions
#' on. If chromCoords is not specified, then predictions will be made on the
#' entire chromosome. Default is NULL.
#' @param tadModel Model object used to obtain predicted probabilities at base-level
#' resolution (examples include \code{gbm}, \code{glmnet},
#' \code{svm}, \code{glm}, etc). For a random forest model, can be obtained
#' using \code{preciseTAD::randomForest}).
#' @param threshold Bases with predicted probabilities that are greater
#' than or equal to this value are labeled as potential TAD boundaries. Values
#' in the range of .95-1.0 are suggested.
#' @param flank Controls how much to flank the final predicted TAD boundaries
#' (necessary for evaluating overlaps, etc.). Default is NULL, i.e., no flanking.
#' @param verbose Option to print progress.
#' @param parallel Option to parallelise the process for obtaining predicted
#' probabilities. Default is FALSE.
#' @param cores Number of cores to use in parallel. Default is NULL.
#' @param splits Number of splits of the test data to speed up the prediction.
#' Default is NULL.
#' @param DBSCAN Whether or not to use \code{\link{dbscan}} instead of
#' agglomerative hierarchical clustering (\code{\link{hclust}}) when clustering
#' pairwise genomic distance. Default is TRUE.
#' @param DBSCAN_params Parameters passed to \code{\link{dbscan}} in list form
#' containing 1) eps and 2) MinPts.
#' @param method.Clust The agglomeration method to be passed to
#' \code{\link{hclust}}. Default is NULL, indicating to use DBSCAN instead.
#' @param PTBR Option to include PTBRs along with predicted boundaries. Default
#' is TRUE.
#' @param CLARA Option to use CLARA (\code{\link{clara}}) instead of PAM
#' (\code{\link{pam}}). Default is TRUE.
#' @param method.Dist Distance metric passed to \code{\link{clara}} or
#' \code{\link{pam}} (if CLARA=FALSE). Default is "euclidean".
#' @param samples Number of subsamples if applying CLARA. Default is 100.
#' Ignored if CLARA=FALSE.
#' @param juicer Option to return predicted boundaries in a format that allows
#' for plotting in Juicebox from Aiden Lab.
#'
#' @return A list object containing at most 3 \code{GRanges} elements including:
#' 1) the genomic coordinates of the called TAD boundaries (CTBP) used to make
#' predictions, 2) the genomic coordinates of preciseTAD predicted
#' boundaries (PTBP) (if `juicer=TRUE`, this will be a data.frame that can be
#' saved and imported into Juicer as a BED file), and 2) the genomic coordinates
#' of preciseTAD predicted regions (PTBRs)if PTBR=TRUE (else NULL)
#' @export
#'
#' @importFrom pROC roc
#' @importFrom pROC coords
#' @importFrom stats predict
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @importFrom stats dist
#' @importFrom dbscan dbscan
#' @import pbapply parallel doSNOW foreach cluster IRanges
#' GenomicRanges
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
#' # Read in GRangesList of 26 TFBS
#' data(tfbsList)
#'
#' # Create the binned data matrix for CHR1 (training) and CHR22 (testing)
#' # using 5 kb binning, distance-type predictors from 26 different TFBS from
#' # the GM12878 cell line, and random under-sampling
#' tadData <- createTADdata(bounds.GR = bounds.GR,
#'                          resolution = 5000,
#'                          genomicElements.GR = tfbsList,
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
#' bounds.GR <- extractBoundaries(domains.mat = arrowhead_gm12878_5kb,
#'                                preprocess = FALSE,
#'                                CHR = "CHR22",
#'                                resolution = 5000)
#'
#' # Apply preciseTAD on a specific 2mb section of CHR22:17000000-19000000
#' pt <- preciseTAD(bounds.GR = bounds.GR,
#'                  genomicElements.GR = tfbsList,
#'                  featureType = "distance",
#'                  CHR = "CHR22",
#'                  chromCoords = list(17000000, 19000000),
#'                  tadModel = tadModel[[1]],
#'                  threshold = .95,
#'                  flank = NULL,
#'                  verbose = TRUE,
#'                  parallel = TRUE,
#'                  cores = 2,
#'                  splits = 2,
#'                  DBSCAN = TRUE,
#'                  DBSCAN_params = list(5000, 3),
#'                  method.Clust = NULL,
#'                  PTBR = TRUE,
#'                  CLARA = TRUE,
#'                  method.Dist = "euclidean",
#'                  samples = 100,
#'                  juicer = FALSE)
preciseTAD = function(bounds.GR, genomicElements.GR, featureType = "distance", CHR,
                      chromCoords = NULL, tadModel, threshold, flank = NULL, verbose = TRUE,
                      parallel = FALSE, cores = NULL, splits = NULL, DBSCAN = TRUE, DBSCAN_params,
                      method.Clust = NULL, PTBR = TRUE, CLARA = TRUE, method.Dist = "euclidean", samples = 100,
                      juicer = FALSE) {

    set.seed(123)

    #ESTABLISHING CHROMOSOME-SPECIFIC SEQINFO#

    #LOADING CHROMOSOME LENGTHS#

    hg19 <- preciseTAD:::hg19

    if(class(chromCoords)=="list") {
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

    #FUNCTIONS FOR DIFFERENT FEATURE TYPES#

    #BINARY OVERLAPS#
    binary_func <- function(binned_data_gr, annot_data_gr) {

        # Finding the total number of overlaps between genomic bins and the specific
        # genomic annotation
        count_binary <- countOverlaps(binned_data_gr, annot_data_gr)

        # Binarizing the overlap (1 or 0)
        count_binary <- ifelse(count_binary > 0, 1, 0)

        return(count_binary)
    }

    #COUNT OVERLAPS#
    count_func <- function(binned_data_gr, annot_data_gr) {

        # Finding the total number of overlaps between genomic bins and the specific
        # genomic annotation
        count_total <- countOverlaps(binned_data_gr, annot_data_gr)

        return(count_total)
    }

    #PERCENT OVERLAPS#
    percent_func <- function(binned_data_gr, annot_data_gr) {

        count_percent <- numeric(length(binned_data_gr))

        # Finding the total number of overlaps between genomic bins and the specific
        # genomic annotation
        c <- countOverlaps(binned_data_gr, annot_data_gr)

        # places where c=0 denotes no overlap places where c>0 denotes some type of
        # overlap, could be partial or within

        # for c=0 assign percent overlap as 0
        count_percent[which(c == 0)] <- 0

        # for c=1:
        count_percent[which(c == 1)] <- width(pintersect(findOverlapPairs(annot_data_gr,
                                                                          binned_data_gr[which(c == 1)])))/(width(binned_data_gr[1]))

        # for c>1: iterate through all bins with multiple overlaps with the annotation of
        # interest find how many and which annotations overlap within each iterate
        # calculate the width of each overlap sum the total widths and divide by bin
        # width bins with some type of overlap
        mo <- which(c > 1)
        count_percent[mo] <- unlist(lapply(lapply(lapply(mo, function(x) {
            width(pintersect(findOverlapPairs(annot_data_gr, binned_data_gr[x])))
        }), sum), function(x) {
            x/width(binned_data_gr[1])
        }))

        return(count_percent)
    }

    #SIGNAL#
    signal_func <- function(binned_data_gr, annot_data_gr){

        if(names(mcols(annot_data_gr))!="coverage"){print("metadata missing coverage column! use annots_to_granges_func function and specify signal parameter!"); return(0)}
        if(class(mcols(annot_data_gr)$coverage)!="numeric"){print("metadata coverage column is not numeric!")}

        count_signal <- numeric(length(binned_data_gr))

        #Finding the total number of overlaps between genomic bins and the specific genomic annotation
        c <- countOverlaps(binned_data_gr, annot_data_gr)

        #places where c=0 denotes no overlap
        #places where c>0 denotes some type of overlap, could be partial or within

        #for c=0 assign signal as 0
        count_signal[which(c==0)] <- 0

        #for c=1:
        count_signal[which(c==1)] <- mcols(annot_data_gr[queryHits(findOverlaps(annot_data_gr,binned_data_gr[which(c==1)]))])$coverage

        #for c>1:
        #iterate through all bins with multiple overlaps with the annotation of interest
        #find how many and which annotations overlap within each iterate
        #calculate the width of each overlap
        #sum the total widths and divide by bin width
        #bins with some type of overlap
        mo <- which(c>1)
        count_signal[mo] <- unlist(lapply(mo, function(x){
            sum(mcols(pintersect(findOverlapPairs(annot_data_gr,binned_data_gr[x])))$coverage)/c[x]
        }))

        return(count_signal)
    }

    #DISTANCE#
    distance_func <- function(binned_data_center_gr, annot_data_center_gr) {

        # distance from center of genomic bin to nearest genomic region of interest
        dist <- mcols(distanceToNearest(binned_data_center_gr, annot_data_center_gr))$distance

        return(dist)
    }

    #CREATING BP RESOLUTION TEST DATA#

    test_data <- matrix(nrow = length(seqDataTest),
                        ncol = length(genomicElements.GR),
                        dimnames = list(NULL,names(genomicElements.GR)))

    g <- split(GRanges(seqnames=tolower(CHR),IRanges(start=seqDataTest,end=seqDataTest)),
               ceiling(seq_along(GRanges(seqnames=tolower(CHR),IRanges(start=seqDataTest,end=seqDataTest)))/10000000))

    if(verbose==TRUE){print(paste0("Establishing bp resolution test data using a ", featureType, " type feature space"))}
    for(j in 1:length(genomicElements.GR)){
        p <-list()
        for(i in 1:length(g)){
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

    if(parallel == TRUE){
        parallel_predictions<-function(fit,testing,c,n){
            cl<-makeCluster(c)
            registerDoSNOW(cl)
            split_testing<-sort(rank(1:nrow(testing)) %% n)
            predictions<-foreach(i=unique(split_testing),
                                 .combine=c,
                                 .packages=c("caret")) %dopar% {
                                     as.numeric(predict(fit,newdata=testing[split_testing==i,],type="prob")[,"Yes"])
                                 }
            stopCluster(cl)
            predictions
        }
        predictions <- parallel_predictions(fit=tadModel,testing=test_data,c=cores,n=splits)
    }else{
        if(class(chromCoords) != "list"){
            array_split <- function(data, number_of_chunks){
                rowIdx <- seq_len(nrow(data))
                lapply(split(rowIdx, cut(rowIdx, pretty(rowIdx, number_of_chunks))), function(x) data[x, ])
            }
            split_ind <- as.integer(nrow(test_data)/10000000) + 1
            test_data <- array_split(test_data, split_ind)
            predictions <- list()
            for(i in 1:length(test_data)){
                predictions[[i]] <- predict(tadModel,newdata=test_data[[i]],type="prob")[,"Yes"]
            }
            predictions <- unlist(predictions)
        }else{
            predictions <- predict(tadModel,newdata=test_data,type="prob")[,"Yes"]
        }
    }

    if (threshold == "roc") {
        test_data_Y <- ifelse(seqDataTest %in% start(bounds.GR), 1, 0)
        t <- pROC::coords(pROC::roc(test_data_Y, predictions, quiet = TRUE),
                          "best", ret = "threshold", transpose = FALSE, best.method = "closest.topleft")
    } else {
        t <- threshold
    }

    if (verbose == TRUE) {
        print(paste0("preciseTAD identified a total of ", length(diff(seqDataTest[which(predictions >= t)])), " base pairs whose predictive probability was equal to or exceeded a threshold of ",
                     t))
    }

    retain <- c(1, cumsum(ifelse(diff(seqDataTest[which(predictions >= t)]) !=
                                     1, 1, 0)) + 1)

    if (verbose == TRUE) {
        print(paste0("preciseTAD identified ", length(unique(retain)), " PTBAs"))
        print("Establishing PTBRs")
    }

    mid <- tapply(seqDataTest[which(predictions >= t)], retain, function(x) {
        return(ceiling((x[1] + x[length(x)])/2))
    })

    mid <- as.vector(mid)

    if (DBSCAN == FALSE) {
        dist_mat <- dist(mid, method = "euclidean")
        if (verbose == TRUE) {
            print("Initializing hierarchical clustering")
        }
        hc1 <- hclust(dist_mat, method = method.Clust)
        k <- pbsapply(2:(length(mid) - 1), function(i) {
            mean(silhouette(cutree(hc1, i), dist = dist_mat)[, "sil_width"])
        })
        k = which.max(k) + 1
        if (verbose == TRUE) {
            print(paste0("preciseTAD identified ", k, " PTBRs"))
        }
    } else {
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
    }

    if (CLARA == TRUE) {
        if (verbose == TRUE) {
            print(paste0("Initializing CLARA with ", k, " clusters"))
            c <- clara(mid, k = k, samples = samples, metric = method.Dist, stand = FALSE,
                       trace = 2, medoids.x = TRUE, rngR = TRUE)
            medoids <- c$medoids
            clustering <- c$clustering
        } else {
            c <- clara(mid, k = k, samples = samples, metric = method.Dist, stand = FALSE,
                       trace = 0, medoids.x = TRUE, keep.data = FALSE, rngR = TRUE)
            medoids <- c$medoids
            clustering <- c$clustering
        }
    } else {
        if (verbose == TRUE) {
            print(paste0("Initializing PAM with ", k, " clusters"))
            c <- pam(x = mid, k = k, diss = FALSE, metric = method.Dist, medoids = NULL,
                     stand = FALSE, cluster.only = FALSE, do.swap = TRUE, keep.diss = FALSE,
                     trace.lev = 2)
            medoids <- c$medoids
            clustering <- c$clustering
        } else {
            c <- pam(x = mid, k = k, diss = FALSE, metric = method.Dist, medoids = NULL,
                     stand = FALSE, cluster.only = FALSE, do.swap = TRUE, keep.diss = FALSE,
                     trace.lev = 0)
            medoids <- c$medoids
            clustering <- c$clustering
        }
    }

    predBound_gr <- GRanges(seqnames = tolower(CHR), IRanges(start = seqDataTest[which(seqDataTest %in%
                                                                                           medoids)], end = seqDataTest[which(seqDataTest %in% medoids)]))

    trueBound_gr <- GRanges(seqnames = tolower(CHR), IRanges(start = seqDataTest[which(seqDataTest %in%
                                                                                           start(bounds.GR))], end = seqDataTest[which(seqDataTest %in% start(bounds.GR))]))

    if (!is.null(flank)) {
        predBound_gr <- flank(predBound_gr, width = flank, both = TRUE)
        trueBound_gr <- flank(trueBound_gr, width = flank, both = TRUE)
    }

    juicer_func <- function(grdat) {
        # n <- length(unique(as.character(seqnames(grdat))))

        chrs <- unique(as.character(seqnames(grdat)))

        mat_list <- list()

        for (i in 1:length(chrs)) {
            grdat_chr <- grdat[which(as.character(seqnames(grdat)) == chrs[i])]

            if (length(grdat_chr) > 1) {
                mymat <- matrix(nrow = (length(grdat_chr) - 1), ncol = 6)
                mymat[, 1] <- mymat[, 4] <- chrs[i]
                for (j in 1:(length(grdat_chr) - 1)) {
                    mymat[j, 2] <- mymat[j, 5] <- start(grdat_chr)[j]
                    mymat[j, 3] <- mymat[j, 6] <- start(grdat_chr)[j + 1]
                }
                mymat <- as.data.frame(mymat)

                mymat$V1 <- as.character(mymat$V1)
                mymat$V4 <- as.character(mymat$V4)
                mymat$V2 <- as.numeric(as.character(mymat$V2))
                mymat$V3 <- as.numeric(as.character(mymat$V3))
                mymat$V5 <- as.numeric(as.character(mymat$V5))
                mymat$V6 <- as.numeric(as.character(mymat$V6))

                mat_list[[i]] <- mymat

            }


        }

        mat_list <- do.call("rbind", mat_list)
        return(mat_list)

    }

    if (PTBR == TRUE) {
        grlist <- GRangesList()
        for (i in 1:length(unique(clustering))) {

            PTBAs = c(1:length(clustering))[which(clustering == i)][1]
            PTBAe = c(1:length(clustering))[which(clustering == i)][length(which(clustering ==
                                                                                     i))]

            PTBRs = mid[PTBAs]
            PTBRe = mid[PTBAe]

            grlist[[i]] <- GRanges(seqnames = tolower(CHR), IRanges(start = PTBRs,
                                                                    end = PTBRe))
        }
        gr <- unlist(grlist)

        if (juicer == TRUE) {
            bp_results <- list(CTBP = trueBound_gr, PTBP = juicer_func(predBound_gr),
                               PTBR = gr)
        } else {
            bp_results <- list(CTBP = trueBound_gr, PTBP = predBound_gr, PTBR = gr)
        }
    } else {
        if (juicer == TRUE) {
            bp_results <- list(CTBP = trueBound_gr, PTBP = juicer_func(predBound_gr),
                               PTBR = NULL)
        } else {
            bp_results <- list(CTBP = trueBound_gr, PTBP = predBound_gr, PTBR = NULL)
        }
    }

    return(bp_results)
}
