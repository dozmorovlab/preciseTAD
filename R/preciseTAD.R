#' Precise TAD boundary prediction at base pair resolution using density-based
#' spatial clustering and partitioning techniques
#'
#' @param bounds.GR \code{GRanges} object with chromosomal coordinates
#' of TAD boundaries used to identify positive cases in a binary classification.
#' framework (can be obtained using \code{\link{extractBoundaries}})
#' @param genomicElements.GR \code{GRangesList} object containing GRanges from
#' each ChIP-seq BED file that was used to train a predictive model (can be
#' obtained using the \code{\link{bedToGRangesList}}).
#' @param featureType Controls how the feature space is constructed (one of
#' either "binary", "oc", "op", or "distance" (log2- transformed). Default is
#' "distance".
#' @param CHR Controls which chromosome to predict boundaries on at base pair
#' resolution.
#' @param chromCoords List containing the starting bp coordinate and ending bp
#' coordinate that defines the region of the linear genome to make predictions
#' on. If chromCoords is not specified then predictions will be made on the
#' entire chromosome. Default is NULL.
#' @param tadModel Model object used to obtain predicted probabilities at base
#' pair resolution (examples include \code{gbm}, \code{glmnet},
#' \code{svm}, \code{glm}, etc). For a random forest model, can be obtained
#' using \code{preciseTAD::randomForest}).
#' @param threshold Base pairs with predicted probabilities that are greater
#' than or equal to this value are labeled as potential TAD boundaries. Values
#' in the range of .95-1.0 are suggested.
#' @param flank Controls how much to flank the final predicted TAD boundaries
#' (necessary for evaluating overlaps, etc.). Default is NULL, i.e. no flanking.
#' @param verbose Option to print progress.
#' @param seed Numeric for reproducibility.
#' @param parallel Option to parallelise the process for obtaining predicted
#' probabilities. Default is FALSE.
#' @param cores Number of cores to use in parallel. Default is NULL.
#' @param splits Number of splits of the test data to speed up prediction.
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
#' for plotting in juicebox from Aiden Lab.
#'
#' @return A list object containing at most 3 \code{GRanges} elements including:
#' 1) the genomic coordinates of preciseTAD predicted regions (PTBRs)if
#' PTBR=TRUE (else NA), 2) the genomic coordinates of preciseTAD predicted
#' boundaries (PTBP) (if `juicer=TRUE`, this will be a data.frame that can be
#' saved and imported into juicer as a txt file), and 3) the genomic coordinates
#' of the called TAD boundaries (CTBP) used to make predictions.
#' @export
#'
#' @importFrom pROC roc
#' @importFrom pROC coords
#' @importFrom stats predict
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @importFrom stats dist
#' @importFrom dbscan dbscan
#' @import pbapply parallel doSNOW foreach cluster bigmemory IRanges
#' GenomicRanges
#'
#' @examples
#' \dontrun{
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
#'                          predictCHR = "CHR22",
#'                          seed = 123)
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
#'                             seed = 123,
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
#'                  seed = 123,
#'                  parallel = TRUE,
#'                  cores = 4,
#'                  splits = 4,
#'                  DBSCAN = TRUE,
#'                  DBSCAN_params = list(5000, 3),
#'                  method.Clust = NULL,
#'                  PTBR = TRUE,
#'                  CLARA = TRUE,
#'                  method.Dist = "euclidean",
#'                  samples = 100,
#'                  juicer = FALSE)
#' }
preciseTAD = function(bounds.GR,
                      genomicElements.GR,
                      featureType="distance",
                      CHR,
                      chromCoords=NULL,
                      tadModel,
                      threshold,
                      flank=NULL,
                      verbose=TRUE,
                      seed=123,
                      parallel=FALSE,
                      cores=NULL,
                      splits=NULL,
                      DBSCAN=TRUE,
                      DBSCAN_params,
                      method.Clust=NULL,
                      PTBR=TRUE,
                      CLARA=TRUE,
                      method.Dist="euclidean",
                      samples=100,
                      juicer=FALSE){

    ##########################
    #CHECKING FUNCTION INPUTS#
    ##########################

    if(class(bounds.GR)!="GRanges"){print("is not a GRanges object!"); return(0)}
    if(class(genomicElements.GR)!="CompressedGRangesList"){print("genomicElements.GR is not a CompressedGRangesList object!"); return(0)}
    for(i in 1:length(genomicElements.GR)){
        if(class(genomicElements.GR[[i]])!="GRanges"){print(paste0(i, "-th object of genomicElements.GR is not a GenomicRanges object!")); return(0)}
    }
    if(class(featureType)!="character"){print("featureType is not a character object!"); return(0)}
    if(!(featureType %in% c("distance","binary","oc","op"))){print("featureType must be one of either 'distance', 'binary', 'oc', or 'op'!"); return(0)}
    if(class(CHR)!="character"){print("CHR is not a character object!"); return(0)}
    for(i in 1:length(CHR)){
        if(grepl("CHR",CHR[i])!=TRUE){print(paste0(i, "-th chromosome for training is not in 'CHR' format!")); return(0)}
    }

    if(FALSE %in% (tolower(CHR) %in% unique(as.character(seqnames(bounds.GR))))){print("CHR parameter is not in the bound.GR parameter!"); return(0)}

    if(class(chromCoords)=="list" & length(chromCoords)!=2){print("chromCoords is a list that is not length 2!"); return(0)}
    if(class(chromCoords)=="list"){
        if((chromCoords[[1]]>chromCoords[[2]])){print("starting base pair coordinate is greater than ending base pair coordinate!"); return(0)}
    }
    if(class(threshold)=="character" & threshold!="roc"){print("threshold is not a number and therefore must be 'roc'!"); return(0)}
    if(class(threshold)=="numeric" & (threshold>1 | threshold<0)){print("if threshold is a number, it must be between 0 and 1!"); return(0)}
    if(class(flank)!="NULL" & class(flank)!="numeric"){print("if flank is not NULL then flank must be a number > 0 representing number of base pairs!"); return(0)}
    if(class(verbose)!="logical"){print("verbose is not a logical object!"); return(0)}
    if(class(seed)!="numeric"){print("seed is not a numeric object!"); return(0)}
    if(class(parallel)!="logical"){print("parallel is not a logical object!"); return(0)}
    if(parallel==TRUE){
        if(class(cores)!="numeric"){print("cores is not a logical object!"); return(0)}
        if(class(splits)!="numeric"){print("splits is not a logical object!"); return(0)}
    }
    if(class(method.Dist)!="character"){print("distance metric is not a character object!"); return(0)}
    if(!(method.Dist %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))){print("method.Dist must be one of either 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'!"); return(0)}
    if(class(DBSCAN)!="logical"){print("DBSCAN is not a logical object!"); return(0)}
    if(DBSCAN==TRUE){
        if(class(DBSCAN_params)!="list"){print("DBSCAN_params is not a list object!"); return(0)}
        if(length(DBSCAN_params)!=2){print("DBSCAN_params is missing either eps or MinPts parameters!"); return(0)}
    }else{
        if(method.Clust!="character"){print("cluster metric is not a character object!"); return(0)}
        if(!(method.Clust %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))){print("method.Clust must be one of either 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'!"); return(0)}
    }
    if(class(CLARA)!="logical"){print("CLARA is not a logical object!"); return(0)}
    if(CLARA==TRUE & class(samples)!="numeric"){print("samples is not a numeric object!"); return(0)}
    if(class(juicer)!="logical"){print("juicer is not a logical object!"); return(0)}

    ##########################################
    #ESTABLISHING CHROMOSOME-SPECIFIC SEQINFO#
    ##########################################

    ############################
    #Loading chromosome lengths#
    ############################

    hg19 <- preciseTAD:::hg19

    if(class(chromCoords)=="list"){
        seqDataTest <- c(chromCoords[[1]]:chromCoords[[2]])

        if("TRUE" %in% table(seqDataTest %in% c(hg19$centromerStart[hg19$chrom==CHR]:hg19$centromerEnd[hg19$chrom==CHR]))){
            centromereTestStart <- hg19$centromerStart[hg19$chrom==CHR]
            centromereTestEnd <- hg19$centromerEnd[hg19$chrom==CHR]
            seqDataTest <- seqDataTest[-which(seqDataTest %in% c(centromereTestStart:centromereTestEnd))]
        }
    }else{
        seqLengthTest <- hg19$length[hg19$chrom==CHR]

        seqDataTest <- as.big.matrix(c(0:seqLengthTest))
        centromereTestStart <- hg19$centromerStart[hg19$chrom==CHR]
        centromereTestEnd <- hg19$centromerEnd[hg19$chrom==CHR]
        seqDataTest <- seqDataTest[,1][-which(seqDataTest[,1] %in% c(centromereTestStart:centromereTestEnd))]
    }

    #######################################
    #Functions for different feature types#
    #######################################

    ####calculating binary overlaps
    binary_func <- function(binned_data_gr,
                            annot_data_gr){

        #Finding the total number of overlaps between genomic bins and the specific genomic annotation
        count_binary <- countOverlaps(binned_data_gr, annot_data_gr)

        #Binarizing the overlap (1 or 0)
        count_binary <- ifelse(count_binary>0,1,0)

        return(count_binary)
    }

    ####calculating count overlaps
    count_func <- function(binned_data_gr,
                           annot_data_gr){

        #Finding the total number of overlaps between genomic bins and the specific genomic annotation
        count_total <- countOverlaps(binned_data_gr, annot_data_gr)

        return(count_total)
    }

    ####calculating percent overlaps
    percent_func <- function(binned_data_gr,
                             annot_data_gr){

        count_percent <- numeric(length(binned_data_gr))

        #Finding the total number of overlaps between genomic bins and the specific genomic annotation
        c <- countOverlaps(binned_data_gr, annot_data_gr)

        #places where c=0 denotes no overlap
        #places where c>0 denotes some type of overlap, could be partial or within

        #for c=0 assign percent overlap as 0
        count_percent[which(c==0)] <- 0

        #for c=1:
        count_percent[which(c==1)] <- width(pintersect(findOverlapPairs(annot_data_gr,binned_data_gr[which(c==1)])))/(width(binned_data_gr[1]))

        #for c>1:
        #iterate through all bins with multiple overlaps with the annotation of interest
        #find how many and which annotations overlap within each iterate
        #calculate the width of each overlap
        #sum the total widths and divide by bin width
        #bins with some type of overlap
        mo <- which(c>1)
        count_percent[mo] <- unlist(
            lapply(
                lapply(
                    lapply(
                        mo, function(x){width(pintersect(findOverlapPairs(annot_data_gr,binned_data_gr[x])))}),
                    sum),
                function(x){x/width(binned_data_gr[1])}))

        return(count_percent)
    }

    ####calculating distance
    distance_func <- function(binned_data_center_gr,
                              annot_data_center_gr){

        #distance from center of genomic bin to nearest genomic region of interest
        dist <- mcols(distanceToNearest(binned_data_center_gr, annot_data_center_gr))$distance

        return(dist)
    }

    ##################################
    #CREATING BP RESOLUTION TEST DATA#
    ##################################

    if(verbose==TRUE){
        print(paste0("Establishing bp resolution test data using a ", featureType, " type feature space"))
        if(featureType=="distance"){
            p <- pblapply(genomicElements.GR, function(x){as.big.matrix(log(distance_func(GRanges(seqnames=tolower(CHR),
                                                                                                  IRanges(start=seqDataTest,
                                                                                                          end=seqDataTest)),
                                                                                          x) + 1, base = 2))})
        }else if(featureType=="binary"){
            p <- pblapply(genomicElements.GR, function(x){as.big.matrix(binary_func(GRanges(seqnames=tolower(CHR),
                                                                                            IRanges(start=seqDataTest,
                                                                                                    end=seqDataTest)),
                                                                                    x))})
        }else if(featureType=="oc"){
            p <- pblapply(genomicElements.GR, function(x){as.big.matrix(count_func(GRanges(seqnames=tolower(CHR),
                                                                                           IRanges(start=seqDataTest,
                                                                                                   end=seqDataTest)),
                                                                                   x))})
        }else if(featureType=="op"){
            p <- pblapply(genomicElements.GR, function(x){as.big.matrix(percent_func(GRanges(seqnames=tolower(CHR),
                                                                                             IRanges(start=seqDataTest,
                                                                                                     end=seqDataTest)),
                                                                                     x))})
        }
    }else{
        if(featureType=="distance"){
            p <- lapply(genomicElements.GR, function(x){as.big.matrix(log(distance_func(GRanges(seqnames=tolower(CHR),
                                                                                                IRanges(start=seqDataTest,
                                                                                                        end=seqDataTest)),
                                                                                        x) + 1, base = 2))})
        }else if(featureType=="binary"){
            p <- lapply(genomicElements.GR, function(x){as.big.matrix(binary_func(GRanges(seqnames=tolower(CHR),
                                                                                          IRanges(start=seqDataTest,
                                                                                                  end=seqDataTest)),
                                                                                  x))})
        }else if(featureType=="oc"){
            p <- lapply(genomicElements.GR, function(x){as.big.matrix(count_func(GRanges(seqnames=tolower(CHR),
                                                                                         IRanges(start=seqDataTest,
                                                                                                 end=seqDataTest)),
                                                                                 x))})
        }else if(featureType=="op"){
            p <- lapply(genomicElements.GR, function(x){as.big.matrix(percent_func(GRanges(seqnames=tolower(CHR),
                                                                                           IRanges(start=seqDataTest,
                                                                                                   end=seqDataTest)),
                                                                                   x))})
        }
    }
    test_data <- big.matrix(nrow = length(seqDataTest),
                            ncol = length(genomicElements.GR),
                            dimnames = list(NULL,
                                            names(genomicElements.GR)))
    for(i in 1:ncol(test_data)){
        test_data[,i] <- p[[i]][,1]
    }

    rm("p")

    ##################################
    #PREDICTING AT BP RESOLUTION     #
    ##################################

    if(parallel==TRUE){
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
        predictions <- as.big.matrix(parallel_predictions(fit=tadModel,testing=as.matrix(test_data),c=cores,n=splits))
    }else{predictions <- as.big.matrix(predict(tadModel,newdata=as.matrix(test_data),type="prob")[,"Yes"])}

    if(threshold=="roc"){
        test_data_Y <- ifelse(seqDataTest %in% start(bounds.GR), 1, 0)
        t <- pROC::coords(pROC::roc(test_data_Y, predictions[,1], quiet=TRUE),
                          "best",
                          ret="threshold",
                          transpose = FALSE,
                          best.method="closest.topleft")
    }else{
        t <- threshold
    }

    if(verbose==TRUE){
        print(paste0("preciseTAD identified a total of ", length(diff(seqDataTest[which(predictions[,1]>=t)])), " base pairs whose predictive probability was equal to or exceeded a threshold of ", t))
    }

    retain <- c(1,cumsum(ifelse(diff(seqDataTest[which(predictions[,1]>=t)])!=1, 1, 0))+1)

    if(verbose==TRUE){
        print(paste0("preciseTAD identified ", length(unique(retain)), " PTBAs"));
        print("Establishing PTBRs")
    }

    mid <- tapply(seqDataTest[which(predictions[,1]>=t)], retain, function(x){
        return(ceiling((x[1]+x[length(x)])/2))
    })

    mid <- as.vector(mid)

    if(DBSCAN==FALSE){
        dist_mat <- dist(mid, method = "euclidean")
        if(verbose==TRUE){print("Initializing hierarchical clustering")}
        hc1 <- hclust(dist_mat, method = method.Clust)
        k <- pbsapply(2:(length(mid)-1), function(i) {
            mean(silhouette(cutree(hc1, i), dist=dist_mat)[,"sil_width"]) })
        k=which.max(k)+1
        if(verbose==TRUE){print(paste0("preciseTAD identified ", k, " PTBRs"))}
    }else{
        if(verbose==TRUE){print("Initializing DBSCAN")}
        res <- dbscan::dbscan(as.matrix(mid), eps=DBSCAN_params[[1]], minPts = DBSCAN_params[[2]])
        if(0 %in% unique(res$cluster)){
            k=length(unique(res$cluster))-1
        }else{k=length(unique(res$cluster))}
        if(verbose==TRUE){print(paste0("preciseTAD identified ", k, " PTBRs"))}
    }

    if(CLARA==TRUE){
        if(verbose==TRUE){
            print(paste0("Initializing CLARA with ", k, " clusters"));
            set.seed(seed)
            c <- clara(mid,
                       k=k,
                       samples=samples,
                       metric = method.Dist,
                       stand = FALSE,
                       trace = 2,
                       medoids.x = TRUE,
                       rngR = TRUE)
            medoids <- c$medoids
            clustering <- c$clustering
        }else{
            set.seed(seed)
            c <- clara(mid,
                       k=k,
                       samples=samples,
                       metric = method.Dist,
                       stand = FALSE,
                       trace = 0,
                       medoids.x = TRUE,
                       keep.data = FALSE,
                       rngR = TRUE)
            medoids <- c$medoids
            clustering <- c$clustering
        }
    }else{
        if(verbose==TRUE){
            print(paste0("Initializing PAM with ", k, " clusters"));
            set.seed(seed)
            c <- pam(x=mid, k=k, diss = FALSE,
                     metric = method.Dist,
                     medoids = NULL,
                     stand = FALSE,
                     cluster.only = FALSE,
                     do.swap = TRUE,
                     keep.diss = FALSE,
                     trace.lev = 2)
            medoids <- c$medoids
            clustering <- c$clustering
        }else{
            set.seed(seed)
            c <- pam(x=mid, k=k, diss = FALSE,
                     metric = method.Dist,
                     medoids = NULL,
                     stand = FALSE,
                     cluster.only = FALSE,
                     do.swap = TRUE,
                     keep.diss = FALSE,
                     trace.lev = 0)
            medoids <- c$medoids
            clustering <- c$clustering
        }
    }

    predBound_gr <- GRanges(seqnames=tolower(CHR),
                            IRanges(start=seqDataTest[which(seqDataTest %in% medoids)],
                                    end=seqDataTest[which(seqDataTest %in% medoids)]))

    trueBound_gr <- GRanges(seqnames=tolower(CHR),
                            IRanges(start=seqDataTest[which(seqDataTest %in% start(bounds.GR))],
                                    end=seqDataTest[which(seqDataTest %in% start(bounds.GR))]))

    if(!is.null(flank)){
        predBound_gr <- flank(predBound_gr, width=flank, both=TRUE)
        trueBound_gr <- flank(trueBound_gr, width=flank, both=TRUE)
    }

    juicer_func <- function(grdat){
        #n <- length(unique(as.character(seqnames(grdat))))

        chrs <- unique(as.character(seqnames(grdat)))

        mat_list <- list()

        for(i in 1:length(chrs)){
            grdat_chr <- grdat[which(as.character(seqnames(grdat))==chrs[i])]

            if(length(grdat_chr)>1){
                mymat <- matrix(nrow=(length(grdat_chr)-1),
                                ncol=6)
                mymat[,1] <- mymat[,4] <- chrs[i]
                for(j in 1:(length(grdat_chr)-1)){
                    mymat[j,2] <- mymat[j,5] <- start(grdat_chr)[j]
                    mymat[j,3] <- mymat[j,6] <- start(grdat_chr)[j+1]
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

        mat_list <- do.call("rbind",mat_list)
        return(mat_list)

    }

    if(PTBR==TRUE){
        grlist <- GRangesList()
        for(i in 1:length(unique(clustering))){

            PTBAs = c(1:length(clustering))[which(clustering==i)][1]
            PTBAe = c(1:length(clustering))[which(clustering==i)][length(which(clustering==i))]

            PTBRs = mid[PTBAs]
            PTBRe = mid[PTBAe]

            grlist[[i]] <- GRanges(seqnames = tolower(CHR),
                                   IRanges(start=PTBRs,
                                           end=PTBRe))
        }
        gr <- unlist(grlist)

        if(juicer==TRUE){
            bp_results <- list(PTBR=gr,
                               PTBP=juicer_func(predBound_gr),
                               CTBP=trueBound_gr)
        }else{
            bp_results <- list(PTBR=gr,
                               PTBP=predBound_gr,
                               CTBP=trueBound_gr)
        }
    }else{
        if(juicer==TRUE){
            bp_results <- list(PTBR=NULL,
                               PTBP=juicer_func(predBound_gr),
                               CTBP=trueBound_gr)
        }else{
            bp_results <- list(PTBR=NULL,
                               PTBP=predBound_gr,
                               CTBP=trueBound_gr)
        }
    }

    return(bp_results)
}
