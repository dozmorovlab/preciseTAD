#' Function to create a data matrix used for building a predictive model to
#' classify boundary regions from functional genomic elements
#'
#' @param bounds.GR a GRanges object with chromosomal coordinates of TAD
#' boundaries used to identify positive cases (can be obtained using
#' \code{\link{extractBoundaries}})
#' @param resolution Numeric, the width to bin the genome at, should match the
#' resolution that TADs were called at
#' @param genomicElements.GR a GRangesList object containing GRanges objects
#' for each ChIP-seq data to leverage in the random forest model (can be
#' obtained using the \code{\link{bedToGRangesList}})
#' @param featureType Character, controls how the feature space is constructed
#' (one of either "binary" (overlap yes/no), "oc" (overlap counts, the number
#' of overlaps), "op" (overlap percent, the percent of bin width covered by the
#' genomic annotation), or "distance" (log2-transformed distance from the center
#' of the nearest genomic annotation to the center of the bin); default is
#' "distance")
#' @param resampling Character, controls if and how the data should be
#' resampled to create balanced classes of boundary vs. nonboundary regions (one
#' of either "none" - no re-sampling, "ros" - Random Over-Sampling, "rus" -
#' Random Under-Sampling, or "smote" - Synthetic Minority Over-sampling
#' TEchnique)
#' @param trainCHR Character vector of chromosomes to use to build the binned
#' data matrix for training
#' @param predictCHR Character vector of chromosomes to use to build the binned
#' data matrix for testing. Default in NULL, indicating no test data is created.
#'  If trainCHR=predictCHR then a 7:3 split is created
#' @param seed Numeric for reproducibility of resampling
#'
#' @return A list object containing two data.frames: 1) the training data, 2)
#' the test data (only if predictCHR is not NULL, otherwise it is NA). "y" is
#' an indicator whether the corresponding bin is a TAD boundary, and the
#' subsequent columns have the association measures between bins and the
#' genomic annotations
#' @export
#'
#
#' @import IRanges GenomicRanges DMwR
#'
#' @examples
#' # Create training data for CHR21 and testing data for CHR22 with
#' # 5 kb binning, oc-type predictors from 26 different transcription factor
#' # binding sites from the GM12878 cell line, and random under-sampling
#'
#' # Read in ARROWHEAD-called TADs at 5kb
#' data(arrowhead_gm12878_5kb)
#'
#' #Extract unique boundaries
#' bounds.GR <- extractBoundaries(domains.mat = arrowhead_gm12878_5kb,
#'                                preprocess = FALSE,
#'                                CHR = c("CHR21", "CHR22"),
#'                                resolution = 5000)
#'
#' # Read in GRangesList of 26 TFBS
#' data(tfbsList)
#'
#' tadData <- createTADdata(bounds.GR = bounds.GR,
#'                          resolution = 5000,
#'                          genomicElements.GR = tfbsList,
#'                          featureType = "oc",
#'                          resampling = "rus",
#'                          trainCHR = "CHR21",
#'                          predictCHR = "CHR22",
#'                          seed = 123)
createTADdata <- function(bounds.GR, resolution, genomicElements.GR, featureType = "distance",
                          resampling, trainCHR, predictCHR = NULL, seed = 123) {

    #CHECKING FUNCTION INPUTS#

    if (class(bounds.GR) != "GRanges") {
        print("is not a GRanges object!")
        return(0)
    }
    if (class(genomicElements.GR) != "CompressedGRangesList") {
        print("genomicElements.GR is not a CompressedGRangesList object!")
        return(0)
    }
    for (i in 1:length(genomicElements.GR)) {
        if (class(genomicElements.GR[[i]]) != "GRanges") {
            print(paste0(i, "-th object of genomicElements.GR is not a GenomicRanges object!"))
            return(0)
        }
    }
    if (class(featureType) != "character") {
        print("featureType is not a character object!")
        return(0)
    }
    if (!(featureType %in% c("distance", "binary", "oc", "op", "signal"))) {
        print("featureType must be one of either 'distance', 'binary', 'oc', or 'op'!")
        return(0)
    }
    if (class(resampling) != "character") {
        print("resampling is not a character object")
        return(0)
    }
    if (!(resampling %in% c("none", "ros", "rus", "smote"))) {
        print("resampling must be one of either 'none','ros','rus', or 'smote'!")
        return(0)
    }
    if (class(resolution) != "numeric") {
        print("resolution is not a numeric object!")
        return(0)
    }

    if (class(trainCHR) != "character") {
        print("trainCHR is not a character object!")
        return(0)
    }
    for (i in 1:length(trainCHR)) {
        if (grepl("CHR", trainCHR[i]) != TRUE) {
            print(paste0(i, "-th chromosome for training is not in 'trainCHR' format!"))
            return(0)
        }
    }

    if (!(is.null(predictCHR)) & class(predictCHR) != "character") {
        print("predictCHR is not a character object!")
        return(0)
    }

    if (!(is.null(predictCHR))) {
        for (i in 1:length(predictCHR)) {
            if (grepl("CHR", predictCHR[i]) != TRUE) {
                print(paste0(i, "-th chromosome for training is not in 'predictCHR' format!"))
                return(0)
            }
        }
    }

    if (length(intersect(trainCHR, predictCHR)) > 0) {
        print("there is a CHR that you are attempting to predict on that you are also training on; a 7:3 split will be implemented on the training:testing data")
    }

    if (!(is.null(predictCHR))) {
        if (length(intersect(trainCHR, predictCHR)) == 0) {
            print(c("The following chromosomes will be used to train:", trainCHR))
            print(c("The following chromosomes will be used to test:", predictCHR))
        }
    } else {
        print(c("The following chromosomes will be used to train:", trainCHR))
    }

    resolution = as.integer(resolution)

    #ESTABLISHING CHROMOSOME-SPECIFIC SEQINFO#

    #LOADING CHROMOSOME LENGTHS#

    hg19 <- preciseTAD:::hg19

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

    #ESTABLISHING DATA MATRIX FOR MODELING#

    if (is.null(predictCHR)) {
        # you are not interested in performance therefore no holdout test data is created

        # building training set
        bounds.GR.flank <- flank(bounds.GR[which(as.character(seqnames(bounds.GR)) %in%
                                                     tolower(trainCHR))], width = (resolution/2), both = TRUE)

        data_mat_list <- list()
        outcome_list <- list()
        train_list <- list()

        for (i in 1:length(trainCHR)) {
            start = (resolution/2)
            chrLength = hg19$length[hg19$chrom %in% trainCHR][i]
            centromereStart <- as.integer(hg19$centromerStart[hg19$chrom==trainCHR[i]])
            centromereEnd <- as.integer(hg19$centromerEnd[hg19$chrom==trainCHR[i]])
            end = chrLength - (chrLength %% resolution) + resolution/2

            data_mat_list[[i]] <- matrix(nrow=length(seq(start,end-1,resolution)[-which(seq(start,end-1,resolution) %in% c(centromereStart:centromereEnd))]),
                                         ncol=length(genomicElements.GR))

            dat_mat_gr <- GRanges(seqnames = tolower(trainCHR[i]),
                                  IRanges(start = seq(start,end-1,resolution)[-which(seq(start,end-1,resolution) %in% c(centromereStart:centromereEnd))],
                                          width=resolution))

            if (featureType == "distance") {
                # for use of distance type features
                dat_mat_gr <- GRanges(seqnames = tolower(trainCHR[i]), IRanges(start = (start(dat_mat_gr) +
                                                                                            end(dat_mat_gr))/2, width = 1))
                for (k in 1:length(genomicElements.GR)) {
                    d <- distance_func(dat_mat_gr, genomicElements.GR[[k]])
                    data_mat_list[[i]][, k] <- d
                }
                for (h in 1:ncol(data_mat_list[[i]])) {
                    data_mat_list[[i]][, h] <- log(data_mat_list[[i]][, h] + 1, base = 2)
                }
            } else if (featureType == "binary") {
                for (k in 1:length(genomicElements.GR)) {
                    cb <- binary_func(dat_mat_gr, genomicElements.GR[[k]])
                    data_mat_list[[i]][, k] <- cb
                }
            } else if (featureType == "oc") {
                for (k in 1:length(genomicElements.GR)) {
                    co <- count_func(dat_mat_gr, genomicElements.GR[[k]])
                    data_mat_list[[i]][, k] <- co
                }
            } else if (featureType == "op") {
                for (k in 1:length(genomicElements.GR)) {
                    cp <- percent_func(dat_mat_gr, genomicElements.GR[[k]])
                    data_mat_list[[i]][, k] <- cp
                }
            } else {
                cs <- signal_func(dat_mat_gr, genomicElements.GR[[k]])
                data_mat_list[[i]][, k] <- cs
            }

            outcome_list[[i]] <- countOverlaps(GRanges(seqnames = tolower(trainCHR[i]),
                                                       IRanges(start = seq(start,end-1,resolution)[-which(seq(start,end-1,resolution) %in% c(centromereStart:centromereEnd))],
                                                               width = resolution)), bounds.GR.flank)
            outcome_list[[i]] <- ifelse(outcome_list[[i]] >= 1, 1, 0)

            train_list[[i]] <- cbind.data.frame(outcome_list[[i]], as.matrix(data_mat_list[[i]]))
            names(train_list[[i]]) <- c("y", names(genomicElements.GR))
            train_list[[i]]$y <- factor(train_list[[i]]$y)
            levels(train_list[[i]]$y) <- c("No", "Yes")

            if (resampling == "ros") {
                # assign sample indeces
                set.seed(seed)
                sampids.train <- sample(x = which(train_list[[i]]$y == "Yes"), size = length(which(train_list[[i]]$y ==
                                                                                                       "No")), replace = TRUE)

                train_list[[i]] <- rbind.data.frame(train_list[[i]][which(train_list[[i]]$y ==
                                                                              "No"), ], train_list[[i]][sampids.train, ])

                # Randomly shuffle the data
                set.seed(seed)
                train_list[[i]] <- train_list[[i]][sample(nrow(train_list[[i]])),
                ]

            } else if (resampling == "rus") {
                set.seed(seed)
                sampids.train <- sample(x = which(train_list[[i]]$y == "No"), size = length(which(train_list[[i]]$y ==
                                                                                                      "Yes")), replace = FALSE)

                train_list[[i]] <- rbind.data.frame(train_list[[i]][which(train_list[[i]]$y ==
                                                                              "Yes"), ], train_list[[i]][sampids.train, ])

                # Randomly shuffle the data
                set.seed(seed)
                train_list[[i]] <- train_list[[i]][sample(nrow(train_list[[i]])),
                ]

            } else if (resampling == "smote") {
                set.seed(seed)
                train_list[[i]] <- SMOTE(y ~ ., data = train_list[[i]], perc.over = 100,
                                         perc.under = 200)

                # Randomly shuffle the data
                set.seed(seed)
                train_list[[i]] <- train_list[[i]][sample(nrow(train_list[[i]])),
                ]
            } else {
                train_list[[i]] = train_list[[i]]
            }


        }

        train_list <- do.call("rbind.data.frame", train_list)

    } else if (isTRUE(all.equal(sort(trainCHR), sort(predictCHR)))) {
        # you are training on the same chr(s) that you are predicting on so we build the
        # full datamatrix (including predictor type) and then split the data into
        # training and testing then perform resampling on training set

        bounds.GR.flank <- flank(bounds.GR, width = (resolution/2), both = TRUE)

        data_mat_list <- list()
        outcome_list <- list()
        full_list <- list()
        train_list <- list()
        test_list <- list()

        for (i in 1:length(trainCHR)) {
            start = (resolution/2)
            chrLength = hg19$length[hg19$chrom %in% trainCHR][i]
            centromereStart <- as.integer(hg19$centromerStart[hg19$chrom==trainCHR[i]])
            centromereEnd <- as.integer(hg19$centromerEnd[hg19$chrom==trainCHR[i]])
            end = chrLength - (chrLength %% resolution) + resolution/2

            data_mat_list[[i]] <- matrix(nrow=length(seq(start,end-1,resolution)[-which(seq(start,end-1,resolution) %in% c(centromereStart:centromereEnd))]),
                                         ncol=length(genomicElements.GR))

            dat_mat_gr <- GRanges(seqnames = tolower(trainCHR[i]),
                                  IRanges(start = seq(start,end-1,resolution)[-which(seq(start,end-1,resolution) %in% c(centromereStart:centromereEnd))],
                                          width=resolution))

            if (featureType == "distance") {
                # for use of distance type features
                dat_mat_gr <- GRanges(seqnames = tolower(trainCHR[i]), IRanges(start = (start(dat_mat_gr) +
                                                                                            end(dat_mat_gr))/2, width = 1))
                for (k in 1:length(genomicElements.GR)) {
                    d <- distance_func(dat_mat_gr, genomicElements.GR[[k]])
                    data_mat_list[[i]][, k] <- d
                }
                for (h in 1:ncol(data_mat_list[[i]])) {
                    data_mat_list[[i]][, h] <- log(data_mat_list[[i]][, h] + 1, base = 2)
                }
            } else if (featureType == "binary") {
                for (k in 1:length(genomicElements.GR)) {
                    cb <- binary_func(dat_mat_gr, genomicElements.GR[[k]])
                    data_mat_list[[i]][, k] <- cb
                }
            } else if (featureType == "oc") {
                for (k in 1:length(genomicElements.GR)) {
                    co <- count_func(dat_mat_gr, genomicElements.GR[[k]])
                    data_mat_list[[i]][, k] <- co
                }
            } else if (featureType == "op") {
                for (k in 1:length(genomicElements.GR)) {
                    cp <- percent_func(dat_mat_gr, genomicElements.GR[[k]])
                    data_mat_list[[i]][, k] <- cp
                }
            }

            outcome_list[[i]] <- countOverlaps(GRanges(seqnames = tolower(trainCHR[i]),
                                                       IRanges(start = seq(start,end-1,resolution)[-which(seq(start,end-1,resolution) %in% c(centromereStart:centromereEnd))],
                                                               width = resolution)), bounds.GR.flank)
            outcome_list[[i]] <- ifelse(outcome_list[[i]] >= 1, 1, 0)

            full_list[[i]] <- cbind.data.frame(outcome_list[[i]], as.matrix(data_mat_list[[i]]))

            names(full_list[[i]]) <- c("y", names(genomicElements.GR))
            full_list[[i]]$y <- factor(full_list[[i]]$y)
            levels(full_list[[i]]$y) <- c("No", "Yes")

            set.seed(seed)
            inTrainingSet <- sample(length(full_list[[i]]$y), floor(length(full_list[[i]]$y) *
                                                                        0.7))
            train_list[[i]] <- full_list[[i]][inTrainingSet, ]
            test_list[[i]] <- full_list[[i]][-inTrainingSet, ]

            if (resampling == "ros") {
                # assign sample indeces
                set.seed(seed)
                sampids.train <- sample(x = which(train_list[[i]]$y == "Yes"), size = length(which(train_list[[i]]$y ==
                                                                                                       "No")), replace = TRUE)

                train_list[[i]] <- rbind.data.frame(train_list[[i]][which(train_list[[i]]$y ==
                                                                              "No"), ], train_list[[i]][sampids.train, ])

                # Randomly shuffle the data
                set.seed(seed)
                train_list[[i]] <- train_list[[i]][sample(nrow(train_list[[i]])),
                ]

            } else if (resampling == "rus") {
                set.seed(seed)
                sampids.train <- sample(x = which(train_list[[i]]$y == "No"), size = length(which(train_list[[i]]$y ==
                                                                                                      "Yes")), replace = FALSE)

                train_list[[i]] <- rbind.data.frame(train_list[[i]][which(train_list[[i]]$y ==
                                                                              "Yes"), ], train_list[[i]][sampids.train, ])

                # Randomly shuffle the data
                set.seed(seed)
                train_list[[i]] <- train_list[[i]][sample(nrow(train_list[[i]])),
                ]

            } else if (resampling == "smote") {
                set.seed(seed)
                train_list[[i]] <- SMOTE(y ~ ., data = train_list[[i]], perc.over = 100,
                                         perc.under = 200)

                # Randomly shuffle the data
                set.seed(seed)
                train_list[[i]] <- train_list[[i]][sample(nrow(train_list[[i]])),
                ]
            } else {
                train_list[[i]] = train_list[[i]]
            }

        }

        train_list <- do.call("rbind.data.frame", train_list)
        test_list <- do.call("rbind.data.frame", test_list)
    } else {
        # you are training on one set of chr(s) and are predicting on another set of
        # chr(s) first build training set on the trainCHR, implement resampling technique
        # then build testing set from the predictCHR, with no resampling technique

        # building training set
        bounds.GR.flank <- flank(bounds.GR[which(as.character(seqnames(bounds.GR)) %in%
                                                     tolower(trainCHR))], width = (resolution/2), both = TRUE)

        data_mat_list <- list()
        outcome_list <- list()
        train_list <- list()

        for (i in 1:length(trainCHR)) {
            start = (resolution/2)
            chrLength = hg19$length[hg19$chrom %in% trainCHR][i]
            centromereStart <- as.integer(hg19$centromerStart[hg19$chrom==trainCHR[i]])
            centromereEnd <- as.integer(hg19$centromerEnd[hg19$chrom==trainCHR[i]])
            end = chrLength - (chrLength %% resolution) + resolution/2

            data_mat_list[[i]] <- matrix(nrow=length(seq(start,end-1,resolution)[-which(seq(start,end-1,resolution) %in% c(centromereStart:centromereEnd))]),
                                         ncol=length(genomicElements.GR))

            dat_mat_gr <- GRanges(seqnames = tolower(trainCHR[i]),
                                  IRanges(start = seq(start,end-1,resolution)[-which(seq(start,end-1,resolution) %in% c(centromereStart:centromereEnd))],
                                          width=resolution))

            if (featureType == "distance") {
                # for use of distance type features
                dat_mat_gr <- GRanges(seqnames = tolower(trainCHR[i]), IRanges(start = (start(dat_mat_gr) +
                                                                                            end(dat_mat_gr))/2, width = 1))
                for (k in 1:length(genomicElements.GR)) {
                    d <- distance_func(dat_mat_gr, genomicElements.GR[[k]])
                    data_mat_list[[i]][, k] <- d
                }
                for (h in 1:ncol(data_mat_list[[i]])) {
                    data_mat_list[[i]][, h] <- log(data_mat_list[[i]][, h] + 1, base = 2)
                }
            } else if (featureType == "binary") {
                for (k in 1:length(genomicElements.GR)) {
                    cb <- binary_func(dat_mat_gr, genomicElements.GR[[k]])
                    data_mat_list[[i]][, k] <- cb
                }
            } else if (featureType == "oc") {
                for (k in 1:length(genomicElements.GR)) {
                    co <- count_func(dat_mat_gr, genomicElements.GR[[k]])
                    data_mat_list[[i]][, k] <- co
                }
            } else if (featureType == "op") {
                for (k in 1:length(genomicElements.GR)) {
                    cp <- percent_func(dat_mat_gr, genomicElements.GR[[k]])
                    data_mat_list[[i]][, k] <- cp
                }
            }

            outcome_list[[i]] <- countOverlaps(GRanges(seqnames = tolower(trainCHR[i]),
                                                       IRanges(start = seq(start,end-1,resolution)[-which(seq(start,end-1,resolution) %in% c(centromereStart:centromereEnd))],
                                                               width = resolution)), bounds.GR.flank)

            outcome_list[[i]] <- ifelse(outcome_list[[i]] >= 1, 1, 0)

            train_list[[i]] <- cbind.data.frame(outcome_list[[i]], as.matrix(data_mat_list[[i]]))
            names(train_list[[i]]) <- c("y", names(genomicElements.GR))
            train_list[[i]]$y <- factor(train_list[[i]]$y)
            levels(train_list[[i]]$y) <- c("No", "Yes")

            if (resampling == "ros") {
                # assign sample indeces
                set.seed(seed)
                sampids.train <- sample(x = which(train_list[[i]]$y == "Yes"), size = length(which(train_list[[i]]$y ==
                                                                                                       "No")), replace = TRUE)

                train_list[[i]] <- rbind.data.frame(train_list[[i]][which(train_list[[i]]$y ==
                                                                              "No"), ], train_list[[i]][sampids.train, ])

                # Randomly shuffle the data
                set.seed(seed)
                train_list[[i]] <- train_list[[i]][sample(nrow(train_list[[i]])),
                ]

            } else if (resampling == "rus") {
                set.seed(seed)
                sampids.train <- sample(x = which(train_list[[i]]$y == "No"), size = length(which(train_list[[i]]$y ==
                                                                                                      "Yes")), replace = FALSE)

                train_list[[i]] <- rbind.data.frame(train_list[[i]][which(train_list[[i]]$y ==
                                                                              "Yes"), ], train_list[[i]][sampids.train, ])

                # Randomly shuffle the data
                set.seed(seed)
                train_list[[i]] <- train_list[[i]][sample(nrow(train_list[[i]])),
                ]

            } else if (resampling == "smote") {
                set.seed(seed)
                train_list[[i]] <- SMOTE(y ~ ., data = train_list[[i]], perc.over = 100,
                                         perc.under = 200)

                # Randomly shuffle the data
                set.seed(seed)
                train_list[[i]] <- train_list[[i]][sample(nrow(train_list[[i]])),
                ]
            } else {
                train_list[[i]] = train_list[[i]]
            }


        }

        train_list <- do.call("rbind.data.frame", train_list)

        # building testing set
        bounds.GR.flank.pred <- flank(bounds.GR[which(as.character(seqnames(bounds.GR)) %in%
                                                          tolower(predictCHR))], width = (resolution/2), both = TRUE)

        data_mat_list <- list()
        outcome_list <- list()
        test_list <- list()

        for (i in 1:length(predictCHR)) {
            start = (resolution/2)
            chrLength = hg19$length[hg19$chrom %in% predictCHR][i]
            centromereStart <- as.integer(hg19$centromerStart[hg19$chrom==predictCHR[i]])
            centromereEnd <- as.integer(hg19$centromerEnd[hg19$chrom==predictCHR[i]])
            end = chrLength - (chrLength %% resolution) + resolution/2

            data_mat_list[[i]] <- matrix(nrow=length(seq(start,end-1,resolution)[-which(seq(start,end-1,resolution) %in% c(centromereStart:centromereEnd))]),
                                         ncol=length(genomicElements.GR))

            dat_mat_gr <- GRanges(seqnames = tolower(predictCHR[i]),
                                  IRanges(start = seq(start,end-1,resolution)[-which(seq(start,end-1,resolution) %in% c(centromereStart:centromereEnd))],
                                          width=resolution))

            if (featureType == "distance") {
                # for use of distance type features
                dat_mat_gr <- GRanges(seqnames = tolower(predictCHR[i]), IRanges(start = (start(dat_mat_gr) +
                                                                                              end(dat_mat_gr))/2, width = 1))
                for (k in 1:length(genomicElements.GR)) {
                    d <- distance_func(dat_mat_gr, genomicElements.GR[[k]])
                    data_mat_list[[i]][, k] <- d
                }
                for (h in 1:ncol(data_mat_list[[i]])) {
                    data_mat_list[[i]][, h] <- log(data_mat_list[[i]][, h] + 1, base = 2)
                }
            } else if (featureType == "binary") {
                for (k in 1:length(genomicElements.GR)) {
                    cb <- binary_func(dat_mat_gr, genomicElements.GR[[k]])
                    data_mat_list[[i]][, k] <- cb
                }
            } else if (featureType == "oc") {
                for (k in 1:length(genomicElements.GR)) {
                    co <- count_func(dat_mat_gr, genomicElements.GR[[k]])
                    data_mat_list[[i]][, k] <- co
                }
            } else if (featureType == "op") {
                for (k in 1:length(genomicElements.GR)) {
                    cp <- percent_func(dat_mat_gr, genomicElements.GR[[k]])
                    data_mat_list[[i]][, k] <- cp
                }
            }

            outcome_list[[i]] <- countOverlaps(GRanges(seqnames = tolower(predictCHR[i]),
                                                       IRanges(start = seq(start,end-1,resolution)[-which(seq(start,end-1,resolution) %in% c(centromereStart:centromereEnd))],
                                                               width = resolution)), bounds.GR.flank.pred)
            outcome_list[[i]] <- ifelse(outcome_list[[i]] >= 1, 1, 0)

            test_list[[i]] <- cbind.data.frame(outcome_list[[i]], as.matrix(data_mat_list[[i]]))
            names(test_list[[i]]) <- c("y", names(genomicElements.GR))
            test_list[[i]]$y <- factor(test_list[[i]]$y)
            levels(test_list[[i]]$y) <- c("No", "Yes")

        }

        test_list <- do.call("rbind.data.frame", test_list)

    }

    if (is.null(predictCHR)) {
        TADdataList <- list(train_list, NA)
    } else {
        TADdataList <- list(train_list, test_list)
    }

    return(TADdataList)
}
