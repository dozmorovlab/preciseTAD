#' Helper function for transforming a GRanges object into matrix form to be
#' saved as .txt or .BED file and imported into juicer
#'
#' @param grdat A GRanges object representing boundary coordinates
#'
#' @return A dataframe that can be saved as a BED file to import into juicer
#'
#' @export
#'
#' @import IRanges GenomicRanges
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
#' pt <- preciseTAD(genomicElements.GR = tfbsList,
#'                  featureType = "distance",
#'                  CHR = "CHR22",
#'                  chromCoords = list(17000000, 19000000),
#'                  tadModel = tadModel[[1]],
#'                  threshold = .75,
#'                  flank = NULL,
#'                  verbose = TRUE,
#'                  parallel = TRUE,
#'                  cores = 2,
#'                  splits = 2,
#'                  DBSCAN_params = list(5000, 3),
#'                  samples = 100)
#'
#' # Transform into juicer format
#' juicer_func(pt[[2]])
#' }
juicer_func <- function(grdat) {
    # n <- length(unique(as.character(seqnames(grdat))))

    chrs <- unique(as.character(seqnames(grdat)))

    mat_list <- list()

    for (i in seq_len(length(chrs))) {
        grdat_chr <- grdat[which(as.character(seqnames(grdat)) == chrs[i])]

        if (length(grdat_chr) > 1) {
            mymat <- data.frame(matrix(nrow = (length(grdat_chr) - 1), ncol = 6))
            mymat[, 1] <- mymat[, 4] <- chrs[i]
            for (j in seq_len(length(grdat_chr) - 1)) {
                mymat[j, 2] <- mymat[j, 5] <- start(grdat_chr)[j]
                mymat[j, 3] <- mymat[j, 6] <- start(grdat_chr)[j + 1]
            }

            mat_list[[i]] <- mymat

        }

    }

    mat_list <- do.call("rbind.data.frame", mat_list)
    return(mat_list)

}
