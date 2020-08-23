#' Helper function used to create (log2) distance type feature space
#'
#' @param binned_data_center_gr A GRanges object of width 1
#' @param annot_data_center_gr A GRanges object of width 1
#'
#' @return A vector of log2 distances to the nearest overlap
#'
#' @import IRanges GenomicRanges
#'
distance_func <- function(binned_data_center_gr, annot_data_center_gr) {

    # distance from center of genomic bin to nearest genomic region of interest
    dist <- mcols(distanceToNearest(binned_data_center_gr, annot_data_center_gr))$distance

    return(dist)
}
