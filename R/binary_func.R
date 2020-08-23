#' Helper function used to create binary overlap type feature space
#'
#' @param binned_data_gr A GRanges object
#' @param annot_data_gr A Granges object
#'
#' @return An indicator vector for whether or not an overlap occured
#'
#' @import IRanges GenomicRanges
#'
#' @examples
binary_func <- function(binned_data_gr, annot_data_gr) {

    # Finding the total number of overlaps between genomic bins and the specific
    # genomic annotation
    count_binary <- countOverlaps(binned_data_gr, annot_data_gr)

    # Binarizing the overlap (1 or 0)
    count_binary <- ifelse(count_binary > 0, 1, 0)

    return(count_binary)
}
