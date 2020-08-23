#' Helper function used to create count overlap type feature space
#'
#' @param binned_data_gr A GRanges object
#' @param annot_data_gr A Granges object
#'
#' @return A vector of counts enumerating the number of overlaps
#'
#' @import IRanges GenomicRanges
#'
#' @examples
count_func <- function(binned_data_gr, annot_data_gr) {

    # Finding the total number of overlaps between genomic bins and the specific
    # genomic annotation
    count_total <- countOverlaps(binned_data_gr, annot_data_gr)

    return(count_total)
}
