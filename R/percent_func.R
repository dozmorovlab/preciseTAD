#' Helper function used to create percent overlap type feature space
#'
#' @param binned_data_gr A GRanges object
#' @param annot_data_gr A Granges object
#'
#' @return A vector of proportions indicating the percentage of overlap
#'
#' @import IRanges GenomicRanges
#'
#' @examples
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
