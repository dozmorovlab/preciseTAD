#' Helper function used to create signal type feature space
#'
#' @param binned_data_gr A GRanges object
#' @param annot_data_gr A GRanges object
#'
#' @return A vector of intensities indicating the signal strength within each
#' overlap
#'
#' @import IRanges GenomicRanges
#'
#' @examples
signal_func <- function(binned_data_gr, annot_data_gr){

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
