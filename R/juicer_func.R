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

juicer_func <- function(grdat) {
    # n <- length(unique(as.character(seqnames(grdat))))
    
    chrs <- unique(as.character(seqnames(grdat)))
    
    mat_list <- list()
    
    for (i in seq_len(length(chrs))) {
        grdat_chr <- grdat[which(as.character(seqnames(grdat)) == chrs[i])]
        
        if (length(grdat_chr) > 1) {
            mymat <- data.frame(matrix(nrow = (length(grdat_chr) - 1), ncol = 6))
            mymat[, 1] <- mymat[, 4] <- gsub("chr", "", chrs[i])
            for (j in seq_len(length(grdat_chr) - 1)) {
                mymat[j, 2] <- mymat[j, 5] <- start(grdat_chr)[j]
                mymat[j, 3] <- mymat[j, 6] <- start(grdat_chr)[j + 1]
            }
            
            mat_list[[i]] <- mymat
            
        }
        
    }
    
    mat_list <- do.call("rbind.data.frame", mat_list)
    names(mat_list) <- c("chr1", "x1", "x2", "chr2", "y1", "y2")
    return(mat_list)
    
}
