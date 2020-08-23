#' Helper function for transforming a GRanges object into matrix form to be
#' saved as .txt or .BED file and imported into juicer
#'
#' @param grdat A GRanges object
#'
#' @return A matrix that can be saved as a BED file to import into juicer
#'
#' @import IRanges GenomicRanges
#'
#' @examples
juicer_func <- function(grdat) {
    # n <- length(unique(as.character(seqnames(grdat))))

    chrs <- unique(as.character(seqnames(grdat)))

    mat_list <- list()

    for (i in seq_len(length(chrs))) {
        grdat_chr <- grdat[which(as.character(seqnames(grdat)) == chrs[i])]

        if (length(grdat_chr) > 1) {
            mymat <- matrix(nrow = (length(grdat_chr) - 1), ncol = 6)
            mymat[, 1] <- mymat[, 4] <- chrs[i]
            for (j in seq_len(length(grdat_chr) - 1)) {
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
