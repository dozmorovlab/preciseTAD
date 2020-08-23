#' Function to create a GRangesList object from functional genomic annotation
#' data derived from ChIP-Seq experiments and in the form of BED files
#'
#' @param filepath Character describing the path to the folder containing the
#' BED files of functional genomic annotations
#' @param pattern Character describing the pattern of the files for the
#' functional genomic annotations. Default is "*.bed"
#' @param signal Numeric referring to the column in the BED files that denotes
#' coverage strength. Must be the same for all files. Default is 5 (fifth
#' column), as is the case with most BED files
#'
#' @return A \code{GRangesList} object where each entry is a \code{GRanges}
#' object specific to each BED file in the path provided
#' @export
#'
#' @importFrom utils read.table
#' @import IRanges GenomicRanges
#'
#' @examples
#' #set path
#' path <- system.file("extdata", package = "preciseTAD")
#' #contains 2 BED files representing YY1 and NFYA
#' #transcription factor binding sites for GM12878
#' tfbsList <- bedToGRangesList(filepath = path, pattern = "*.bed", signal=4)
bedToGRangesList <- function(filepath, pattern = "*.bed", signal = 5) {

    #CREATING GRangesList#

    annots_gr_list <- GRangesList()
    # if(grepl(paste0(names, collapse = '|'), list.files(filepath)))
    for (i in seq_len(length(list.files(filepath, pattern = pattern)))) {
        dat <- read.table(paste0(filepath, "/", list.files(filepath, pattern = pattern)[i]),
                          header = FALSE)

        annots_gr_list[[i]] <- GRanges(seqnames = dat$V1, IRanges(start = dat$V2,
                                                                  end = dat$V3))

        if (!(is.null(signal))) {
            mcols(annots_gr_list[[i]])$coverage <- dat[, signal]
        }

    }
    names(annots_gr_list) <- gsub(pattern, "", list.files(filepath, pattern = pattern))

    return(annots_gr_list)
}
