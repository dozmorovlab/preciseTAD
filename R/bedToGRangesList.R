#' Function to create a GRangesList object from functional genomic annotation
#' data derived from ChIP-Seq experiments and in the form of .BED files
#'
#' @param filepath Character describing the path to the folder containing the
#' .BED files of functional genomic annotations
#' @param pattern Character describing the pattern of the files for the
#' functional genomic annotations. Default is "*.bed"
#' @param signal Numeric refering to the column in the .BED files that denotes
#' the coverage strength. Must be the same for all files.
#'
#' @return A \code{GRangesList} object where each entry is a \code{GRanges}
#' object specific to each .BED file in the path provided
#' @export
#'
#' @importFrom utils read.table
#' @import IRanges GenomicRanges
#'
#' @examples
#' #set path
#' path = system.file("extdata", package = "preciseTAD")
#' #contains 2 .BED files representing YY1 and NFYA
#' #transcription factor binding sites for GM12878
#' tfbsList <- bedToGRangesList(filepath=path, pattern = "*.bed", signal=4)
bedToGRangesList <- function(filepath,
                             pattern = "*.bed",
                             signal){

    #############
    #STOP CHECKS#
    #############

    if(class(filepath)!="character"){
        print("filepath is not a character object!")
        return(0)
    }
    if(class(pattern)!="character"){
        print("pattern is not a character object!")
        return(0)
    }
    if(!is.null(signal)){
        if(class(signal)!="numeric"){
            print("signal is not a numeric object!")
            return(0)
        }
    }

    ######################
    #Creating GRangesList#
    ######################

    annots_gr_list <- GRangesList()
    #if(grepl(paste0(names, collapse = "|"), list.files(filepath)))
    for(i in 1:length(list.files(filepath, pattern = pattern))){
        dat <- read.table(paste0(filepath,
                                 "/",
                                 list.files(filepath,
                                            pattern = pattern)[i]),
                          header = FALSE)

        annots_gr_list[[i]] <- GRanges(seqnames=dat$V1,
                                       IRanges(start=dat$V2,
                                               end=dat$V3))

        if(!(is.null(signal))){
            mcols(annots_gr_list[[i]])$coverage <- dat[,signal]
        }

    }
    names(annots_gr_list) <- gsub(pattern, "", list.files(filepath, pattern = pattern))

    return(annots_gr_list)
}
