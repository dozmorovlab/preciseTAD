#' Function to extract using boundaries from domain data
#'
#' @param domains.mat either a \code{matrix} or \code{data.frame} with at least
#' 3 columns. First column is either of class "numeric" or "integer". The second
#' and third columns are the start and end coordinates of domains, respectively.
#' Required.
#' @param preprocess logical, indicating whether or not domains exceeding 2mb
#' in width or smaller than 2*(the specified resolution) should be filtered out
#' (default is FALSE, all boundaries will be used). Required.
#' @param CHR character, specifying which chromosome(s) to extract domain
#' boundaries on (ex: "CHR22"). Required.
#' @param resolution numeric, the Hi-C data resolution that domains were called
#' at. Ignored if preprocess is FALSE.
#'
#' @return A \code{GRanges} object
#' @export
#'
#' @import IRanges GenomicRanges
#' @examples
#' #Read in domain data from ARROWHEAD at 5 kb for GM12878
#' data("arrowhead_gm12878_5kb")
#' #Extract unique boundaries for CHRs 1-8 and 10-22
#' bounds.GR <- extractBoundaries(domains.mat=arrowhead_gm12878_5kb,
#'                                preprocess=FALSE,
#'                                CHR=paste0("CHR",c(1:8,10:22)),
#'                                resolution=5000)
extractBoundaries <- function(domains.mat, preprocess = FALSE, CHR, resolution) {

    #EXTRACTING UNIQUE BOUNDARIES#

    resolution <- as.integer(resolution)
    domains.mat <- domains.mat[, c(1,2,3)]
    domains.mat[, 1] <- paste0("chr", domains.mat[, 1])
    colnames(domains.mat) <- c("Chromosome", "Start", "End")
    if (preprocess == TRUE) {
        ## restricting domain data to TADs > 2*resolution and < 2,000,000
        domains.mat <- domains.mat[which((domains.mat$End - domains.mat$Start) >
                                             (2 * resolution) & (domains.mat$End - domains.mat$Start) < 2e+06), ]
    }

    ## concatenating boundary coordinates into one long vector
    coords <- domains.mat
    colnames(coords)[2:3] <- c("coordinate", "coordinate")
    coords <- rbind.data.frame(coords[, c(1, 2)], coords[, c(1, 3)])
    coords$Chromosome <- as.character(coords$Chromosome)

    ## sorting by chromosome and coordinate
    coords <- coords[order(as.numeric(substr(coords$Chromosome, 4, 5)), coords$coordinate,
                           decreasing = FALSE), ]

    ## remove duplicates for coordinates that are conjoined
    coords <- coords[!duplicated(coords), ]

    bounds <- GRanges(seqnames = coords$Chromosome, ranges = IRanges(start = coords$coordinate,
                                                                     width = 1))

    # chromsome specific boundary data
    bounds_chr_specific <- bounds[which(as.character(seqnames(bounds)) %in% tolower(CHR))]

    return(bounds_chr_specific)
}
