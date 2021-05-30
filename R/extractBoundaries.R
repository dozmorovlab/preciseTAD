#' Function to extract boundaries from domain data.
#'
#' @param domains.mat either a \code{matrix} or \code{data.frame} with at least
#' 3 columns. First column is chromosome number/character (1, 2, 3, X) or ID 
#' (chr1, chr2).  Non-autosomal/sex chromosomes will be filtered.
#' The second and third columns are the start and end coordinates 
#' of the domains, respectively. Note these are coordinates of the domain anchor
#' centers, not anchors. Only the first three columns are used. Required.
#' @param filter logical, indicating whether or not domains exceeding 2mb
#' in width or smaller than 2*(the specified resolution) should be filtered out
#' (default is FALSE, all boundaries will be used). Required.
#' @param CHR character vector, specifying which chromosome(s) to extract domain
#' boundaries on (ex: "chr22", case ignored). Unused seqnames are dropped. Required.
#' @param resolution numeric, the Hi-C data resolution that domains were called
#' at. Ignored if filter is FALSE, required otherwise.
#'
#' @return A \code{GRanges} object
#' @export
#'
#' @import IRanges GenomicRanges 
#' @importFrom gtools mixedorder
#' @examples
#' #Read in domain data from ARROWHEAD at 5 kb for GM12878
#' data("arrowhead_gm12878_5kb")
#' #Extract unique boundaries for CHRs 1-8 and 10-22
#' bounds.GR <- extractBoundaries(domains.mat=arrowhead_gm12878_5kb,
#'                                filter=FALSE,
#'                                CHR=paste0("CHR",c(1:8,10:22)),
#'                                resolution=5000)
extractBoundaries <- function(domains.mat, filter = FALSE, CHR, resolution) {
    
    resolution <- as.integer(resolution)
    #EXTRACTING UNIQUE BOUNDARIES#
    domains.mat <- domains.mat[, c(1,2,3)]
    # Filter out non-autosomal chromosomes that contain suffixes and underscores
    domains.mat <- domains.mat[!grepl("PATCH|CTG|ALT|GL|KI|UN|RND|_", domains.mat[, 1], ignore.case = TRUE, perl = TRUE), ]
    colnames(domains.mat) <- c("Chromosome", "Start", "End")
    domains.mat$Chromosome <- as.character(domains.mat$Chromosome)
    # Check if chromosomes have "chr" prefix, first three characters
    chr_prefix <- unique(substr(domains.mat$Chromosome, 1, 3))
    if (all(chr_prefix == "chr")) { 
        use_chr_prefix <- FALSE 
    } else { 
        use_chr_prefix <- TRUE 
    }
    
    if (filter == TRUE) {
        ## restricting domain data to TADs > 2*resolution and < 2,000,000
        domains.mat <- domains.mat[which((domains.mat$End - domains.mat$Start) >
                                             (2 * resolution) & (domains.mat$End - domains.mat$Start) < 2e+06), ]
    }
    
    ## concatenating boundary coordinates into one long vector
    coords <- domains.mat
    colnames(coords)[2:3] <- c("coordinate", "coordinate")
    coords <- rbind.data.frame(coords[, c(1, 2)], coords[, c(1, 3)])
    ## create a data frame for future GRanges    
    coords <- data.frame(chrom = paste0(ifelse(use_chr_prefix, "chr", ""), coords$Chromosome), start = coords$coordinate, end = coords$coordinate + 1)
    
    ## sorting by coordinate (smallest to largest), then chromosome (natural sort)
    coords <- coords[order(coords$start), ]
    coords <- coords[gtools::mixedorder(coords$chrom, decreasing = FALSE), ]
    ## remove duplicates for coordinates that are conjoined
    coords <- coords[!duplicated(coords), ]
    
    bounds <- GenomicRanges::makeGRangesFromDataFrame(coords)
    
    # chromsome specific boundary data
    bounds_chr_specific <- bounds[as.character(bounds@seqnames) %in% tolower(CHR)]

    return(bounds_chr_specific)
}
