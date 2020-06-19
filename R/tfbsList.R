#' A list of the chromosomal coordinates for 26 transcription factor binding
#' sites from the Gm12878 cell line
#'
#' A GRangesList containing 26 GRanges objects each with the following elements
#' \describe{
#' \item{seqnames}{The chromosome number}
#' \item{ranges}{IRanges object with start and end coordinates for each TFBS}
#' \item{strand}{empty column}
#' \item{coverage}{a metadata column with peak strength values at each site}
#' }
#'
#' @source Data obtained from ENCODE.
#' ENCODE integrative analysis (PMID: 22955616; PMCID: PMC3439153)
#' ENCODE portal (PMID: 29126249; PMCID: PMC5753278)
#' Available at
#' \url{https://www.encodeproject.org/}
#' @return A GrangeList
#'
"tfbsList"
