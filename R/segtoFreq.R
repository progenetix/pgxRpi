#' Calculate CNV frequency data from given segment data
#'
#' Thie function calculates the frequency of deletions and duplications
#'
#' @param data Segment data with CNV states. The first four columns should specify sample ID, chromosome, start position, and end position, respectively. The column representing CNV states should contain either "DUP" for duplications or "DEL" for deletions.
#' @param cnv_column_idx Index of the column specifying CNV state. Default is 6, following the "pgxseg" format used in Progenetix. 
#' If the input segment data uses the general `.seg` file format, it might need to be set differently.
#' @param cohort_name A string specifying the cohort name. Default is "unspecified cohort".
#' @param assembly A string specifying the genome assembly version for CNV frequency calculation. Allowed options are "hg19" or "hg38". Default is "hg38".
#' @param bin_size Size of genomic bins used to split the genome, in base pairs (bp). Default is 1,000,000.
#' @param overlap Numeric value defining the amount of overlap between bins and segments considered as bin-specific CNV, in base pairs (bp). Default is 1,000.
#' @param soft_expansion Fraction of `bin_size` to determine merge criteria.
#' During the generation of genomic bins, division starts at the centromere and expands towards the telomeres on both sides.
#' If the size of the last bin is smaller than `soft_expansion` * bin_size, it will be merged with the previous bin. Default is 0.1.
#' @return The binned CNV frequency stored in "pgxfreq" format
#' @export
#' @examples
#' ## load necessary data (this step can be skipped in real implementation)
#' data("hg38_cytoband")
#' ## get pgxseg data
#' seg <- read.table(system.file("extdata", "example.pgxseg",package = 'pgxRpi'),header=TRUE)
#' ## calculate frequency data
#' freq <- segtoFreq(seg)
#' ## visualize
#' pgxFreqplot(freq)

segtoFreq <- function(data,cnv_column_idx = 6,cohort_name = "unspecified cohort",assembly="hg38",bin_size= 1000000,overlap=1000,soft_expansion = 0.1){
    if (is.numeric(cnv_column_idx)) data <- data[,c(1,2,3,4,5,cnv_column_idx)]
	if (!all(unique(data[,6]) %in% c("DUP","DEL"))) stop ("calling states are invalid")
	
    bin.data <- extract.bin.feature(data,genome=assembly,bin.size=bin_size,overlap=overlap,soft.expansion = soft_expansion)
    bin.dup.data <- bin.data[['dup']]
    bin.del.data <- bin.data[['del']]
    nsamples <- length(unique(data[,1]))
    if (nsamples == 1){
        gain.freq <- bin.dup.data*100
        loss.freq <- bin.del.data*100
    } else{
        gain.freq <- (colSums(bin.dup.data)/nsamples)*100
        loss.freq <- (colSums(bin.del.data)/nsamples)*100
    }
  
    freq.seg <- cbind(group=cohort_name,bin.data[['bin']][,c(2,3,4)])
    freq.seg$gain_frequency <- as.numeric(gain.freq)
    freq.seg$loss_frequency <- as.numeric(loss.freq)
    frequency <- GenomicRanges::makeGRangesListFromDataFrame(freq.seg,split.field = "group",keep.extra.columns=TRUE)
    S4Vectors::mcols(frequency) <- data.frame(filter = "", label="", sample_count=nsamples)
    return(frequency)
}