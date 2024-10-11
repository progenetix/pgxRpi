#' Calculate CNV frequency data from given segment data
#'
#' Thie function calculates the frequency of deletions and duplications
#'
#' @param data Segment data containing CNV states. The first four columns should represent sample ID, chromosome, start position, and end position, respectively. 
#' The fifth column can contain the number of markers or other relevant information. 
#' The column representing CNV states (with a column index of 6 or higher) should either contain "DUP" for duplications and "DEL" for deletions, 
#' or level-specific CNV states such as "EFO:0030072", "EFO:0030071", "EFO:0020073", and "EFO:0030068", which correspond to high-level duplication, low-level duplication, high-level deletion, and low-level deletion, respectively.
#' @param cnv_column_idx Index of the column specifying the CNV state. Default is 6, based on the "pgxseg" format used in Progenetix. 
#' If the input segment data follows the general `.seg` file format, this index may need to be adjusted accordingly.
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
    if (!is(cnv_column_idx,"numeric") | cnv_column_idx < 6) stop("The parameter cnv_column_idx is invalid")
    data <- data[,c(1,2,3,4,5,cnv_column_idx)]
    if (all(grepl("EFO",data[,6]))){
      efo_mapping <- data.frame(efo_code = c("EFO:0030068","EFO:0030071","EFO:0030072","EFO:0020073","EFO:0030070","EFO:0030067"),
                                cnv_state = c("-1","+1","+2","-2","DUP","DEL"))
      data[,6] <- efo_mapping$cnv_state[match(data[,6],efo_mapping$efo_code)]
    }
	  if (!all(unique(data[,6]) %in% c("DUP","DEL")) & !all(unique(data[,6]) %in% c("+2","+1","-2","-1"))) stop ("CNV states are invalid")

    bin.data <- extract.bin.feature(data,genome=assembly,bin.size=bin_size,overlap=overlap,soft.expansion = soft_expansion)
    bin.dup.data1 <- bin.data[['dup1']]
    bin.del.data1 <- bin.data[['del1']]
    level_specific <- all(unique(data[,6]) %in% c("+2","+1","-2","-1"))
    if (level_specific){
        bin.dup.data2 <- bin.data[['dup2']]
        bin.del.data2 <- bin.data[['del2']]
    }
    nsamples <- length(unique(data[,1]))
    if (nsamples == 1){
        gain.freq1 <- bin.dup.data1*100
        loss.freq1 <- bin.del.data1*100
        if (level_specific){
            gain.freq2 <-bin.dup.data2*100
            loss.freq2 <-bin.del.data2*100
        }
    } else{
        gain.freq1 <- (colSums(bin.dup.data1)/nsamples)*100
        loss.freq1 <- (colSums(bin.del.data1)/nsamples)*100
        if (level_specific){
            gain.freq2 <- (colSums(bin.dup.data2)/nsamples)*100
            loss.freq2 <- (colSums(bin.del.data2)/nsamples)*100
        }
    }
  
    freq.seg <- cbind(group=cohort_name,bin.data[['bin']][,c(2,3,4)])

    if (level_specific){
        freq.seg$low_gain_frequency <- as.numeric(gain.freq1)
        freq.seg$low_loss_frequency <- as.numeric(loss.freq1)
        freq.seg$high_gain_frequency <- as.numeric(gain.freq2)
        freq.seg$high_loss_frequency <- as.numeric(loss.freq2)
    } else{
        freq.seg$gain_frequency <- as.numeric(gain.freq1)
        freq.seg$loss_frequency <- as.numeric(loss.freq1)
    }
    frequency <- GenomicRanges::makeGRangesListFromDataFrame(freq.seg,split.field = "group",keep.extra.columns=TRUE)
    S4Vectors::mcols(frequency) <- data.frame(filter = "", label="", sample_count=nsamples)
    return(frequency)
}