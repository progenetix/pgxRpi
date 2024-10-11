#' Plot CNV frequency data
#'
#' Thie function plots the frequency of deletions and duplications
#'
#' @param data CNV frequency object returned by the `pgxLoader` or `segtoFreq` functions.
#' @param chrom A vector specifying which chromosomes to plot. If NULL, the plot will cover the entire genome. 
#' If specified, the frequencies are plotted with one panel for each chromosome. Default is NULL.
#' @param layout Number of columns and rows in plot. Only used in plot by chromosome. Default is c(1,1).
#' @param filters Index or string value indicating which filter to plot. The length of filters
#' is limited to one if the parameter `circos` is FALSE. Default is the first filter.
#' @param circos A logical value indicating whether to return a circos plot. If TRUE, it returns a circos plot 
#' that can display and compare multiple filters. Default is FALSE.
#' @param assembly A string specifying the genome assembly version to apply to CNV frequency plotting. 
#' Allowed options are "hg19" and "hg38". Default is "hg38".
#' @return The binned CNV frequency plot
#' @importFrom grDevices dev.cur dev.new devAskNewPage
#' @importFrom graphics abline axis mtext par polygon rect text title
#' @importFrom utils data
#' @importFrom methods is
#' @export
#' @examples
#' ## load necessary data (this step can be skipped in real implementation)
#' data("hg38_cytoband")
#' ## get frequency data
#' freq <- pgxLoader(type="cnv_frequency", output ='pgxfreq', filters="NCIT:C3512")
#' ## visualize
#' pgxFreqplot(freq)

pgxFreqplot <- function(data,chrom=NULL,layout=c(1,1),filters=NULL,circos = FALSE,assembly = 'hg38'){
    if (circos){
        return(cplotpgxFreq(data,filters,assembly))
    }

    if (length(filters) > 1){
        stop("The length of filters exceeds limit")
    }

    type <- ifelse(is.null(chrom),"genome","bychrom")
    switch(type,
            genome = genomeFreq(data=data,id=filters,assembly = assembly),
            bychrom = chromosomeFreq(data=data,chrom=chrom,layout=layout,id=filters,assembly = assembly))
}


cplotpgxFreq  <- function(data,filters,assembly){
    if (is.null(filters)) filters <- 1
    data_list <- list()
    # input is pgxfreq
    if (is(data, "CompressedGRangesList")){
        meta <- S4Vectors::mcols(data)
        id_name <- rownames(meta[filters,])
        for (i in id_name){
            ind_data <- data[[i]]
            ind_data_freq <- S4Vectors::mcols(ind_data)
            if ("low_gain_frequency" %in% names(ind_data_freq)){
                ind_df <- data.frame(chr=paste0('chr',as.character(GenomicRanges::seqnames(ind_data))),
                                    start=GenomicRanges::start(ind_data),
                                    end=GenomicRanges::end(ind_data),
                                    lgain=ind_data_freq$low_gain_frequency,
                                    lloss=ind_data_freq$low_loss_frequency,
                                    hgain=ind_data_freq$high_gain_frequency,
                                    hloss=ind_data_freq$high_loss_frequency)
            } else if ("gain_frequency" %in% names(ind_data_freq)){
                ind_df <- data.frame(chr=paste0('chr',as.character(GenomicRanges::seqnames(ind_data))),
                                    start=GenomicRanges::start(ind_data),
                                    end=GenomicRanges::end(ind_data),
                                    lgain=ind_data_freq$gain_frequency,
                                    lloss=ind_data_freq$loss_frequency)
            } else{
                stop("Input is invalid")
            }
            data_list[[i]] <- ind_df
        } 
    } else if (is(data, "RangedSummarizedExperiment")){
        meta <- SummarizedExperiment::colData(data)
        id_name <- rownames(meta[filters,])
        range.info <- unique(GenomicRanges::granges(SummarizedExperiment::rowRanges(data)))
        for (i in id_name){
            ind_data_freq <- SummarizedExperiment::assays(data)
            if ("lowlevel_cnv_frequency" %in% names(ind_data_freq) & "highlevel_cnv_frequency" %in% names(ind_data_freq)){                
                ind_data_row <- SummarizedExperiment::rowData(data)
                gain.freq1 <- ind_data_freq$lowlevel_cnv_frequency[which(ind_data_row$type == "gain"),i]
                loss.freq1 <- ind_data_freq$lowlevel_cnv_frequency[which(ind_data_row$type == "loss"),i]
                gain.freq2 <- ind_data_freq$highlevel_cnv_frequency[which(ind_data_row$type == "gain"),i]
                loss.freq2 <- ind_data_freq$highlevel_cnv_frequency[which(ind_data_row$type == "loss"),i]
                ind_df <- data.frame(chr=paste0('chr',as.character(GenomicRanges::seqnames(range.info))),
                                    start=GenomicRanges::start(range.info),
                                    end=GenomicRanges::end(range.info),
                                    lgain=gain.freq1,
                                    lloss=loss.freq1,
                                    hgain=gain.freq2,
                                    hloss=loss.freq2)
            }else{
                stop("Input is invalid")
            }                   
            data_list[[i]] <- ind_df
        }
    } else{
        stop("Input is invalid")
    }
 
  # set text size and color
    cex_text <- ifelse(length(data_list) == 1, 0.8, 0.6)
    col.lowgain <-'#FFC633'
    col.lowloss <- '#33A0FF'
    col.highgain="#d00000"
    col.highloss="#0000FF"

      # plot
    for (i in seq_len(length(id_name))){
        id <- id_name[i]
        data_cir <- data_list[[id]]
        # Initialize
        if (i == 1){
            circlize::circos.clear()
            circlize::circos.par(start.degree=90, gap.degree=c(rep(1,23),10))
            circlize::circos.initializeWithIdeogram(species = assembly)
        }

        if (!is.null(data_cir$hgain) & !is.null(data_cir$hloss)){
            circlize::circos.genomicTrack(data_cir, ylim = c(-100, 100),
                                  numeric.column = c("lgain", "lloss","hgain","hloss"),
                                  panel.fun = function(region, value, ...){
                                      value1 <- unlist(value[1])
                                      value2 <- unlist(value[2])
                                      value3 <- unlist(value[3])
                                      value4 <- unlist(value[4])
                                      pos <- unlist((region[2]+region[1])/2)              
                                      circlize::circos.barplot(value=value1,pos=pos, border = col.lowgain, col = col.lowgain)
                                      circlize::circos.barplot(value=-value2,pos=pos, border = col.lowloss, col = col.lowloss)
                                      circlize::circos.barplot(value=value3,pos=pos, border = col.highgain, col = col.highgain)
                                      circlize::circos.barplot(value=-value4,pos=pos, border = col.highloss, col = col.highloss)
                                      y <- seq(-100,100,by=20)
                                      circlize::circos.segments(0,y,circlize::CELL_META$xlim[2], y,col=ifelse(y>0,col.lowgain,ifelse(y == 0,'#909090',col.lowloss)), lwd=0.5)
                                  })

        } else{
            circlize::circos.genomicTrack(data_cir, ylim = c(-100, 100),
                      numeric.column = c("lgain", "lloss"),
                      panel.fun = function(region, value, ...){
                          value1 <- unlist(value[1])
                          value2 <- unlist(value[2])
                          pos <- unlist((region[2]+region[1])/2)              
                          circlize::circos.barplot(value=value1,pos=pos, border = col.lowgain, col = col.lowgain)
                          circlize::circos.barplot(value=-value2,pos=pos, border = col.lowloss, col = col.lowloss)
                          y <- seq(-100,100,by=20)
                          circlize::circos.segments(0,y,circlize::CELL_META$xlim[2], y,col=ifelse(y>0,col.lowgain,ifelse(y == 0,'#909090',col.lowloss)), lwd=0.5)
                      })

        }

        if (i <= 2 & length(id_name) <= 2){text(0,-0.1*(i-1),id,cex = cex_text)}
    }

}
