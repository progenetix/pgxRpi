#' Plot CNV frequency data
#'
#' Thie function plots the frequency of deletions and duplications
#'
#' @param data The frequency object returned by `pgxLoader` function.
#' @param chrom A vector with chromosomes to be plotted. If NULL, return the plot
#' by genome. If specified the frequencies are plotted with one panel for each
#' chromosome. Default is NULL.
#' @param layout Number of columns and rows in plot. Only used in plot by chromosome.
#' Default is c(1,1).
#' @param filters Index or string value to indicate which filter to be plotted, such as 1
#' (the first filters in `data` slot of object ) or 'NCIT:C4038' (specific filter name). The length of filters
#' is limited to one if the parameter `circos` is False. Default is 1.
#' @param circos A logical value to indicate if return a circos plot. If TRUE, it
#' can return a circos plot with multiple filters for display and comparison.
#' Default is FALSE.
#' @param highlight Indices of genomic bins to be highlighted with red color.
#' @param assembly A string specifying which genome assembly version should be applied to CNV frequency plotting.
#' Allowed options are "hg19", "hg38". Default is "hg38" (genome version used in Progenetix).
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
#' freq <- pgxLoader(type="frequency", output ='pgxfreq', filters="NCIT:C3512")
#' ## visualize
#' pgxFreqplot(freq)

pgxFreqplot <- function(data,chrom=NULL,layout=c(1,1),filters=NULL,circos = FALSE,highlight=NULL,assembly = 'hg38'){
    if (circos){
        return(cplotpgxFreq(data,filters,highlight,assembly))
    }

    if (length(filters) > 1){
        stop("The length of filters exceeds limit")
    }

    type <- ifelse(is.null(chrom),"genome","bychrom")
    switch(type,
            genome = genomeFreq(data=data,id=filters,highlight=highlight,assembly = assembly),
            bychrom = chromosomeFreq(data=data,chrom=chrom,layout=layout,id=filters,highlight=highlight,assembly = assembly))
}


cplotpgxFreq  <- function(data,filters,highlight,assembly){
    if (is.null(filters)) filters <- 1
    data_list <- list()
    # input is pgxfreq
    if (is(data, "CompressedGRangesList")){
        meta <- S4Vectors::mcols(data)
        id_name <- rownames(meta[filters,])
        for (i in id_name){
          ind_data <- data[[i]]
          data_list[[i]] <- data.frame(chr=paste0('chr',as.character(GenomicRanges::seqnames(ind_data))),
                                       start=GenomicRanges::start(ind_data),
                                       end=GenomicRanges::end(ind_data),
                                       gain=S4Vectors::mcols(ind_data)$gain_frequency,
                                       loss=S4Vectors::mcols(ind_data)$loss_frequency)
        } 
    } else if (is(data, "RangedSummarizedExperiment")){
        meta <- SummarizedExperiment::colData(data)
        id_name <- rownames(meta[filters,])
        range.info <- unique(GenomicRanges::granges(SummarizedExperiment::rowRanges(data)))
        for (i in id_name){
            gain.freq <- SummarizedExperiment::assay(data)[which(SummarizedExperiment::rowData(data)$type == "DUP"),i]
            loss.freq <- SummarizedExperiment::assay(data)[which(SummarizedExperiment::rowData(data)$type == "DEL"),i]
            data_list[[i]] <-  data.frame(chr=paste0('chr',as.character(GenomicRanges::seqnames(range.info))),
                                             start=GenomicRanges::start(range.info),
                                             end=GenomicRanges::end(range.info),
                                             gain=gain.freq,
                                             loss=-loss.freq)
        }
    } else{
        stop("\n The input is invalid \n")
    }
 

  # Identify highlight location
    if (!is.null(highlight)){
        h_chr <- data_list[[1]][['chr']][highlight]
        h_start <- data_list[[1]][['start']][highlight]
    }


  # set text size and color
    cex_text <- ifelse(length(data_list) == 1, 0.8, 0.6)
    col.gain <-'#FFC633'
    col.loss <- '#33A0FF'
    bar.col.gain <- col.gain
    bar.col.loss <- col.loss

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
        circlize::circos.genomicTrack(data_cir, ylim = c(-100, 100),
                                      numeric.column = c("gain", "loss"),
                                      panel.fun = function(region, value, ...){
                                          value1 <- unlist(value[1])
                                          value2 <- unlist(value[2])
                                          pos <- unlist((region[2]+region[1])/2)
                                          # check if highlight exists
                                          if (!is.null(highlight)){
                                          # check if the plotting chr is highlight chr
                                              idx <- as.numeric(which(h_chr == circlize::CELL_META$sector.index))
                                          if (length(idx) != 0){
                                              bar.col.gain <- rep(col.gain,length(value1))
                                              bar.col.gain[which(unlist(region[1]) %in% h_start[idx])] <- 'red'
                                              bar.col.loss <- rep(col.loss,length(value1))
                                              bar.col.loss[which(unlist(region[1]) %in% h_start[idx])] <- 'red'
                                          }
                                        }
                                        circlize::circos.barplot(value=value1,pos=pos, border = bar.col.gain, col = bar.col.gain)
                                        circlize::circos.barplot(value=value2,pos=pos, border = bar.col.loss, col = bar.col.loss)
                                        y <- seq(-100,100,by=20)
                                        circlize::circos.segments(0,y,circlize::CELL_META$xlim[2], y,col=ifelse(y>0,col.gain,ifelse(y == 0,'#909090',col.loss)))
                                      })
        if (i <= 2 & length(id_name) <= 2){text(0,-0.1*(i-1),id,cex = cex_text)}
    }

}
