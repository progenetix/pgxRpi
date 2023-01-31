#' Plot CNV frequency data
#'
#' Thie function plots the frequency of deletions and amplifications
#'
#' @param data Frequency object returned by `pgxLoader` function or individual frequency matrix with 'pgxfreq' format.
#' In the returned object by `pgxLoader`, the frequency matrices in `data` slot must be stored as 'pgxfreq' format,
#' which can be specified by `output` parameter of `pgxLoader` function.
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
#' @export


pgxFreqplot <- function(data,chrom=NULL,layout=c(1,1),filters=NULL,circos = FALSE,highlight=NULL,assembly = 'hg38'){
    if (is.null(filters)){
        filters <- names(data$data)[1]
    }

    if (any(!(filters %in% names(data$data))) & any(is.na(match(filters,seq(length(data$data)))))){
        stop(paste("The filter is not contained in data:", filters[!(filters %in% names(data$data))]))
    }

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

  if (all(names(data) %in% c("data","meta"))){
    data_pgxfreq <- data$data[filters]
    id_name <-  names(data_pgxfreq)
  } else{
    attempt::stop_if(.x= length(colnames(data)) == 0, msg="\n The input is invalid \n")
    attempt::stop_if(.x= any(colnames(data)[c(5,6)] != c('gain_frequency','loss_frequency')), msg="\n The input is invalid \n")
    data_pgxfreq <- list()
    id_name <- unique(data[,1])
    for (i in id_name){
      data_pgxfreq[[i]] <- data[data[,1] %in% i,]
    }
  }

  # Identify highlight location
  if (!is.null(highlight)){
    h_chr <- paste0('chr',data_pgxfreq[[1]][highlight,2])
    h_start <- data_pgxfreq[[1]][highlight,3]
  }


  # set text size and color
  cex_text <- ifelse(length(data_pgxfreq) == 1, 0.8, 0.6)
  col.gain <-'#FFC633'
    col.loss <- '#33A0FF'
      bar.col.gain <- col.gain
      bar.col.loss <- col.loss

      # plot
      for (i in seq(length(id_name))){
        id <- id_name[i]
        data_cir <- data_pgxfreq[[id]][,c(2,3,4,5,6)]
        data_cir[,5] <- -data_cir[,5]
        data_cir[,1] <- paste0('chr',data_cir[,1])
        # Initialize
        if (i == 1){
          circlize::circos.clear()
          circlize::circos.par(start.degree=90, gap.degree=c(rep(1,23),10))
          circlize::circos.initializeWithIdeogram(species = assembly)
        }
        circlize::circos.genomicTrack(data_cir, ylim = c(-100, 100),numeric.column = c("gain_frequency", "loss_frequency"),
                                      panel.fun = function(region, value, ...) {
                                        value1 = unlist(value[1])
                                        value2 = unlist(value[2])
                                        pos=unlist((region[2]+region[1])/2)
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
                                        y = seq(-100,100,by=20)
                                        circlize::circos.segments(0,y,circlize::CELL_META$xlim[2], y,col=ifelse(y>0,col.gain,ifelse(y == 0,'#909090',col.loss)))
                                      })
        if (i <= 2 & length(id_name) <= 2){text(0,-0.1*(i-1),id,cex = cex_text)}
      }

}
