#' Extract, analyse and visualize "pgxseg" files
#'
#' This function extracts segments, CNV frequency, and metadata from local "pgxseg" files and
#' supports survival data visualization  
#'
#' @param file A string specifying the path and name of the "pgxseg" file where the data is to be read. 
#' @param group_id A string specifying which id is used for grouping in KM plot or CNV frequency calculation. Default is "group_id".
#' @param show_KM_plot A logical value determining whether to return the Kaplan-Meier plot based on metadata. Default is FALSE.
#' @param return_metadata A logical value determining whether to return metadata. Default is FALSE.
#' @param return_seg A logical value determining whether to return segment data. Default is FALSE.
#' @param return_frequency A logical value determining whether to return CNV frequency data. The frequency calculation is based on segments in segment data and specified group id in metadata. Default is FALSE.
#' @param assembly A string specifying which genome assembly version should be applied to CNV frequency calculation and plotting. Allowed options are "hg19", "hg38". Default is "hg38".
#' @param bin_size Size of genomic bins used in CNV frequency calculation to split the genome, in base pairs (bp). Default is 1,000,000.
#' @param overlap Numeric value defining the amount of overlap between bins and segments considered as bin-specific CNV, in base pairs (bp). Default is 1,000.
#' @param soft_expansion Fraction of `bin_size` to determine merge criteria.
#' During the generation of genomic bins, division starts at the centromere and expands towards the telomeres on both sides.
#' If the size of the last bin is smaller than `soft_expansion` * bin_size, it will be merged with the previous bin. Default is 0.1.
#' @param ... Other parameters relevant to KM plot. These include `pval`, `pval.coord`, `pval.method`, `conf.int`, `linetype`, and `palette` (see ggsurvplot from survminer)
#' @importFrom utils URLencode modifyList read.table write.table data
#' @importFrom methods show
#' @return Segments data, CNV frequency object, meta data or KM plots from local "pgxseg" files
#' @export
#' @examples
#' file_path <- system.file("extdata", "example.pgxseg",package = 'pgxRpi')
#' info <- pgxSegprocess(file=file_path,show_KM_plot = TRUE, return_seg = TRUE, return_metadata = TRUE)
pgxSegprocess <- function(file,group_id = 'group_id', show_KM_plot=FALSE,return_metadata=FALSE,return_seg=FALSE,return_frequency=FALSE,assembly='hg38',bin_size= 1000000,overlap=1000,soft_expansion = 0.1,...){
    if(!any(show_KM_plot, return_metadata, return_seg, return_frequency)){return()}
    full.data <- readLines(file)
    idx <- grep("#sample=>",full.data)
    meta.lines <- full.data[idx]
    meta <- lapply(meta.lines,function(x){
        info <- unlist(strsplit(x,split=';'))
        col_name <- sub('=.*','',info)[-1]
        col_data <- sub('.*=','',info)[-1]
        df <- as.data.frame(t(data.frame(col_data,row.names = col_name)))
        return(df)
    })
    meta <- dplyr::bind_rows(meta)
    rownames(meta) <- seq(dim(meta)[1])
    if (!group_id %in% colnames(meta)){
        stop('group_id: ',group_id,' is not in metadata')
    }
    if (show_KM_plot){
        # remove samples without survival data
        meta.sel.df <- meta[!is.na(meta$followup_time) & meta$followup_state_id != 'EFO:0030039',]
        if(dim(meta.sel.df)[1] == 0) stop("\n No available survival data \n")
          # transform ISOdate interal to days
        meta.sel.df$followup_time <- lubridate::time_length(meta.sel.df$followup_time,unit='day')
          # Map EFO code to numerical representation for plotting
        meta.sel.df$group_id <- meta.sel.df[,group_id]
        meta.sel.df$followup_state_id[meta.sel.df$followup_state_id == "EFO:0030041"] <- 0 # alive
        meta.sel.df$followup_state_id[meta.sel.df$followup_state_id == "EFO:0030049"] <- 1 # death
        meta.sel.df$followup_state_id <- as.numeric(meta.sel.df$followup_state_id )
        sfit <- survival::survfit(survival::Surv(followup_time, followup_state_id)~group_id , data=meta.sel.df)
    
        splot.op <- list(
          fit=sfit,
          data=meta.sel.df,
          palette=NULL,
          linetype=1,
          conf.int=FALSE,
          pval = FALSE,
          pval.method=FALSE,
          pval.coord = NULL,
          ggtheme=ggplot2::theme_light()
        )
        splot.op <- modifyList(splot.op,list(...))

        ggp <- do.call(survminer::ggsurvplot,splot.op)
        show(ggp$plot)
    }
   
    if (return_metadata | return_frequency){
        meta$followup_time <- ifelse(length(meta$followup_time) == 0,NA,
                                     round(lubridate::time_length(meta$followup_time,unit='day')))
        meta <- dplyr::mutate(meta,followup_state_label = NA,.after = followup_state_id)
        meta$followup_state_label[meta$followup_state_id == "EFO:0030041"] <- 'alive'
        meta$followup_state_label[meta$followup_state_id == "EFO:0030049"] <- 'dead'
        meta$followup_state_label[meta$followup_state_id == "EFO:0030039"] <- 'unknown'
        meta <- dplyr::rename(meta,'followup_time(days)' = followup_time)
    }
  
    if (return_seg  | return_frequency){
        seg <- read.table(file,sep='\t',header=TRUE)
        # sort chr and start per sample
        seg <- lapply(unique(seg[,1]),FUN = function(x){
        indseg <- seg[seg[,1] == x,]
        indseg <- indseg[order(as.numeric(indseg[,3])),]
        chr <- indseg[,2]
        chr[chr == 'X'] <- 23
        chr[chr == 'Y'] <- 24
        indseg <- indseg[order(as.numeric(chr)),]
        return(indseg)
        })
        seg <- do.call(rbind,seg)
        rownames(seg) <- seq_len(dim(seg)[1])
        if (return_frequency){
            bin.data <- extract.bin.feature(seg,genome=assembly,bin.size=bin_size,overlap=overlap,soft.expansion = soft_expansion)
            bin.dup.data <- bin.data[['dup']]
            bin.del.data <- bin.data[['del']]
  
            freq_data <- list()
            freq_meta <- list()
            for (id in unique(meta[,group_id])){
                plot.samples <- meta[,1][meta[,group_id] %in% id]
                # select samples from samples in pgxseg data
                plot.idx <- rownames(bin.dup.data) %in%  plot.samples
                nsamples <- sum(plot.idx)
                if (nsamples == 0){
                    freq.seg <- NULL
                }else {
                    if (nsamples == 1){
                        gain.freq <-  (bin.dup.data[plot.idx,]/nsamples)*100
                        loss.freq <-  (bin.del.data[plot.idx,]/nsamples)*100
                    } else{
                        gain.freq <-  (colSums(bin.dup.data[plot.idx,])/nsamples)*100
                        loss.freq <-  (colSums(bin.del.data[plot.idx,])/nsamples)*100
                    }
                    freq.seg.group <- data.frame(group_id=id)
                    freq.seg <- cbind(freq.seg.group,bin.data[['bin']][,c(2,3,4)])
                    freq.seg$gain_frequency <- gain.freq
                    freq.seg$loss_frequency <- loss.freq
                }
        
                freq_data[[id]] <- freq.seg
                group_label <- sub("_id", "_label",group_id)
                freq_meta[[id]] <-  data.frame(group_id = id, group_label=ifelse(group_label %in% colnames(meta),unique(meta[,group_label][meta[,group_id] == id]),NA),
                                              sample_count=nsamples)
            }
            
            freq_meta <- do.call(rbind,freq_meta)
            freq_data <- do.call(rbind,freq_data)
            
            frequency <- GenomicRanges::makeGRangesListFromDataFrame(freq_data,split.field = 'group_id',keep.extra.columns=TRUE)
            S4Vectors::mcols(frequency) <- freq_meta
        }
    }
  
    result <- list()
    if (return_seg){
        result[['seg']] <- seg
    } 
  
    if (return_metadata){
        result[['meta']] <- meta
    } 
  
    if (return_frequency){
        result[['frequency']] <- frequency
    }
  
    if (length(result) == 1){
        return(get(names(result)))
    }
  
    if (length(result) > 0){
        return(result)
    }
}
    
  



