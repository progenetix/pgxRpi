#' Extract, analyse and visualize "pgxseg" files
#'
#' This function extracts segments, CNV frequency, and metadata from "pgxseg" files and
#' supports survival data visualization  
#'
#' @param file A string specifying the path and name of the "pgxseg" file where the data is to be read. 
#' @param show_KM_plot A logical value determining whether to return the Kaplan-Meier plot based on metadata. Default is FALSE.
#' @param return_metadata A logical value determining whether to return metadata. Default is FALSE.
#' @param return_seg A logical value determining whether to return segment data. Default is FALSE.
#' @param return_frequency A logical value determining whether to return frequency data. The frequency calculation is based on segments in segment data and group id in metadata. Default is FALSE.
#' @param assembly A string specifying which genome assembly version should be applied to CNV frequency calculation and plotting. Allowed options are "hg19", "hg38". Default is "hg38".
#' @param ... Other parameters relevant to KM plot. These include `pval`, `pval.coord`, `pval.method`, `conf.int`, `linetype`, and `palette` (see ggsurvplot from survminer)
#' @return Segments data, CNV frequency object, meta data or KM plots from "pgxseg" files
#' @export

pgxSegprocess <- function(file,show_KM_plot=F,return_metadata=F,return_seg=F,return_frequency=F,assembly='hg38',...){
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
  if (show_KM_plot == T){
        # remove samples without survival data
    meta.sel.df <- meta[!is.na(meta$followup_time) & meta$followup_state_id != 'EFO:0030039',]
    attempt::stop_if(.x= dim(meta.sel.df)[1] == 0, msg="\n No available survival data \n")
      # transform ISOdate interal to days
    meta.sel.df$followup_time <- lubridate::time_length(meta.sel.df$followup_time,unit='day')
      # Map EFO code to numerical representation for plotting
    meta.sel.df$followup_state_id[meta.sel.df$followup_state_id == "EFO:0030041"] <- 0 # alive
    meta.sel.df$followup_state_id[meta.sel.df$followup_state_id == "EFO:0030049"] <- 1 # death
    meta.sel.df$followup_state_id <- as.numeric(meta.sel.df$followup_state_id )
      
    sfit <- survival::survfit(survival::Surv(followup_time, followup_state_id)~group_id , data=meta.sel.df)
    
    splot.op <- list(
      fit=sfit,
      data=meta.sel.df,
      palette=NULL,
      linetype=1,
      conf.int=F,
      pval = F,
      pval.method=F,
      pval.coord = NULL,
      ggtheme=ggplot2::theme_light()
    )
    splot.op <- modifyList(splot.op,list(...))

    ggp <- do.call(survminer::ggsurvplot,splot.op)
    print(ggp$plot)
    }
   
  if (return_metadata == T){
    meta$followup_time <- ifelse(length(meta$followup_time) == 0,NA,lubridate::time_length(meta$followup_time,unit='day'))
    meta$followup_state_id[meta$followup_state_id == "EFO:0030041"] <- 'alive'
    meta$followup_state_id[meta$followup_state_id == "EFO:0030049"] <- 'dead'
    meta$followup_state_id[meta$followup_state_id == "EFO:0030039"] <- 'unknown'
    meta <- dplyr::rename(meta,'followup_time(days)' = followup_time, 'followup_state'=followup_state_id)
  }
  
  if (return_seg == T | return_frequency == T){
    seg <- read.table(file,sep='\t',header=T)
    if (return_frequency == T){
      bin.data <- extract.bin.feature(seg,genome=assembly)
      bin.dup.data <- bin.data[['dup']]
      bin.del.data <- bin.data[['del']]
      frequency <- list()
      frequency[['data']] <- list()
      frequency[['meta']] <- c()
      total_sample_count <- 0
      for (id in unique(meta$group_id)){
        plot.samples <- meta[,1][meta$group_id %in% id]
        # select samples from samples in pgxseg data
        plot.idx <- rownames(bin.dup.data) %in%  plot.samples
        nsamples <- sum(plot.idx)
        if (nsamples == 0){
          freq.seg <- NULL
        } else {
          if (nsamples == 1){
            gain.freq <-  as.numeric(bin.dup.data[plot.idx,])/nsamples
            loss.freq <-  as.numeric(bin.del.data[plot.idx,])/nsamples
          } else{
          gain.freq <-  colSums(bin.dup.data[plot.idx,])/nsamples
          loss.freq <-  colSums(bin.del.data[plot.idx,])/nsamples
          }
          freq.seg.group <- data.frame(group_id=id)
          freq.seg <- cbind(freq.seg.group,bin.data[['bin']][,c(2,3,4)])
          freq.seg$gain_frequency <- gain.freq
          freq.seg$loss_frequency <- loss.freq
        }
        
        frequency$data[[id]] <- freq.seg
        frequency$meta <- rbind(frequency$meta, data.frame(group_id = id, 
                                                           group_label=unique(meta$group_label[meta$group_id == id]),
                                                           sample_count=nsamples))
        total_sample_count <-  total_sample_count+nsamples
      }
      frequency$meta <- rbind(frequency$meta, data.frame(group_id = 'total', 
                                                         group_label='',
                                                         sample_count=total_sample_count))
      frequency$data[['total']] <- Reduce(rbind,frequency$data)

      }
    }
  
  result <- list()
  if (return_seg == T){
    result[['seg']] <- seg
  } 
  
  if (return_metadata == T){
    result[['meta']] <- meta
  } 
  
  if (return_frequency == T){
    result[['frequency']] <- frequency
  }
  
  if (length(result) == 1){
    return(get(names(result)))
  }
  
  if (length(result) > 0){
    return(result)
  }
}
    
  



