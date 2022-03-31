#' @import attempt

check_internet <- function(){
    stop_if_not(.x = curl::has_internet(), msg = "Please check your internet connexion")
}

transform_id <- function(id){
    if (length(id) > 1) {
        filter <- id[1]
        for (tag in id[2:length(id)]) {
            filter  <-paste(filter, tag, sep=',')
        } }else{
            filter <- id
        }
    return(filter)
}

read_sample <- function(url){
    header <- read.table(url, nrow=1,stringsAsFactors = FALSE, sep = "\t",fill=TRUE,quote="",header=F)
    data <- read.table(url, skip=1,stringsAsFactors = FALSE, sep = "\t",fill=TRUE,quote="",header=F)
    header <- gsub('___','_',header)
    colnames(data) <- header
    data <- data[grepl('pgxbs',data$biosample_id),]
    suppressWarnings(data$death <- as.numeric(data$death))
    suppressWarnings(data$followup_months <- as.numeric(data$followup_months))
    return(data.frame(data))
}

read_variant_pgxseg <- function(url){
    result <- read.table(url, header = T, sep="\t")
    colnames(result)[2] <- 'chromosome'
    col <- c("start","end","log2")
    suppressWarnings(result[,col] <- sapply(result[,col], as.numeric))
    result <- result[order(result$start),]
    chr <- result$chromosome
    chr[which(chr == 'X')] <- 23
    chr[which(chr == 'Y')] <- 24
    chr <- as.integer(chr)
    result <- result[order(chr),]
    result <- result[order(result$biosample_id),]
    return(result)
}

read_variant_pgxseg_meta <- function(url){
    meta <- readLines(url)
    idx <- length(grep("#",meta))
    return(meta[c(1:idx)])
}

pgidCheck <- function(id){
    id_prefix <- id
    idx <- grep("pgx",id_prefix)
    if (length(idx) != 0){
      id_prefix[idx] <- gsub(".*:","",id_prefix[idx])
      id_prefix[idx] <- gsub("-.*","",id_prefix[idx])
    }
    id_prefix <- unique(sub(":.*", "", id_prefix))
    total_url <- url(description=paste0("https://progenetix.org/services/collations?filters=",id_prefix[1],"&output=text"),
                     open='r')
    info <- read.table(total_url, header=F, sep="\t", na="NA",fill=TRUE,quote = "")
    close(total_url)
    if (length(id_prefix) > 1){
        for (i in c(2:length(id_prefix))){
            total_url <- url(description=paste0("https://progenetix.org/services/collations?filters=",id_prefix[i],"&output=text"),
                             open='r')
            temp <- read.table(total_url, header=F, na="NA",sep = '\t',fill=TRUE,quote = "")
            close(total_url)
            info <- rbind(info, temp)
        }
    }

    return(id %in% unlist(info[1]))
}

pgxFreqLoader <- function(output, codematches, filters) {
    # check filters again
   stop_if(.x=length(filters) < 1, msg="\n At least one valid filter has to be provided \n")
    # check output format
    stop_if(.x=!(output %in% c('pgxseg','pgxmatrix')),msg="\n Output is invalid. Only support 'pgxseg' or 'pgxmatrix' \n")
    # check if filters exists
    idcheck <- pgidCheck(filters)

    if (!all(idcheck)){
        cat(" No results for id", filters[!idcheck], "in progenetix database.\n","\n Only query id:", filters[idcheck],'\n')
        filters <- filters[idcheck]
    }
    # start query
    pg.url <- paste0("http://www.progenetix.org/services/intervalFrequencies/?output=",output)


    if (length(filters)>1){
        filter <- transform_id(filters)
        pg.url  <- paste(pg.url,'&filters=',filter,sep="")
        } else{
        pg.url  <- paste(pg.url, '&', 'id', '=', filters,sep="")
    }

    if (codematches){
        pg.url  <- paste0(pg.url, '&method=codematches')
    }

    meta <- readLines(pg.url)[c(1:(length(filters)+1))]
    meta_lst <- unlist(strsplit(meta,split = ';'))
    label <- meta_lst[grep("label",meta_lst)]
    label<- gsub('.*=','',label)
    count <- meta_lst[grep("sample_count",meta_lst)]
    count <- as.numeric(gsub('.*=','',count))
    meta <- data.frame(code = c(filters,'total'), label=c(label,''), sample_count=c(count,sum(count)))
    pg.data  <- read.table(pg.url, header=T, sep="\t", na="NA")
    colnames(pg.data)[1] <- 'filters'
    data_lst <- list()
    for (i in filters){
        data_lst[[i]] <- pg.data[pg.data$filters == i,]
    }
    data_lst[['total']] <- pg.data
    result <- list(meta = meta, data = data_lst)

    return(result)
}

pgxSampleLoader <- function(biosample_id,filters,codematches,skip,limit){
    if (!(is.null(filters))){
        for (i in c(1:length(filters))) {
            if_next <- FALSE
            url <- paste0("http://progenetix.org/cgi/bycon/beaconServer/biosamples.py?filters=",
                              filters[i],"&output=table")
            url  <- ifelse(is.null(limit), url, paste0(url,"&limit=",limit))
            url  <- ifelse(is.null(skip), url, paste0(url,"&skip=",skip))
            if (!(exists('res_1'))){
                try_catch(res_1 <- read_sample(url),.e= function(e){
                    cat(paste("No samples with the filter", filters[i],"\n"))
                    if_next <<- TRUE})
                if (if_next){ next }
            }else {
                try_catch(temp <- read_sample(url),.e= function(e){
                    cat(paste("No samples with the filter", filters[i],"\n"))
                    if_next <<- TRUE
                    })
                if (if_next){ next }
                res_1 <- rbind(res_1, temp)}
            }}
    if (!(is.null(biosample_id))){
        len <- length(biosample_id)
        count <- ceiling(len/50)
        for (i in c(1:count)){
            if_next <- FALSE
            j = i*50
            if (j > len){
                j<-len
            }
            filter <- transform_id(biosample_id[c(((i-1)*50+1):j)])
            url <- paste0("http://progenetix.org/cgi/bycon/beaconServer/biosamples.py?biosampleIds=",
                              filter,"&output=table")
            if (!(exists('res_2'))){
                try_catch(res_2 <- read_sample(url),.e= function(e){
                    if_next <<- TRUE})
                if (if_next){ next }
            }else {
                try_catch(temp <- read_sample(url),.e= function(e){
                    if_next <<- TRUE
                    })
                if (if_next){ next }
                res_2 <- rbind(res_2, temp)}
            }}
    stop_if(.x= !(exists('res_1') | exists('res_2')) , msg='No samples retrieved')

    if (exists('res_1') & exists('res_2')){
        res <- rbind(res_1, res_2)
    } else if (exists('res_1')){
        res <- res_1
    } else{
        res <- res_2
    }


    if (codematches){
        idx <- res$histological_diagnosis__id %in% filters | res$sampled_tissue__id %in% filters | res$icdo_morphology__id %in% filters |
        res$icdo_topography__id %in% filters | res$external_references__id_PMID %in% filters | res$biosample_id %in% biosample_id |
        res$external_references__id_geo.GSM %in% filters | res$external_references__id_geo.GSE %in% filters | 
        res$external_references__id_geo.GPL %in% filters | res$external_references__id_cellosaurus %in% filters | 
        res$external_references__id_arrayexpress %in% filters
        
        res <- res[idx,]
        if (dim(res)[1] == 0){
            cat("\n WARNING: the option `codematches=TRUE` filters out all samples \n")
        }}
    
    colnames(res)[c(1,3)] <- c('biosample_id','callset_id')
    return(res)
}

pgxVariantLoader <- function(biosample_id, output, save_file,filename){
    
    if (save_file){
        stop_if(.x= !(output %in% c("pgxseg","seg")), msg = "The parameter 'output' is invalid (available: 'seg' or 'pgxseg')")
    }

    len <- length(biosample_id)
    count <- ceiling(len/50)
    for (i in c(1:count)){
        if_next <- FALSE
        j = i*50
        if (j > len){
            j<-len
        }
        filter <- transform_id(biosample_id[c(((i-1)*50+1):j)])
        url <- paste0("http://progenetix.org/cgi/bycon/beaconServer/variants.py?biosampleIds=",
                      filter)
        if (!(is.null(output))){
            if (output == 'seg' | output == 'pgxseg'){
                url <- paste0(url, "&output=pgxseg")
                try_catch(temp <- read_variant_pgxseg(url), .e = function(e){if_next <<- TRUE}, 
                          .w = function(w){if_next <<- TRUE})
                if (if_next){ next }
            } else{
                stop("output is invalid (NULL, 'seg' or 'pgxseg')")
            }
        }else{
            try_catch(temp <- rjson::fromJSON(file = url), .e = function(e){if_next <<- TRUE}, .w = function(w){if_next <<- TRUE})
            if (if_next){  next }
            temp <- lapply(temp$response$resultSets[[1]]$results,unlist)
            temp <- lapply(temp,function(x){
                var_meta <- x[c('variation.location.sequenceId','variation.variantInternalId','variation.relativeCopyClass','variation.updated')]
                id_ind <-  which(names(x) =='caseLevelData.id')
                temp_row <- c()
                for (i in id_ind){
                    temp_row <- rbind(temp_row,c(x[i],x[i-1],x[i-2]))}
                var_meta <- matrix(rep(var_meta,dim(temp_row)[1]),nrow=dim(temp_row)[1],
                                   byrow = T)
                temp_row <- as.data.frame(cbind(temp_row,var_meta))
                colnames(temp_row) <- c("variant_id","biosample_id","analysis_id",
                                        "reference_genome","variant","variant_status","updated_time")
                return(temp_row)
            })
            temp <- Reduce(rbind,temp)
            # change interval order 
            location <- strsplit(temp$variant,split = ':')
            chr <-  sapply(location, function(x){x[1]})
            chr[chr == 'X'] <- 23
            chr[chr == 'Y'] <- 24
            chr <- as.numeric(chr)
            start <- sapply(location, function(x){x[2]})
            end <- as.numeric(gsub('.*-','',start))
            start <- as.numeric(gsub('-.*','',start))
            temp$chr <- chr
            temp$start <- start
            temp$end <- end
            temp <- temp[order(temp$start),]
            temp <- temp[order(temp$chr),]
            temp <- temp[order(temp$biosample_id),]
            temp$chr[temp$chr == 23] <- 'X'
            temp$chr[temp$chr == 24] <- 'Y'
            temp <- temp[,c(1,2,3,5,8,9,10,4,6,7)]
        }

        if (!(exists("res"))){
            res <- temp
            if (save_file){
                if (output == 'pgxseg'){
                    meta <- read_variant_pgxseg_meta(url)
                }} 
            }else{
                res <- rbind(res, temp)
                if (save_file){
                    if (output == 'pgxseg'){
                    meta <- c(meta, meta[-c(1,2)])
                }}
            }
    }
    if (!(exists("res"))){
        cat("No variants retrieved")
        return(NULL)}
    
    if (!(is.null(output))){
        if (output == 'seg'){
            res <- res[,c(1,2,3,4,6,5)]
        }}
    
    # quality check
    res <- res[res$biosample_id %in% biosample_id,]
    rownames(res) <- seq(nrow(res))
    if (save_file){
        if (output=='pgxseg'){
            if (is.null(filename)){
                filename <- "variants.pgxseg"
            }

            write.table(meta, file=filename,row.names = FALSE,col.names = FALSE, quote = FALSE)
            suppressWarnings(write.table(res, append=TRUE, sep='\t',file=filename,row.names = FALSE,col.names = TRUE, quote = FALSE))
            } else if(output == 'seg'){
            if (is.null(filename)){
                filename <- "variants.seg"
            }
            
            write.table(res, file=filename, sep='\t',row.names = FALSE,col.names = TRUE, quote = FALSE)
        } else{
            stop("Output is null. Please specify output format")
        }
        return(cat("The file is saved"))
    }

    return(res)
}


pgxcallsetLoader <- function(filters,limit,skip,codematches){
  # check filters 
  stop_if(.x=length(filters) < 1, msg="\n at least one valid filter has to be provided \n")
  
  idcheck <- pgidCheck(filters)
  if (!all(idcheck)){
    cat(" No results for id", filters[!idcheck], "in progenetix database.\n","\n Only query id:", filters[idcheck],'\n')
    filters <- filters[idcheck]
  }
  # start query
  pg.url <- "http://www.progenetix.org/beacon/callsets/?output=pgxmatrix"
  
  filter <- transform_id(filters)
  pg.url  <- paste(pg.url,'&filters=',filter,sep="")
  
  pg.url  <- ifelse(is.null(limit), pg.url, paste0(pg.url,"&limit=",limit))
  pg.url  <- ifelse(is.null(skip), pg.url, paste0(pg.url,"&skip=",skip))
  
  meta <- readLines(pg.url,n=7)
  meta <-  gsub("#meta=>\"","",meta[7])
  meta <-  gsub("\"","",meta)
  if (length(grep('WARNING',meta)) == 1){
    cat(paste("\n", meta, "\n"))
  }
  pg.data  <- read.table(pg.url, header=T, sep="\t", na="NA")
  
  if (codematches){
    pg.data <- pg.data[pg.data$group_id %in% filters,]
    if (dim(pg.data)[1] == 0){
      cat("\n Warning: the option `codematches=TRUE` filters out all samples \n")
    }}
  
  return(pg.data)
}


pgxCovLoader <- function(filters,codematches,skip,limit){
  stop_if(.x=length(filters) < 1, msg="\n One valid filter has to be provided \n")
  stop_if(.x=length(filters) > 1, msg="\n This query only supports one filter \n")
  if (codematches){
    biosample_id <- pgxSampleLoader(biosample_id = NULL,filters = filters,codematches = T,skip=NULL,limit=NULL)
    biosample_id <- biosample_id$biosample_id
    count <- pgxCount(filters)[,3]
    if (count > 2000){ 
      count <- floor(count/2000)
      for (i in seq(count)){
        temp_biosample_id <- pgxSampleLoader(biosample_id = NULL,filters = filters,codematches = T,skip=i,limit=NULL)
        temp_biosample_id <-  temp_biosample_id$biosample_id
        biosample_id <- c(biosample_id, temp_biosample_id)
      }
    }
    
    if (length(biosample_id) == 0){
      return()
    }}
  
  url <- paste0('https://progenetix.org/beacon/callsets/?filters=',filters)
  url  <- ifelse(is.null(limit), url, paste0(url,"&limit=",limit))
  url  <- ifelse(is.null(skip), url, paste0(url,"&skip=",skip))
  data <- rjson::fromJSON(file = url)
  sample <- unlist(lapply(data$response$resultSets[[1]]$results, function(x){
    return(x$biosampleId)
  }))
  
  if (codematches){
    sel_sample_idx <- sample %in% biosample_id 
  } else{
    sel_sample_idx <- seq(length(sample))
  }

  
  data_2 <- lapply(data$response$resultSets[[1]]$results[sel_sample_idx], function(x){
    stats <- lapply(x$cnvChroStats, function(y){
      as.data.frame(y)
    })
    cnvstats <- Reduce(rbind,stats)
    rownames(cnvstats) <- names(stats)
    res <- list()
    res[[1]] <- cnvstats
    res[[2]] <- as.data.frame(x$cnvStats)
    return(res)
  })
  
  
  total_frac <- lapply(data_2, function(x){
    x[[2]]
  })
  total_frac <- Reduce(rbind, total_frac)
  rownames(total_frac) <- make.unique(sample[sel_sample_idx])
  total_frac <- total_frac[,c(2,6,4)]
  
  data_2 <- lapply(data_2,function(x){
    x[[1]]
  })
  
  chrom_arm_list <- c()
  for ( i in c(seq(22),'x','y')){
    chrom_arm_list <- c(chrom_arm_list, paste0(i,'p'), paste0(i,'q'))
  }
  chrom_arm_idx <- match(chrom_arm_list, rownames(data_2[[1]]))
  
  chrom_list <- c(seq(22),'X','Y')
  chrom_idx <- match(chrom_list, rownames(data_2[[1]]))
  
  data_3 <- lapply(data_2, function(x){
    val <- c(x$dupfraction[chrom_arm_idx],x$delfraction[chrom_arm_idx],x$dupfraction[chrom_idx],x$delfraction[chrom_idx])
    names <- c(paste0(chrom_arm_list,'.dup'),paste0(chrom_arm_list,'.del'),paste0(chrom_list,'.dup'),paste0(chrom_list,'.del'))
    return(  t(data.frame(val,row.names =  names)))
  }) 
  
  data_3 <- Reduce(rbind,data_3)
  rownames(data_3) <- make.unique(sample[sel_sample_idx])
  
  arm_frac <- data_3[,c(1:96)]
  chrom_frac <- data_3[,c(97:144)]

  result <- list()
  result$chrom_arm_coverage <- arm_frac
  result$whole_chrom_coverage <- chrom_frac
  result$whole_genome_coverage <- total_frac
  return(result)
}
