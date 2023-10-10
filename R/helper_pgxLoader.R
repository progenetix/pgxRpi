transform_id <- function(id){
    filter <- paste(id,collapse = ",")
    return(filter)
}

read_variant_pgxseg <- function(url){
    result <- read.table(url, header = TRUE, sep="\t")
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
    return(meta[seq_len(idx)])
}

disease_code_check <- function(id, remain.id, url){
  total.id <- c()
  for (code in c('NCIT','icdom','icdot','UBERON')){
    idx <- grep(code,id)
    if (length(idx) > 0){
      encoded_url <- URLencode(paste0(url,code))
      res <- rjson::fromJSON(file = encoded_url)
      total.id <- c(total.id,unlist(lapply(res$response$results, function(x){x$id})))
      remain.id <- remain.id [!remain.id %in% id[idx]]
    }
  }
  return(list(total=total.id,remain=remain.id))
}

pgidCheck <- function(id,domain){
    # this query doesn't work for individual GSM id and age filters
    geogsm.idx <- grep('geo:GSM',id)
    age.idx <- grep('age:',id)
    pass.idx <- c(geogsm.idx,age.idx)
    if (length(pass.idx > 0)){
      remain.id <- id[-pass.idx]
    } else{
      remain.id <- id
    }
    
    url <- paste0(domain,"/services/collations?filters=")
    
    # check if disease_code associated filters exist in Progenetix
    check.res <- disease_code_check(id, remain.id, url)
    remain.id <- check.res$remain
    total.id <- check.res$total
    # check if non disease_code associated filters exist in Progenetix
    if (length(remain.id) > 0){
      encoded_url <- URLencode(paste0(url,transform_id(remain.id)))
      res <- rjson::fromJSON(file = encoded_url)
      total.id <- c(total.id, unlist(lapply(res$response$results, function(x){x$id})))
    }
    
    return(id %in% c(total.id,id[pass.idx]))
}

pgxFreqLoader <- function(output, codematches, filters, domain) {
    # check if filters exists
    idcheck <- pgidCheck(filters,domain)
    if (!all(idcheck)){
        warning("\n No results for filters ", filters[!idcheck], " in progenetix database.\n")
        filters <- filters[idcheck]
    }
    # start query
    url <- paste0(domain,"/services/intervalFrequencies/?output=",output)
  
    filter <- transform_id(filters)
    url  <- paste0(url,'&filters=',filter)
  
    if (codematches){
        url  <- paste0(url, '&method=codematches')
    }

    encoded_url <- URLencode(url)
    # access metadata from JSON output
    meta <- readLines(encoded_url)[seq_len(length(filters))+1]
    meta_lst <- unlist(strsplit(meta,split = ';'))
    label <- meta_lst[grep("label",meta_lst)]
    label<- gsub('.*=','',label)
    count <- meta_lst[grep("sample_count",meta_lst)]
    count <- as.numeric(gsub('.*=','',count))
    meta <- data.frame(filter = filters, label, sample_count=count)
    # access data 
    pg.data  <- read.table(encoded_url, header=TRUE, sep="\t")
    colnames(pg.data)[1] <- 'filters'
    if (output == 'pgxfreq'){
      colnames(pg.data)[2] <- 'chr'
      data_lst <- GenomicRanges::makeGRangesListFromDataFrame(pg.data,split.field = 'filters',keep.extra.columns=TRUE)
      S4Vectors::mcols(data_lst) <- meta
    } else if (output == "pgxmatrix"){
      frequency <- as.data.frame(t(pg.data[,-1]))
      rowRanges <- rownames(frequency)
      rowRanges <- lapply(rowRanges,function(x){
        range_str <- unlist(strsplit(x,split = '\\.'))
        chr <- range_str[1]
        chr <- gsub("X","",chr)
        if (chr == "") chr <- 'X'
        return(data.frame(chr=chr,start=range_str[2],end = range_str[3],type=range_str[4]))
      })
      rowRanges <- do.call(rbind, rowRanges)
      rowRanges <- GenomicRanges::GRanges(rowRanges)
      names(rowRanges) <- seq_len(dim(frequency)[1])
      
      colnames(frequency) <- pg.data[,1]
      rownames(frequency) <- seq_len(dim(frequency)[1])
      
      colData <- meta
      data_lst <- SummarizedExperiment::SummarizedExperiment(assays=list(frequency=frequency),
                           rowRanges=rowRanges, colData=meta)
    }
    
    return(data_lst)
}

pgxmetaLoader <- function(type, biosample_id, individual_id, filters, codematches, skip, limit, filterLogic, domain){
    if (type == "biosample"){
        query_key <- "biosamples"
    }

    if (type == "individual"){
        query_key <- "individuals"
        # avoid to query upper level info based on lower level info
        biosample_id <- NULL
    }

    if (!(is.null(filters))){
        # check if filters exists
        idcheck <- pgidCheck(filters,domain)
        if (!all(idcheck)){
            warning("\n No results for filters ", filters[!idcheck], " in progenetix database.\n")
            filters <- filters[idcheck]
        }
        
        if (filterLogic == "AND"){
            filters <- transform_id(filters)
            url <- paste0(domain,"/beacon/",query_key,"/?filters=",filters,"&output=datatable")
            url  <- ifelse(is.null(limit), url, paste0(url,"&limit=",limit))
            url  <- ifelse(is.null(skip), url, paste0(url,"&skip=",skip))
            encoded_url <- URLencode(url)
            attempt::try_catch(
              res_1 <- read.table(encoded_url,stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header=TRUE),
              .e= function(e){
                warning("\n Query fails for the filters ", filters,"\n")
              })
         
        } else if (filterLogic == "OR"){
          # pick age filter
          age.idx <- grep('age:',filters)
          trans.filters <- filters
          if (length(age.idx)>0){
              age.filter <- filters[age.idx]
              trans.filters <- filters[-age.idx]
              trans.filters <- sapply(trans.filters, FUN=function(x){transform_id(c(x,age.filter))})
          } 

          res_1 <- list()
          for (i in seq_len(length(trans.filters))){
              url <- paste0(domain,"/beacon/",query_key,"/?filters=",trans.filters[i],"&output=datatable")
              url  <- ifelse(is.null(limit), url, paste0(url,"&limit=",limit))
              url  <- ifelse(is.null(skip), url, paste0(url,"&skip=",skip))
              encoded_url <- URLencode(url)
      
              attempt::try_catch(
                res_1[[i]] <- read.table(encoded_url,stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header=TRUE),
                .e= function(e){
                  warning("\n Query fails for the filter ", trans.filters[i],"\n")
                })
          }

          res_1 <- do.call(rbind, res_1)
          
        } else{
          stop("filterLogic is invalid ('AND' or 'OR')")
        }
    }


    if (!(is.null(biosample_id))){
        len <- length(biosample_id)
        count <- ceiling(len/50)
        res_2 <- list()
        for (i in seq_len(count)){
            j <- i*50
            if (j > len){
                j<-len
            }
            filter <- transform_id(biosample_id[c(((i-1)*50+1):j)])
            url <- paste0(domain,"/beacon/",query_key,"/?biosampleIds=",filter,"&output=datatable")
            encoded_url <- URLencode(url)
            attempt::try_catch(res_2[[i]] <- read.table(encoded_url,stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header=TRUE),.e= function(e){
                warning("\n Query fails for biosample_id ", filter, "\n")
            })           
        }
        res_2 <- do.call(rbind, res_2)
    }

    if (!(is.null(individual_id))){
        len <- length(individual_id)
        count <- ceiling(len/50)
        res_3 <- list()
        for (i in seq_len(count)){
            j <- i*50
            if (j > len){
                j<-len
            }
            filter <- transform_id(individual_id[c(((i-1)*50+1):j)])
            url <- paste0(domain,"/beacon/",query_key,"/?individualIds=",filter,"&output=datatable")
            encoded_url <- URLencode(url)
            attempt::try_catch(res_3[[i]] <- read.table(encoded_url,stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header=TRUE),.e= function(e){
                warning("\n Query fails for individual_id ", filter, "\n")
            })
        }
        res_3 <- do.call(rbind, res_3)
    }

    attempt::stop_if(.x= !(exists('res_1') | exists('res_2')| exists('res_3')) , msg='No data retrieved')

    res <- c()
    if (exists('res_1')){
        res <- res_1 
    } 
    
    if (exists('res_2')){
        res <- plyr::rbind.fill(res, res_2)
    }

    if (exists('res_3')){
        res <- plyr::rbind.fill(res, res_3)
    }

    rownames(res) <- seq(dim(res)[1])
    if (codematches){
        if (type == "biosample"){
            idx <- res$biosample_id %in% biosample_id | res$individual_id %in% individual_id | 
            res$histological_diagnosis_id %in% filters | res$sampled_tissue_id %in% filters | 
            res$icdo_morphology_id %in% filters | res$icdo_topography_id %in% filters 
        } else if (type == "individual"){
            idx <- res$histological_diagnosis_id %in% filters | res$individual_id %in% individual_id
        }     
    
        res <- res[idx,]
        if (dim(res)[1] == 0){
            warning("\n the option `codematches=TRUE` filters out all samples \n")
        }
    }
    
    res <- res[!duplicated(res),]
    return(res)
}

pgxVariantLoader <- function(biosample_id, output, save_file, filename, domain){ 
    if (save_file & is.null(output)){
      stop("The parameter 'output' is invalid when 'save_file=TRUE' (available: \"seg\" or \"pgxseg\")") 
    } 
             
    len <- length(biosample_id)
    count <- ceiling(len/50)
    res <- list()
    meta <- c()
    for (i in seq_len(count)){
        j <- i*50
        if (j > len){
            j<-len
        }
        filter <- transform_id(biosample_id[c(((i-1)*50+1):j)])
        url <- paste0(domain,"/beacon/variants/?biosampleIds=",
                      filter,"&limit=0")
        if (!(is.null(output))){
            url <- paste0(url, "&output=pgxseg")
            encoded_url <- URLencode(url)
            attempt::try_catch(res[[i]] <- read_variant_pgxseg(encoded_url), .e = function(e){
                warning("\n Query fails for biosample_id ", filter, "\n")
            })
        }else{
            encoded_url <- URLencode(url)
            attempt::try_catch(res[[i]] <- rjson::fromJSON(file = encoded_url), .e = function(e){
                 warning("\n Query fails for biosample_id ", filter, "\n")
            })
            
            if (length(res) < i){
                res[[i]] <- NA
                next
            } 
            
            res[[i]] <- lapply(res[[i]]$response$resultSets[[1]]$results,unlist)
                          
            if (length(res[[i]]) == 0){
                warning("\n Query fails for biosample_id ", filter, "\n")
                res[[i]] <- NA
                next
            }
            
            temp <- res[[i]]              
            temp <- lapply(temp,function(x){
                var_meta <- x[c('variation.location.sequence_id','variantInternalId','variation.relative_copy_class','variation.updated')]
                id_ind <-  which(names(x) =='caseLevelData.id')
                temp_row <- c()
                for (i in id_ind){
                    temp_row <- rbind(temp_row,c(x[i],x[i-1],x[i-2]))}
                var_meta <- matrix(rep(var_meta,dim(temp_row)[1]),nrow=dim(temp_row)[1],
                                   byrow = TRUE)
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
            res[[i]] <- temp
        }

        if (save_file){
            if (output == 'pgxseg'){
                if (length(meta) == 0){
                     meta <- read_variant_pgxseg_meta(encoded_url)
                }else{
                     meta <- c(meta, meta[-c(1,2)])
                }
            }
        }

    }

    res[is.na(res)] <- NULL

    if (length(res) == 0){
            warning("\n No variants retrieved \n")
            return()
    }

    res <- do.call(rbind, res)




    if (!(is.null(output))){
        if (output == 'seg'){
            res <- res[,c(1,2,3,4,6,5)]
        }
    }
    
    # quality check
    res <- res[res$biosample_id %in% biosample_id,]
    rownames(res) <- seq(nrow(res))

    if (save_file){
        if (output=='pgxseg'){
            if (is.null(filename)) filename <- "variants.pgxseg"            
            write.table(meta, file=filename,row.names = FALSE,col.names = FALSE, quote = FALSE)
            suppressWarnings(write.table(res, append=TRUE, sep='\t',file=filename,row.names = FALSE,col.names = TRUE, quote = FALSE))
        } else if(output == 'seg'){
            if (is.null(filename))  filename <- "variants.seg"  
            write.table(res, file=filename, sep='\t',row.names = FALSE,col.names = TRUE, quote = FALSE)
        } 
        message("\n The file is saved \n")
        return()
    }

    return(res)
}


pgxcallsetLoader <- function(filters,limit,skip,codematches,domain){
    # start query
    url  <- paste0(domain,"/beacon/analyses/?output=pgxmatrix",'&filters=',filters)
    url  <- ifelse(is.null(limit), url, paste0(url,"&limit=",limit))
    url  <- ifelse(is.null(skip), url, paste0(url,"&skip=",skip))
    
    encoded_url <- URLencode(url)
    meta <- readLines(encoded_url,n=7)
    meta <-  gsub("#meta=>\"","",meta[7])
    meta <-  gsub("\"","",meta)
    if (length(grep('WARNING',meta)) == 1 & limit != 0){
        warning("\n", meta, "\n")
    }

    pg.data  <- read.table(encoded_url, header=TRUE, sep="\t")
    # remove automatic prefix X 
    colnames(pg.data) <- gsub("X","",colnames(pg.data))
    # add chr prefix to avoid colnames with numeric start
    colnames(pg.data) <- paste0("chr",colnames(pg.data) )
    # recover X chr    
    colnames(pg.data) <- gsub("chr\\.","chrX\\.",colnames(pg.data))

    if (codematches){
        pg.data <- pg.data[pg.data$group_id %in% filters,]
        if (dim(pg.data)[1] == 0){
            warning("\n the option `codematches=TRUE` filters out all samples \n")
    }}
  
    return(pg.data)
}


pgxCovLoader <- function(filters,codematches,skip,limit,domain){
    if (codematches){
        biosample_id <- pgxmetaLoader(type = 'biosample', biosample_id = NULL,individual_id = NULL,filters = filters,codematches = TRUE,skip=NULL,limit=0,filterLogic="AND",domain=domain)
        biosample_id <- biosample_id$biosample_id    
        if (length(biosample_id) == 0){
            warning("\n the option `codematches=TRUE` filters out all samples \n")
            return()
        }
    }
  
    url <- paste0(domain,'/beacon/analyses/?output=cnvstats&filters=',filters)
    url  <- ifelse(is.null(limit), url, paste0(url,"&limit=",limit))
    url  <- ifelse(is.null(skip), url, paste0(url,"&skip=",skip))
    encoded_url <- URLencode(url)
    data <- rjson::fromJSON(file = encoded_url)
    sample <- unlist(lapply(data$response$resultSets[[1]]$results, function(x){
        return(x$biosampleId)
    }))
  
    if (codematches){
        sel_sample_idx <- sample %in% biosample_id 
        if (sum(sel_sample_idx) == 0 & length(sample) > 0){
            warning("\n the option `codematches=TRUE` filters out all samples \n")
            return()
        }
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
  
  
    total_frac <- lapply(data_2, function(x){x[[2]]})
    # some results are SNV instead of CNV
    rm_idx <- which(sapply(total_frac,length) == 0)
    if (length(rm_idx) > 0){
        sel.sample <- sample[sel_sample_idx][-rm_idx]
        data_2 <- data_2[-rm_idx]
    } else{
    sel.sample <- sample[sel_sample_idx]
    }
    total_frac <- Reduce(rbind, total_frac)
    rownames(total_frac) <- make.unique(sel.sample)
    total_frac <- total_frac[,c(2,6,4)]
  
    data_2 <- lapply(data_2,function(x){x[[1]]})
  
    chrom_arm_list <- c()
    for ( i in c(seq(22),'x','y')){
        chrom_arm_list <- c(chrom_arm_list, paste0(i,'p'), paste0(i,'q'))
    }
    chrom_arm_idx <- match(chrom_arm_list, rownames(data_2[[1]]))
  
    chrom_list <- c(seq(22),'X','Y')
    chrom_idx <- match(chrom_list, rownames(data_2[[1]]))
  
    data_3 <- lapply(data_2, function(x){
        val <- c(x$dupfraction[chrom_arm_idx],
                 x$delfraction[chrom_arm_idx],
                 x$dupfraction[chrom_idx],
                 x$delfraction[chrom_idx])
        names <- c(paste0('chr',chrom_arm_list,'.dup'),
                   paste0('chr',chrom_arm_list,'.del'),
                   paste0('chr',chrom_list,'.dup'),
                   paste0('chr',chrom_list,'.del'))
        return(  t(data.frame(val,row.names =  names)))
    }) 
  
    data_3 <- Reduce(rbind,data_3)
    rownames(data_3) <- make.unique(sel.sample)
  
    arm_frac <- data_3[,seq_len(96)]
    chrom_frac <- data_3[,c(97:144)]

    result <- list()
    result$chrom_arm_coverage <- arm_frac
    result$whole_chrom_coverage <- chrom_frac
    result$whole_genome_coverage <- total_frac
    return(result)
}



