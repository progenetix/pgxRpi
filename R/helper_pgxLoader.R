checkUnusedParameters <- function(param, param_name, rightparam_name) {
  if (!is.null(param)) {
    rightparam <- paste(rightparam_name,collapse = " and ")
    warning("\n The parameter ", param_name, " is not used in this query. Only ",rightparam," are accepted. \n")
  }
}

checkMissingParameters <- function(param, param_name){
  if (length(param) < 1) {
    param_name <- paste(param_name,collapse = " or ")
    stop("\n The parameter ",param_name," is missing. \n")
  }
}

transform_id <- function(id){
    filter <- paste(id,collapse = ",")
    return(filter)
}

read_variant_pgxseg <- function(biosample_id, domain, dataset){
    url <- paste0(domain,"/services/pgxsegvariants/?datasetIds=",dataset,"&biosampleIds=",biosample_id)
    encoded_url <- URLencode(url)
    result <- read.table(encoded_url, header = TRUE, sep="\t",quote="")
    if (dim(result)[1] == 0) return(NA)
    col <- c("start","end","log2")
    suppressWarnings(result[,col] <- sapply(result[,col], as.numeric))
    result <- result[order(result$start),]
    chr <- result[,2]
    chr[which(chr == 'X')] <- 23
    chr[which(chr == 'Y')] <- 24
    chr <- as.integer(chr)
    result <- result[order(chr),]
    result <- result[order(result$biosample_id),]
    return(result)
}

read_variant_beacon <- function(biosample_id, domain, dataset){
    url <- paste0(domain,"/beacon/biosamples/",biosample_id,"/g_variants","?datasetIds=",dataset)
    encoded_url <- URLencode(url)
    response <- NULL
    tryCatch({
      response <- GET(encoded_url)
      httr::stop_for_status(response)
    }, error = function(e) {
      NULL
    })
    
    if (is.null(response)) return(NULL)
    
    content <- content(response)
    result <- lapply(content$response$resultSets[[1]]$results,unlist)
    
    # when no variants for this id              
    if (length(result) == 0) return(NA)
    
    result <- lapply(result,function(x){
      var_meta <- x[c('caseLevelData.id','caseLevelData.biosampleId','caseLevelData.analysisId',
                      'variation.subject.sequence_id','variantInternalId','caseLevelData.info.cnvValue','variation.copyChange')]
      temp_row <- as.data.frame(matrix(var_meta,nrow=1,byrow = TRUE))
      colnames(temp_row) <- c("variant_id","biosample_id","analysis_id",
                              "reference_genome","variant","variant_log2","variant_copychange")
      return(temp_row)
    })
    result <- Reduce(rbind,result)
    
    # change interval order 
    location <- strsplit(result$variant,split = ':')
    start <- sapply(location, function(x){x[2]})
    start <- as.numeric(gsub('-.*','',start))
    result <- result[order(start),]
    location <- strsplit(result$variant,split = ':')
    chr <-  sapply(location, function(x){x[1]})
    chr[chr == 'X'] <- 23
    chr[chr == 'Y'] <- 24
    chr <- as.numeric(chr)
    result <- result[order(chr),]
    
    # in case one sample corresponds to multiple analysis
    result <- result[order(result$analysis_id),]
    return(result)
}


read_variant_pgxseg_meta <- function(biosample_id,domain, dataset){
    url <- paste0(domain,"/services/pgxsegvariants/?datasetIds=",dataset,"&biosampleIds=",biosample_id)
    encoded_url <- URLencode(url)
    meta <- readLines(encoded_url)
    idx <- length(grep("#",meta))
    return(meta[seq_len(idx)])
}
  
read_cov_json <- function(url,codematches=FALSE,all_biosample_id=NULL){
  encoded_url <- URLencode(url)
  data <- content(GET(encoded_url))
  sample <- unlist(lapply(data$response$resultSets[[1]]$results, function(x){
    return(x$biosampleId)
  }))
  
  if (codematches){
    sel_sample_idx <- sample %in% all_biosample_id 
    if (sum(sel_sample_idx) == 0 & length(sample) > 0){
      warning("\n The option `codematches=TRUE` filters out all samples accessed by filters \n")
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
    return( t(data.frame(val,row.names=names)))
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



disease_code_check <- function(id, remain.id, url){
  total.id <- c()
  for (code in c('NCIT','pgx:icdom','pgx:icdot','UBERON')){
    idx <- grep(code,id)
    if (length(idx) > 0){
      encoded_url <- URLencode(paste0(url,code))
      res <- content(GET(encoded_url))
      total.id <- c(total.id,unlist(lapply(res$response$results, function(x){x$id})))
      remain.id <- remain.id [!remain.id %in% id[idx]]
    }
  }
  return(list(total=total.id,remain=remain.id))
}

pgidCheck <- function(id,domain,dataset){
    # this query doesn't work for individual GSM id and age filters
    geogsm.idx <- grep('geo:GSM',id)
    age.idx <- grep('age:',id)
    pass.idx <- c(geogsm.idx,age.idx)
    if (length(pass.idx > 0)){
      remain.id <- id[-pass.idx]
    } else{
      remain.id <- id
    }
    
    url <- paste0(domain,"/services/collations?datasetIds=",dataset,"&filters=")
    
    # check if disease_code associated filters exist in Progenetix
    check.res <- disease_code_check(id, remain.id, url)
    remain.id <- check.res$remain
    total.id <- check.res$total
    # check if non disease_code associated filters exist in Progenetix
    if (length(remain.id) > 0){
      encoded_url <- URLencode(paste0(url,transform_id(remain.id)))
      res <- content(GET(encoded_url))
      total.id <- c(total.id, unlist(lapply(res$response$results, function(x){x$id})))
    }
    
    return(id %in% c(total.id,id[pass.idx]))
}

pgxFreqLoader <- function(output, codematches, filters, domain, dataset) {
    # check if filters exists
    idcheck <- pgidCheck(filters,domain,dataset)
    if (!all(idcheck)){
        databasename <- ifelse(dataset=="cellz","cancercelllines","Progenetix")
        warning("\n No results for filters ", filters[!idcheck], " in ", databasename, " database.","\n")
        filters <- filters[idcheck]
    }
    # start query
    url <- paste0(domain,"/services/intervalFrequencies/?datasetIds=",dataset,"&output=",output)
  
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
    pg.data  <- read.table(encoded_url, header=TRUE, sep="\t",,quote="")
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

pgxmetaLoader <- function(type, biosample_id, individual_id, filters, codematches, skip, limit, filterLogic, domain, dataset){
    if (!(is.null(filters))){
        # check if filters exists
        idcheck <- pgidCheck(filters,domain,dataset)
        if (!all(idcheck)){
            databasename <- ifelse(dataset=="cellz","cancercelllines","Progenetix")
            warning("\n No results for filters ", filters[!idcheck], " in ", databasename, " database.","\n")
            filters <- filters[idcheck]
        }
        
        if (filterLogic == "AND"){
            filters <- transform_id(filters)
            url <- paste0(domain,"/services/sampletable/?datasetIds=",dataset,"&filters=",filters,"&responseEntityPathId=",type,"s")
            url  <- ifelse(is.null(limit), url, paste0(url,"&limit=",limit))
            url  <- ifelse(is.null(skip), url, paste0(url,"&skip=",skip))
            encoded_url <- URLencode(url)
            attempt::try_catch(
              res_1 <- read.table(encoded_url,stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header=TRUE,quote=""),
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
              url <- paste0(domain,"/services/sampletable/?datasetIds=",dataset,"&filters=",trans.filters[i],"&responseEntityPathId=",type,"s")
              url  <- ifelse(is.null(limit), url, paste0(url,"&limit=",limit))
              url  <- ifelse(is.null(skip), url, paste0(url,"&skip=",skip))
              encoded_url <- URLencode(url)
      
              attempt::try_catch(
                res_1[[i]] <- read.table(encoded_url,stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header=TRUE,quote=""),
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
        filter <- transform_id(biosample_id)
        url <- paste0(domain,"/services/sampletable/?datasetIds=",dataset,"&biosampleIds=",filter,"&responseEntityPathId=",type,"s")
        encoded_url <- URLencode(url)
        attempt::try_catch(res_2 <- read.table(encoded_url,stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header=TRUE,quote=""),.e= function(e){
                warning("\n Query fails for biosample_id ", filter, "\n")
        })           
    }

    if (!(is.null(individual_id))){
        filter <- transform_id(individual_id)
        url <- paste0(domain,"/services/sampletable/?datasetIds=",dataset,"&individualIds=",filter,"&responseEntityPathId=",type,"s")
        encoded_url <- URLencode(url)
        attempt::try_catch(res_3 <- read.table(encoded_url,stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header=TRUE,quote=""),.e= function(e){
                warning("\n Query fails for individual_id ", filter, "\n")
        })
    }

    if(!(exists('res_1') | exists('res_2')| exists('res_3'))) stop("No data retrieved")

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
    
    if (dim(res)[1] == 0) stop("No data retrieved")
    
    rownames(res) <- seq(dim(res)[1])
    if (codematches){
        if (type == "biosample"){
            idx <- res$biosample_id %in% biosample_id | res$individual_id %in% individual_id | 
            res$histological_diagnosis_id %in% filters | res$sampled_tissue_id %in% filters | 
            res$icdo_morphology_id %in% filters | res$icdo_topography_id %in% filters 
        } else if (type == "individual"){
            idx <- res$histological_diagnosis_id %in% filters | res$individual_id %in% individual_id
            if (!is.null(biosample_id)){
              warning("\n The option `codematches=TRUE` filters out samples accessed by biosample_id \n")
            }
        }     
    
        res <- res[idx,]
        if (dim(res)[1] == 0){
            warning("\n The option `codematches=TRUE` filters out all samples \n")
        }
    }
    
    res <- res[!duplicated(res),]
    return(res)
}

retry_query <- function(results,retry_query_idx,biosample_id,domain,dataset,num_cores){
    retry_count <- 1
    while (length(retry_query_idx) > 0 & retry_count <= 5){
      re_try_result <- parallel::mclapply(biosample_id[retry_query_idx],
                                        FUN = function(i){read_variant_beacon(i,domain, dataset)}, 
                                        mc.cores = num_cores)
      results[retry_query_idx] <-  re_try_result
      retry_query_idx <-  which(sapply(results,is.null))
      retry_count <- retry_count + 1
    }
    return(list(results=results,retry_query_idx= retry_query_idx))
}

pgxVariantLoader <- function(biosample_id, output, save_file, filename, domain, dataset, num_cores){
    if (save_file & is.null(output)){
      stop("The parameter 'output' is invalid when 'save_file=TRUE' (available: \"seg\" or \"pgxseg\")") 
    } 
    
    num_cores <- min(num_cores, parallelly::availableCores())
    
    if (!(is.null(output))){
        results <- parallel::mclapply(biosample_id,
                                      FUN = function(i){read_variant_pgxseg(i,domain, dataset)}, 
                                      mc.cores = num_cores)
    }else{
        results <- parallel::mclapply(biosample_id,
                                      FUN = function(i){read_variant_beacon(i,domain, dataset)}, 
                                      mc.cores = num_cores)
        
        
        retry_query_idx <- which(sapply(results,is.null))
        
        if (length(retry_query_idx) > 0){
            retry_result <- retry_query(results,retry_query_idx,biosample_id,domain,dataset,num_cores)
            retry_query_idx <- retry_result[['retry_query_idx']]
            results <- retry_result[['results']]
            if (length(retry_query_idx) > 0) warning("\n Query fails for biosample_id ", paste(biosample_id[retry_query_idx],collapse = ','), "\n")
        }
    }
    fail_idx <- which(is.na(results))
    if (length(fail_idx > 0)) warning("\n Query fails for biosample_id ", paste(biosample_id[fail_idx],collapse = ','), "\n")
    
    results[is.na(results)] <- NULL
    results <- do.call(rbind, results)
        
    if (dim(results)[1] == 0){
        warning("\n No variants retrieved \n")
        return()
    }
        
    # quality check
    results <- results[results$biosample_id %in% biosample_id,]
    rownames(results) <- seq(nrow(results))
    
    # format conversion
    if (!(is.null(output))){
      if (output == 'seg'){
        results <- results[,c(1,2,3,4,6,5)]
      }
    }
    
    if (save_file){
        if (output=='pgxseg'){
           # query metadata
            meta <-  parallel::mclapply(biosample_id,
                                      FUN = function(i){read_variant_pgxseg_meta(i,domain, dataset)}, 
                                      mc.cores = numCores)
            head <- meta[[1]][c(1:2)]
            meta <- lapply(meta, FUN = function(x){return(x[-c(1,2,3)])})
            meta <- do.call(c,meta)
            meta <- c(head,meta)
            # write result
            if (is.null(filename)) filename <- "variants.pgxseg"            
            write.table(meta, file=filename,row.names = FALSE,col.names = FALSE, quote = FALSE)
            suppressWarnings(write.table(results, append=TRUE, sep='\t',file=filename,row.names = FALSE,col.names = TRUE, quote = FALSE))
        } else if(output == 'seg'){
            if (is.null(filename))  filename <- "variants.seg"  
            write.table(results, file=filename, sep='\t',row.names = FALSE,col.names = TRUE, quote = FALSE)
        } 
        message("\n The file is saved \n")
        return()
    }

    return(results)
}


pgxcallsetLoader <- function(biosample_id, individual_id, filters, limit, skip, codematches, domain, dataset){
    pg.data <- list()
    if (!is.null(filters)){
      url  <- paste0(domain,"/services/samplematrix/?datasetIds=",dataset,'&filters=',filters)
      url  <- ifelse(is.null(limit), url, paste0(url,"&limit=",limit))
      url  <- ifelse(is.null(skip), url, paste0(url,"&skip=",skip))
      encoded_url <- URLencode(url)
      pg.data[[1]] <- read.table(encoded_url, header=TRUE, sep="\t",quote="")
      ori.dim <- dim(pg.data[[1]])[1]
      if (codematches){
        pg.data[[1]] <- pg.data[[1]][pg.data[[1]]$group_id %in% filters,]
        if (ori.dim > 0 & dim(pg.data[[1]])[1] == 0){
          warning("\n The option `codematches=TRUE` filters out all samples accessed by filters \n")
        }}
    } 
    
    if (!is.null(biosample_id)){
      url  <- paste0(domain,"/services/samplematrix/?datasetIds=",dataset,"&biosampleIds=",transform_id(biosample_id))
      encoded_url <- URLencode(url)
      pg.data[[2]] <- read.table(encoded_url, header=TRUE, sep="\t",quote="")
    }  
    
    if (!is.null(individual_id)){
      url  <- paste0(domain,"/services/samplematrix/?datasetIds=",dataset,"&individualIds=",transform_id(individual_id))
      encoded_url <- URLencode(url)
      pg.data[[3]] <- read.table(encoded_url, header=TRUE, sep="\t",quote="")
    }

    pg.data <- do.call(rbind,pg.data)
    # remove automatic prefix X 
    colnames(pg.data) <- gsub("X","",colnames(pg.data))
    # add chr prefix to avoid colnames with numeric start
    colnames(pg.data)[4:ncol(pg.data)] <- paste0("chr",colnames(pg.data)[4:ncol(pg.data)])
    # recover X chr    
    colnames(pg.data) <- gsub("chr\\.","chrX\\.",colnames(pg.data))

    return(pg.data)
}


pgxCovLoader <- function(biosample_id, individual_id, filters, codematches, skip, limit, domain, dataset){
    pg.data <- list()
    if (!is.null(filters)){
        url <- paste0(domain,"/beacon/analyses/?datasetIds=",dataset,"&output=cnvstats&filters=",filters)
        url  <- ifelse(is.null(limit), url, paste0(url,"&limit=",limit)) 
        url  <- ifelse(is.null(skip), url, paste0(url,"&skip=",skip))
        if (codematches){
          suppressWarnings(all_biosample_id <- pgxmetaLoader(type = 'biosample', 
                                                             biosample_id = NULL,individual_id = NULL,
                                                             filters = filters,codematches = TRUE,
                                                             skip=NULL,limit=0,filterLogic="AND",domain=domain,dataset=dataset))
          all_biosample_id <- all_biosample_id$biosample_id    
        }
        pg.data[[1]] <- read_cov_json(url,codematches,all_biosample_id)
    } 
  
    if (!is.null(biosample_id)){
        url  <- paste0(domain,"/beacon/analyses/?datasetIds=",dataset,"&output=cnvstats&biosampleIds=",transform_id(biosample_id))
        pg.data[[2]] <- read_cov_json(url)
    }
    
    if (!is.null(individual_id)){
        url  <- paste0(domain,"/beacon/analyses/?datasetIds=",dataset,"&output=cnvstats&&individualIds=",transform_id(individual_id))
        pg.data[[3]] <- read_cov_json(url)
    }
    
    result <- list()
    result$chrom_arm_coverage <- do.call(rbind,lapply(pg.data,function(x){return(x$chrom_arm_coverage)}))
    result$whole_chrom_coverage <- do.call(rbind,lapply(pg.data,function(x){return(x$whole_chrom_coverage)}))
    result$whole_genome_coverage <- do.call(rbind,lapply(pg.data,function(x){return(x$whole_genome_coverage)}))
    
    return(result)
}



