# general utility function ------------------------------------------------

check_unused_parameters <- function(param, param_name, rightparam_name) {
  if (!is.null(param)) {
    rightparam <- paste(rightparam_name,collapse = " and ")
    warning("\n The parameter ", param_name, " is not used in this query. Only ",rightparam," are accepted. \n")
  }
}

check_missing_parameters <- function(param, param_name){
  if (length(param) < 1) {
    param_name <- paste(param_name,collapse = " or ")
    stop("\n The parameter ",param_name," is missing. \n")
  }
}

transform_id <- function(id){
    filter <- paste(id,collapse = ",")
    return(filter)
}

add_parameter <- function(url, param_name, param_value){
    url  <- ifelse(is.null(param_value), url, paste0(url,"&",param_name,"=",param_value))
    return(url)
}

transform_dataset_parameter <- function(domain, dataset){
    if (domain %in% c("http://progenetix.org","progenetix.org") & !is.null(dataset)){
         # actual dataset name for cancercelllines
         if (dataset == "cancercelllines") return("cellz")
    }
    return(dataset) 
}

# utility function for transforming beacon response -----------------------

extract_general_results <- function(data, mapping){
    keys <- unlist(strsplit(mapping, "\\."))
    result <- data
    for (key in keys){
        # if key doesn't exist
        if (!key %in% names(result)){
          # in case a list
            if (length(result) == 1){
                result <- result[[1]]
            } else{
          # in case an array
                result <- unlist(result)
            }
        if (!key %in% names(result)) return(NA)
        result <- result[which(names(result) == key)]
        } else{
            result <- result[[key]]
        }
    }
    if (length(result) == 0) return(NA)
    result <- paste0(result,collapse = ",")
    return(result)
}

extract_special_results <- function(data, mapping, prefix){
    keys <- unlist(strsplit(mapping, "\\."))
    result <- data[[keys[1]]]
    result <- sapply(result,function(x){x[[keys[2]]]})
    if (prefix == "pubmed") prefix <- "PMID"
    if (keys[1] == "externalReferences") result <- result[grepl(prefix,result,ignore.case = TRUE)]
    if (length(result) == 0) return(NA)
    result <- paste0(result,collapse = ",")
    return(result)
}

extract_all_results <- function(data, mapping, type){
    result <- list()
    for (column in names(mapping)){
        column_prefix <- unlist(strsplit(column, "_"))[1]
        # special cohort information for biosamples
        if (column_prefix %in% c("pubmed","cellosaurus","cbioportal","tcga","cohort")){
            result[[column]] <- extract_special_results(data,mapping[[column]][[1]],column_prefix)
        # general information for all entity types
        } else{
            result[[column]] <- extract_general_results(data, mapping[[column]][[1]])
        }
    }  
  return(as.data.frame(result))
}

extract_beacon_query <- function(url, type){
    mapping_rules <- yaml::read_yaml(system.file("config", "datatable_mappings.yaml", package = "pgxRpi"))
    mapping_rules <- mapping_rules[[type]]
    query <- GET(url)
    data <- content(query)
    if (!data$responseSummary$exists) return(NULL)
    data_list <- data$response$resultSets[[1]]$results
    data <- lapply(data_list, function(x){extract_all_results(x, mapping_rules, type)})
    data_df <- do.call(rbind, data)
    return(data_df)
}

# function to query metadata for biosamples, individuals, analyses --------

pgxmetaLoader <- function(type, biosample_id, individual_id, filters, codematches, skip, limit, domain, entry_point, dataset){
    res <- c()
    dataset <- transform_dataset_parameter(domain, dataset)
    # query by filters
    if (!(is.null(filters))){               
        url <- paste0(domain,"/",entry_point, "/", type, "?filters=",transform_id(filters))
        url <- add_parameter(url,"datasetIds", dataset)
        url <- add_parameter(url,"limit",limit)
        url <- add_parameter(url,"skip",skip)
        encoded_url <- URLencode(url)
        attempt::try_catch(
            res_1 <- extract_beacon_query(encoded_url, type),.e= function(e){
                warning("\n Query fails for the filters ", paste(filters,collapse=",") ,"\n")
            }
        )
        if (exists('res_1')){res <- rbind(res,res_1)}               
    }

     # query by biosample_id
    if (!(is.null(biosample_id))){
        biosample_ids <- transform_id(biosample_id)
        url <- paste0(domain,"/",entry_point, "/", type, "?biosampleIds=",biosample_ids)
        url <- add_parameter(url,"datasetIds",dataset)
        encoded_url <- URLencode(url)
        attempt::try_catch(
            res_2 <- extract_beacon_query(encoded_url, type),.e= function(e){
                warning("\n Query fails for biosample_id ", paste(biosample_ids,collapse = ","), "\n")
            }
        )
        if (exists('res_2')){res <- rbind(res,res_2)} 
    }

    # query by individual_id
    if (!(is.null(individual_id))){
        individual_ids <- transform_id(individual_id)
        url <- paste0(domain,"/",entry_point, "/", type,"?individualIds=",individual_ids)
        url <- add_parameter(url,"datasetIds",dataset)
        encoded_url <- URLencode(url)
        attempt::try_catch(
            res_3 <- extract_beacon_query(encoded_url, type),.e= function(e){
                warning("\n Query fails for individual_id ", paste(individual_ids,collapse = ","), "\n")
            }
        )
        if (exists('res_3')){res <- rbind(res,res_3)} 
    }
    
    if (length(res) == 0) stop("No data retrieved")
    
    rownames(res) <- seq(dim(res)[1])
    if (codematches){
        idx <- rownames(res)
        if (type == "biosamples"){
            idx <- res$biosample_id %in% biosample_id | res$individual_id %in% individual_id | 
            res$histological_diagnosis_id %in% filters | res$sampled_tissue_id %in% filters | 
            res$icdo_morphology_id %in% filters | res$icdo_topography_id %in% filters 
        } else if (type == "individuals"){
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

# utility function for transforming variants data -------------------------

## beacon response

read_variant_beacon <- function(biosample_id, domain, entry_point, dataset){
    url <- paste0(domain,"/",entry_point, "/g_variants","?biosampleIds=",biosample_id)
    url <- add_parameter(url,"datasetIds",dataset)
    encoded_url <- URLencode(url)

    # make query not broken
    result <- NA
    attempt::try_catch(
    result <- extract_beacon_query(encoded_url, type="g_variants"),.e= function(e){NULL}
    )        
    if (is.null(result)) return(result)
    if (all(is.na(result))) return(result)

    # variantInternalId in progenetix is interval-based and used for sorting
    if (domain %in% c("http://progenetix.org","progenetix.org")){
    # change interval order 
    location <- strsplit(result$variant_internal_id,split = ':')
    start <- sapply(location, function(x){x[2]})
    start <- as.numeric(gsub('-.*','',start))
    result <- result[order(start),]
    location <- strsplit(result$variant_internal_id,split = ':')
    chr <-  sapply(location, function(x){x[1]})
    chr[chr == 'X'] <- 23
    chr[chr == 'Y'] <- 24
    chr <- as.numeric(chr)
    result <- result[order(chr),]
    }

    # in case one sample corresponds to multiple analysis
    result <- result[order(result$analysis_id),]
    return(result)
}

## exported pgxseg data by bycon service 

read_variant_pgxseg <- function(biosample_id, domain, dataset){
    url <- paste0(domain,"/services/pgxsegvariants/?biosampleIds=",biosample_id)
    url <- add_parameter(url,"datasetIds",dataset)
    encoded_url <- URLencode(url)
    
    seg <- read.table(encoded_url, header = TRUE, sep="\t",quote="")
    if (dim(seg)[1] == 0) return(NA)
    col <- c("start","end","log2")
    suppressWarnings(seg[,col] <- sapply(seg[,col], as.numeric))
    seg <- seg[order(seg$start),]
    chr <- seg[,2]
    chr[which(chr == 'X')] <- 23
    chr[which(chr == 'Y')] <- 24
    chr <- as.integer(chr)
    seg <- seg[order(chr),]
    seg <- seg[order(seg$biosample_id),]

    meta <- readLines(encoded_url)
    idx <- length(grep("#",meta))
    meta <- meta[seq_len(idx)]

    return(list(seg=seg, meta=meta))
}

# function to query variants ----------------------------------------------

pgxVariantLoader <- function(biosample_id, output, save_file, filename, domain, entry_point, dataset, num_cores){
    dataset <- transform_dataset_parameter(domain, dataset)
    num_cores <- min(num_cores, parallel::detectCores() - 1)
    future::plan(future::multisession,workers = num_cores)
    
    if (!(is.null(output))){
        results <- future.apply::future_lapply(biosample_id,FUN = function(i){read_variant_pgxseg(i, domain, dataset)})
    }else{
        results <- future.apply::future_lapply(biosample_id,FUN = function(i){read_variant_beacon(i, domain, entry_point, dataset)})
    }
 
    fail_idx <- which(is.na(results))
    if (length(fail_idx) == length(results)) stop("Query fails for all samples")
    if (length(fail_idx > 0)) warning("\n Query fails for biosample_id ", paste(biosample_id[fail_idx],collapse = ','), "\n")
    
    results[is.na(results)] <- NULL
    
    if (!(is.null(output))){
      meta <- lapply(results,FUN= function(x){x[["meta"]]})
      head <- meta[[1]][1:2]
      meta <- lapply(meta, FUN = function(x){return(x[-c(1,2,3)])})
      meta <- do.call(c,meta)
      meta <- c(head,meta)
      results <- lapply(results,FUN= function(x){x[["seg"]]})
    }

    results <- do.call(rbind, results)
    # if the query succeed but no data in database
    if (is.null(results)) stop("No data retrieved")
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
      # pgxseg format
        if (!is.null(output)){
          if (output=='pgxseg'){
            # write result
            write.table(meta, file=filename,row.names = FALSE,col.names = FALSE, quote = FALSE)
            suppressWarnings(write.table(results, append=TRUE, sep='\t',file=filename,row.names = FALSE,col.names = TRUE, quote = FALSE))
          } 
      # tsv format
        } else {
            write.table(results, file=filename, sep='\t',row.names = FALSE,col.names = TRUE, quote = FALSE)
        } 
        message("\n The file is saved \n")
        return()
    }

    return(results)
}

# function to query cnv frequency -----------------------------------------

pgxFreqLoader <- function(output, filters, domain, dataset) {
    if (!domain %in% c("http://progenetix.org","progenetix.org")) stop("CNV frequency data can only be accessed from progenetix.org")

    # start query
    url <- paste0(domain,"/services/intervalFrequencies?output=",output)
    url <- add_parameter(url,"datasetIds",transform_dataset_parameter(domain,dataset))
    url <- add_parameter(url,"filters",transform_id(filters))
    encoded_url <- URLencode(url)

    # access data 
    pg.data  <- read.table(encoded_url, header=TRUE, sep="\t",,quote="")
    colnames(pg.data)[1] <- 'filters'

    filter_exist <- filters %in% unique(pg.data[,1])
    if (all(!filter_exist)) stop("No data retrieved")
    if (any(!filter_exist)) warning("\n Query fails for the filters ", paste(filters[!filter_exist],collapse = ','),"\n")
    filters <- filters[filter_exist]

    # access metadata from JSON output
    meta <- readLines(encoded_url)[seq_len(length(filters))+1]
    meta_lst <- unlist(strsplit(meta,split = ';'))
    label <- meta_lst[grep("label",meta_lst)]
    label<- gsub('.*=','',label)
    count <- meta_lst[grep("sample_count",meta_lst)]
    count <- as.numeric(gsub('.*=','',count))
    meta <- data.frame(filter = filters, label, sample_count=count)
    
    # convert to GRangesList
    if (output == 'pgxfreq'){
      colnames(pg.data)[2] <- 'chr'
      data_lst <- GenomicRanges::makeGRangesListFromDataFrame(pg.data,split.field = 'filters',keep.extra.columns=TRUE)
      S4Vectors::mcols(data_lst) <- meta
    
    # convert to RangedSummarizedExperiment
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
      
      data_lst <- SummarizedExperiment::SummarizedExperiment(assays=list(cnv_frequency=frequency),
                           rowRanges=rowRanges, colData=meta)
    }
    
    return(data_lst)
}

# utility function for transforming cnv fraction data ---------------------

read_cnvstats_json <- function(url,codematches=FALSE,all_biosample_id=NULL){
    encoded_url <- URLencode(url)
    data <- content(GET(encoded_url))
    sample <- unlist(lapply(data$response$resultSets[[1]]$results, function(x){
    return(x$biosampleId)
    }))

    # filter with exact code match
    if (codematches){
    sel_sample_idx <- sample %in% all_biosample_id 
    if (sum(sel_sample_idx) == 0 & length(sample) > 0){
      warning("\n The option `codematches=TRUE` filters out all samples accessed by filters \n")
      return()
    }
    } else{
    sel_sample_idx <- seq(length(sample))
    }

    # extract fraction
    data_2 <- lapply(data$response$resultSets[[1]]$results[sel_sample_idx], function(x){
    stats <- lapply(x$cnvChroStats, function(y){
      as.data.frame(y)
    })
    cnvstats <- Reduce(rbind,stats)
    rownames(cnvstats) <- names(stats)
    res <- list()
    res[[1]] <- cnvstats
    res[[2]] <- as.data.frame(x$cnvStats)
    res[[3]] <- x$id
    return(res)
    })

    # some results are SNV instead of CNV
    total_frac <- lapply(data_2, function(x){x[[2]]})
    rm_idx <- which(sapply(total_frac,length) == 0)
    if (length(rm_idx) > 0) data_2 <- data_2[-rm_idx]
   
    analyses_ids <- sapply(data_2, function(x){x[[3]]})

    # whole genome CNV fraction

    total_frac <- Reduce(rbind, total_frac)
    rownames(total_frac) <- analyses_ids
    total_frac <- total_frac[,c(2,6,4)]

    # chromosome & chromosome arm CNV fraction
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
    rownames(data_3) <- analyses_ids

    arm_frac <- data_3[,seq_len(96)]
    chrom_frac <- data_3[,c(97:144)]

    result <- list()
    result$arm_cnv_frac <- arm_frac
    result$chr_cnv_frac <- chrom_frac
    result$genome_cnv_frac <- total_frac
    return(result)
}

# function to query cnv fraction ------------------------------------------

pgxFracLoader <- function(biosample_id, individual_id, filters, codematches, skip, limit, domain, dataset){
    if (!domain %in% c("http://progenetix.org","progenetix.org")) stop("CNV fraction data can only be accessed from progenetix.org")
    dataset <- transform_dataset_parameter(domain, dataset)
    
    pg.data <- list()
    if (!is.null(filters)){
        url <- paste0(domain,"/beacon/analyses/?output=cnvstats&filters=",transform_id(filters))
        url <- add_parameter(url,"datasetIds",dataset)
        url <- add_parameter(url,"limit",limit)
        url <- add_parameter(url,"skip",skip)

        if (codematches){
          suppressWarnings(all_biosample_id <- pgxmetaLoader(type = 'biosamples', 
                                                             biosample_id = NULL,individual_id = NULL,
                                                             filters = filters,codematches = TRUE,
                                                             skip = NULL,limit = 0,domain = domain,entry_point = "beacon",dataset = dataset))
          all_biosample_id <- all_biosample_id$biosample_id    
        }
        pg.data[[1]] <- read_cnvstats_json(url,codematches,all_biosample_id)
    } 
  
    if (!is.null(biosample_id)){
        url  <- paste0(domain,"/beacon/analyses/?output=cnvstats&biosampleIds=",transform_id(biosample_id))
        url <- add_parameter(url,"datasetIds",dataset)
        pg.data[[2]] <- read_cnvstats_json(url)
    }
    
    if (!is.null(individual_id)){
        url <- paste0(domain,"/beacon/analyses/?output=cnvstats&individualIds=",transform_id(individual_id))
        url <- add_parameter(url,"datasetIds",dataset)
        pg.data[[3]] <- read_cnvstats_json(url)
    }
    
    result <- list()
    result$arm_cnv_frac <- do.call(rbind,lapply(pg.data,function(x){return(x$arm_cnv_frac)}))
    result$chr_cnv_frac <- do.call(rbind,lapply(pg.data,function(x){return(x$chr_cnv_frac)}))
    result$genome_cnv_frac <- do.call(rbind,lapply(pg.data,function(x){return(x$genome_cnv_frac)}))
    
    return(result)
}

# function to query sample callset -----------------------------

pgxcallsetLoader <- function(biosample_id, individual_id, filters, limit, skip, codematches, domain, dataset){
    if (!domain %in% c("http://progenetix.org","progenetix.org")) stop("pgxmatrix data can only be accessed from progenetix.org")
    dataset <- transform_dataset_parameter(domain, dataset)
    
    pg.data <- list()

    if (!is.null(filters)){
        url  <- paste0(domain,"/services/samplematrix/?filters=",transform_id(filters))
        url <- add_parameter(url,"datasetIds",dataset)
        url <- add_parameter(url,"limit",limit)
        url <- add_parameter(url,"skip",skip)        
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
        url  <- paste0(domain,"/services/samplematrix/?biosampleIds=",transform_id(biosample_id))
        url <- add_parameter(url,"datasetIds",dataset)
        encoded_url <- URLencode(url)
        pg.data[[2]] <- read.table(encoded_url, header=TRUE, sep="\t",quote="")
    }  
    
    if (!is.null(individual_id)){
      url  <- paste0(domain,"/services/samplematrix/?individualIds=",transform_id(individual_id))
      url <- add_parameter(url,"datasetIds",dataset)
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

     # convert to RangedSummarizedExperiment
   
    callset <- as.data.frame(t(pg.data[,-c(1,2,3)]))
    rowRanges <- colnames(pg.data)[4:ncol(pg.data)]
    rowRanges <- lapply(rowRanges,function(x){
        range_str <- unlist(strsplit(x,split = '\\.'))
        chr <- range_str[1]
        chr <- gsub("X","",chr)
        if (chr == "") chr <- 'X'
        return(data.frame(chr=chr,start=range_str[2],end = range_str[3],type=range_str[4]))
        })
    rowRanges <- do.call(rbind, rowRanges)
    rowRanges <- GenomicRanges::GRanges(rowRanges)
    names(rowRanges) <- seq_len(dim(callset)[1])
      
    colnames(callset) <- pg.data$analysis_id
    rownames(callset) <- seq_len(dim(callset)[1])
      
    result <- SummarizedExperiment::SummarizedExperiment(assays=list(cnv_matrix=callset),
                           rowRanges=rowRanges, colData=pg.data[,c(1,2,3)])
    
    return(result)
}

# function to query sample count -----------------------------

pgxCount <- function(filters=NULL,domain="http://progenetix.org",dataset=NULL){
    if (!domain %in% c("http://progenetix.org","progenetix.org")) stop("This function accesses sample count data from progenetix.org.")

    dataset <- transform_dataset_parameter(domain, dataset)
    
    filter <- transform_id(filters)
    url <- paste0(domain,"/services/collations?filters=",filter)
    url <- add_parameter(url,"datasetIds",dataset)
    encoded_url <- URLencode(url)
    info <-  content(GET(url))
    res <- lapply(info$response$results, function(x){
        if (is.null(x$label)) x$label <- NA
        df <- data.frame(filters=x$id,label=x$label,total_count=x$count,exact_match_count=x$codeMatches)
    })
    
    res <- Reduce(rbind,res)
    if (length(res) == 0) stop("No data retrieved")
  
    return (res)
}


