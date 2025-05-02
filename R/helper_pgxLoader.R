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

check_pgx_domain <- function(domain, data_type){
    if (!domain %in% c("http://progenetix.org","progenetix.org","https://cancercelllines.org","cancercelllines.org")){
      stop(data_type," can only be accessed from progenetix.org or cancercelllines.org")
    }
}

camel_2_snake <- function(camelstring) {
    # Use gsub to replace uppercase letters with _ followed by the lowercase version
    snake_case <- gsub("([a-z])([A-Z])", "\\1_\\L\\2", camelstring, perl = TRUE)
    # Convert the entire string to lowercase for consistency
    return(tolower(snake_case))
}

# utility function for transforming beacon response -----------------------

extract_list_data <- function(original_column, data){
    result <- list()
  
    if (length(data) == 0){
        result[[original_column]] <- NA
        return(result)
    }
  
    if (length(data) == 1){
        if(is.null(names(unlist(data[[1]])))){
        result[[original_column]] <- paste(data[[1]],collapse = ",")
        return(result)
        } else{
            data <- data[[1]]
        }
    }
    
    for (i in seq_len(length(data))){
        subdata <- data[i]
        # expand column name
        single_colname <- original_column
        if (!is.null(names(subdata))){
            single_colname <- camel_2_snake(paste(single_colname,names(subdata),sep = "_"))
        }
        result <- append(result,extract_list_data(single_colname,subdata))
    }
  
    # merge duplicated fields
    dup_names <- unique(names(result)[duplicated(names(result))])
    if (length(dup_names) > 0){
        rm_idx <- c()
        for (name in dup_names){
            idx <- which(names(result) == name)
            result[[idx[1]]] <- paste0(unlist(result[idx]),collapse = ",")
            rm_idx <- c(rm_idx,idx[2:length(idx)])
        }
        result <- result[-rm_idx]
    }
  
    return(result)
}

extract_general_results <- function(data, column, mapping){
    keys <- unlist(strsplit(mapping, "\\."))
    result <- data
    for (key in keys){
        # if key doesn't exist
        if (!key %in% names(result)){
            # in case JSON object
            if (length(result) == 1){
                result <- result[[1]]
            } else{
            # in case JSON array
                result <- unlist(result)
            }
          
            if (!key %in% names(result)){
                result <- NULL
                next
            }
            result <- unlist(lapply(which(names(result) == key), function(i){result[[i]]}))
        } else{
            result <- result[[key]]
        }
    }
    result <- extract_list_data(column,result)
    return(result)
}

extract_all_results <- function(data, mapping, type){
    result <- list()
    for (column in names(mapping)){
        col_data <- extract_general_results(data, column, mapping[[column]][[1]])
        result <- append(result,  col_data)
    }  
    return(as.data.frame(result))
}

extract_beacon_query <- function(url, type, dataset){
    mapping_rules <- yaml::read_yaml(system.file("config", "datatable_mappings.yaml", package = "pgxRpi"))
    mapping_rules <- mapping_rules[[type]]
    query <- GET(url)
    data <- content(query)
    data_df <- c()
    
    # beaconFilteringTermsResults response
    if (type == "filtering_terms"){
      data_info <- lapply(data$response$filteringTerms, function(x){extract_all_results(x, mapping_rules, type)})
      data_df <- do.call(rbind, data_info)
    } else if(type == "counts"){
    # beaconCountResponseSection
      if (!"responseSummary" %in% names(data)) return(NA)
      data_df <- extract_all_results(data$responseSummary, mapping_rules, type)
      data_df <- data_df$count
    } else{
    # beaconResultsets response 
      if (!data$responseSummary$exists) return(NULL)
      if (!is.null(dataset)){
        datasetids <- sapply(data$response$resultSets,function(x){x$id})
        id_idx <- which(datasetids %in% dataset)
        if (length(id_idx) == 0) stop("Data not found for the specified dataset ", dataset)
      } else{
        id_idx <- seq(length(data$response$resultSets))
      }
      
      for (i in id_idx){
        data_list <- data$response$resultSets[[i]]$results
        data_info <- lapply(data_list, function(x){extract_all_results(x, mapping_rules, type)})
        data_df <- dplyr::bind_rows(data_df,do.call(dplyr::bind_rows, data_info))
      }
    }

    return(data_df)
}

# function to query metadata for biosamples, individuals, analyses --------

metaLoader <- function(type, biosample_id, individual_id, filters, codematches, skip, limit, domain, entry_point, dataset){
  res <- c()
  # query by filters
  if (!(is.null(filters))){               
    url <- paste0(domain,"/",entry_point, "/", type, "?filters=",transform_id(filters),"&limit=",limit,"&skip=",skip)
    encoded_url <- URLencode(url)
    attempt::try_catch(
      res_1 <- extract_beacon_query(encoded_url, type, dataset),.e= function(e){
        warning("\n Query to ", domain, " fails for the filters ", paste(filters,collapse=",") ,"\n")
      }
    )
    if (exists('res_1')){res <- dplyr::bind_rows(res,res_1)}               
  }
  
  # query by biosample_id
  if (!(is.null(biosample_id))){
    biosample_ids <- transform_id(biosample_id)
    url <- paste0(domain,"/",entry_point, "/biosamples/", biosample_ids, "/", type)
    encoded_url <- URLencode(url)
    attempt::try_catch(
      res_2 <- extract_beacon_query(encoded_url, type, dataset),.e= function(e){
        warning("\n Query to ", domain, " fails for biosample_id ", paste(biosample_ids,collapse = ","), "\n")
      }
    )
    if (exists('res_2')){res <- dplyr::bind_rows(res,res_2)} 
  }
  
  # query by individual_id
  if (!(is.null(individual_id))){
    individual_ids <- transform_id(individual_id)
    url <- paste0(domain,"/",entry_point, "/individuals/", individual_ids, "/", type)
    encoded_url <- URLencode(url)
    attempt::try_catch(
      res_3 <- extract_beacon_query(encoded_url, type, dataset),.e= function(e){
        warning("\n Query to ", domain, " fails for individual_id ", paste(individual_ids,collapse = ","), "\n")
      }
    )
    if (exists('res_3')){res <- dplyr::bind_rows(res,res_3)} 
  }
  
  # query with no conditions
  if (is.null(filters) & is.null(individual_id) & is.null(biosample_id)){
    url <- paste0(domain,"/",entry_point, "/", type)
    encoded_url <- URLencode(url)
    attempt::try_catch(
      res_4 <- extract_beacon_query(encoded_url, type, dataset),.e= function(e){
        warning("\n Query to ", domain, " fails for ", type, "\n")
      }
    )
    if (exists('res_4')){res <- dplyr::bind_rows(res,res_4)} 
  }
  
  if (length(res) == 0) return(NA)
  
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

pgxmetaLoader <- function(type, biosample_id, individual_id, filters, codematches, skip, limit, domain, entry_point, dataset, num_cores){
    if (length(entry_point) == 1) entry_point <- rep(entry_point,length(domain))
    if (length(entry_point) != length(domain)) stop("The parameters 'domain' and 'entry_point' do not match")
        
    num_cores <- min(num_cores, parallel::detectCores() - 1)
    future::plan(future::multisession,workers = num_cores)
  
    results <- future.apply::future_lapply(seq_len(length(domain)),FUN = function(i){metaLoader(type, biosample_id, individual_id, filters, codematches, skip, limit, domain[i], entry_point[i], dataset)})
    
    fail_idx <- which(is.na(results))
    if (length(fail_idx) == length(results)) stop("No data retrieved")
    if (length(fail_idx > 0)) {
        warning("\n Query fails for domain ", paste(domain[fail_idx],collapse = ','), "\n")
        domain <- domain[-fail_idx]
        results <- results[-fail_idx]
    }
   
    if (length(results) == 1){
        results <- results[[1]]
    }else{
        names(results) <- domain
    }
    return(results)     
}

# utility function for transforming variants data -------------------------

## beacon response

read_variant_beacon <- function(biosample_id, limit, domain, entry_point, dataset){
    if (is.null(biosample_id)){
        url <- paste0(domain,"/",entry_point, "/g_variants")
    }else{
        url <- paste0(domain,"/",entry_point,"/biosamples/",biosample_id,"/g_variants","?limit=",limit)
    }
    
    encoded_url <- URLencode(url)
    # make query not broken
    result <- NA
    attempt::try_catch(
    result <- extract_beacon_query(encoded_url, type="g_variants", dataset),.e= function(e){NULL}
    )        
    if (is.null(result)) return(result)
    if (all(is.na(result))) return(result)

    # variantInternalId in progenetix is interval-based and used for sorting
    if (domain %in% c("http://progenetix.org","progenetix.org","https://cancercelllines.org","cancercelllines.org")){
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

read_variant_pgxseg <- function(biosample_id, output, limit, domain){
    url <- switch(output,
                  pgxseg= paste0(domain,"/services/pgxsegvariants/?biosampleIds=",biosample_id,"&limit=",limit),
                  seg= paste0(domain,"/services/variantsbedfile/?output=igv&biosampleIds=",biosample_id, "&limit=",limit))
    
    encoded_url <- URLencode(url)
    
    seg <- read.table(encoded_url, header = TRUE, sep="\t",quote="")
    if (dim(seg)[1] == 0) return(NA)
    col <- switch(output,
                  pgxseg=c(3,4,5),
                  seg=c(3,4,6))
    suppressWarnings(seg[,col] <- sapply(seg[,col], as.numeric))
    seg <- seg[order(seg[,3]),]
    chr <- seg[,2]
    chr[which(chr == 'X')] <- 23
    chr[which(chr == 'Y')] <- 24
    chr <- as.integer(chr)
    seg <- seg[order(chr),]
    seg <- seg[order(seg[,1]),]
    meta <- NULL
    if (output=="pgxseg"){
      meta <- readLines(encoded_url)
      idx <- length(grep("#",meta))
      if (length(idx) != 0)  meta <- meta[seq_len(idx)]
    }
    return(list(seg=seg, meta=meta))
}

# function to query variants ----------------------------------------------

pgxVariantLoader <- function(biosample_id, output, limit, save_file, filename, domain, entry_point, dataset, num_cores){
    if (length(domain) > 1 | length(entry_point) > 1) stop("This query only supports one domain")
    if (!(is.null(output))) check_pgx_domain(domain, "Variant data in non-beacon output format")

    # query with no condition
    if (is.null(biosample_id)){
        results <- read_variant_beacon(biosample_id, domain, entry_point, dataset)
        if (all(is.na(results))) stop("No data retrieved")
    # query by biosample id
    }else{
        num_cores <- min(num_cores, parallel::detectCores() - 1)
        future::plan(future::multisession,workers = num_cores)

        if (!(is.null(output))){
            results <- future.apply::future_lapply(biosample_id,FUN = function(i){read_variant_pgxseg(i, output, limit, domain)})
        }else{
            results <- future.apply::future_lapply(biosample_id,FUN = function(i){read_variant_beacon(i, limit, domain, entry_point, dataset)})
        }

        fail_idx <- which(is.na(results))
        if (length(fail_idx) == length(results)) stop("No data retrieved")
        if (length(fail_idx > 0)) warning("\n Query fails for biosample_id ", paste(biosample_id[fail_idx],collapse = ','), "\n")

        results[is.na(results)] <- NULL

        if (!(is.null(output))){
            meta <- lapply(results,FUN= function(x){x[["meta"]]})
            meta <- do.call(c,meta)
            results <- lapply(results,FUN= function(x){x[["seg"]]})
        }

        results <- do.call(dplyr::bind_rows, results)
        # if the query succeed but no data in database
        if (length(results) == 0) stop("No data retrieved")
        # quality check
        id_col <- which(colnames(results) == "biosample_id")
        if (length(id_col) == 0) id_col <- 1
        results <- results[results[,id_col] %in% biosample_id,]
        rownames(results) <- seq(nrow(results))
    }
    
    if (save_file){
        append <- FALSE
        # pgxseg format
        if (!is.null(output)){
            if (output=='pgxseg'){
                write.table(meta, file=filename,row.names = FALSE,col.names = FALSE, quote = FALSE)
                append <- TRUE
            }
        }
        suppressWarnings(write.table(results, append=append, sep='\t',file=filename,row.names = FALSE,col.names = TRUE, quote = FALSE))
        message("\n The file is saved \n")
        return()
    }

    return(results)
}

# function to query sample count -----------------------------

countLoader <- function(filters,domain,entry_point){
    individual_url <- paste0(domain,"/",entry_point, "/individuals?filters=")
    sample_url <- paste0(domain,"/",entry_point, "/biosamples?filters=")
    analyses_url <- paste0(domain,"/",entry_point, "/analyses?filters=")
  
    res <- list()
    for (i in seq_len(length(filters))){
        info1 <- extract_beacon_query(paste0(individual_url,filters[i]),"counts",NULL)
        info2 <- extract_beacon_query(paste0(sample_url,filters[i]),"counts",NULL)
        info3 <- extract_beacon_query(paste0(analyses_url,filters[i]),"counts",NULL)
        res[[i]] <- data.frame(filter=filters[i],entity=c("individuals","biosamples","analyses"),count=c(info1,info2,info3))
    } 
  
    res <- do.call(rbind,res)
    if (length(res) == 0) return(NA)
  
    return (res)
}

pgxCount <- function(filters,domain,entry_point,num_cores){
    if (length(entry_point) == 1) entry_point <- rep(entry_point,length(domain))
    if (length(entry_point) != length(domain)) stop("The parameters 'domain' and 'entry_point' do not match")
  
    num_cores <- min(num_cores, parallel::detectCores() - 1)
    future::plan(future::multisession,workers = num_cores)
  
    results <- future.apply::future_lapply(seq_len(length(domain)),FUN = function(i){countLoader(filters,domain[i],entry_point[i])})
  
    fail_idx <- which(is.na(results))
    if (length(fail_idx) == length(results)) stop("No data retrieved")
    if (length(fail_idx > 0)) {
        warning("\n Query fails for domain ", paste(domain[fail_idx],collapse = ','), "\n")
        domain <- domain[-fail_idx]
        results <- results[-fail_idx]
    }
  
    if (length(results) == 1){
        results <- results[[1]]
    }else{
        names(results) <- domain
    }
    return(results)     
}

# function to query cnv frequency -----------------------------------------

pgxFreqLoader <- function(output, filters, domain) {
    if (length(domain) > 1) stop("This query only supports one domain")
    check_pgx_domain(domain, "CNV frequency data")

    # start query
    url <- paste0(domain,"/services/intervalFrequencies?filters=",transform_id(filters))
    encoded_url <- URLencode(url)

    # access data 
    pg_data  <- content(GET(encoded_url))

    if (!pg_data$responseSummary$exists) stop("No data retrieved")
    if (pg_data$responseSummary$numTotalResults != length(filters)){
        pg_allids <- sapply(pg_data$response$results,function(x){x$groupId})
        warning("\n Query to ", domain, " fails for filter ", paste(filters[!filters %in% pg_allids],collapse = ','), "\n")
    }

    pg_data_lst <- lapply(pg_data$response$results, function(x){
        indmeta <- as.data.frame(x[c("groupId","label","sampleCount")])
        colnames(indmeta) <- c("filter","label","sample_count")
        indfreq <- lapply(x$intervalFrequencies,function(x){
            as.data.frame(x[c("referenceName","start","end","gainFrequency","gainHlfrequency","lossFrequency","lossHlfrequency")])
        })
        indfreq <- do.call(rbind,indfreq)
        return(list(meta=indmeta,data=indfreq))
    })
    meta <- do.call(rbind,lapply(pg_data_lst,function(x){x$meta}))
    
    # convert to RangedSummarizedExperiment
    if (output == 'pgxmatrix'){
        result_freq_lst <- lapply(pg_data_lst, function(x){
            lfrequency <- data.frame(c(x$data$gainFrequency,x$data$lossFrequency))
            colnames(lfrequency) <- x$meta$filter
            hfrequency <- data.frame(c(x$data$gainHlfrequency,x$data$lossHlfrequency))
            colnames(hfrequency) <- x$meta$filter
            return(list(lfreq=lfrequency,hfreq=hfrequency))
        })
        
        lfrequency <- do.call(cbind,lapply(result_freq_lst, function(x){x$lfreq}))
        hfrequency <- do.call(cbind,lapply(result_freq_lst, function(x){x$hfreq}))
        
        rowRanges <- pg_data_lst[[1]]$data[c("referenceName","start","end")]
        colnames(rowRanges)[1] <- 'chr'
        rowRanges <- rbind(rowRanges,rowRanges) 
        rowRanges$type <- rep(c("gain","loss"),each = nrow(pg_data_lst[[1]]$data))
        rowRanges <- GenomicRanges::GRanges(rowRanges)

        result <- SummarizedExperiment::SummarizedExperiment(assays=list(lowlevel_cnv_frequency=lfrequency,
                                                                         highlevel_cnv_frequency = hfrequency),
                                                             rowRanges=rowRanges, colData=meta)
    # convert to GRangesList
    } else if (output == "pgxfreq"){
        result_freq <- do.call(rbind,lapply(pg_data_lst,function(x){cbind(filter = x$meta$filter,x$data)}))
        colnames(result_freq)[2] <- 'chr'
        colnames(result_freq)[5:8] <- c("low_gain_frequency","high_gain_frequency","low_loss_frequency","high_loss_frequency")  
        result <- GenomicRanges::makeGRangesListFromDataFrame(result_freq,split.field = 'filter',keep.extra.columns=TRUE)   
        S4Vectors::mcols(result) <- meta
    }
    
    return(result)
}

# utility function for transforming cnv fraction data ---------------------

read_cnvstats_json <- function(url,codematches=FALSE,all_biosample_id=NULL){
    encoded_url <- URLencode(url)
    data <- content(GET(encoded_url))
    # this is progenetix/cellz specific data format, so the length of dataset is 1 and this function doesn't require dataset parameter  
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

    data_3 <- as.data.frame(Reduce(rbind,data_3))
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

pgxFracLoader <- function(biosample_id, individual_id, filters, codematches, skip, limit, domain){
    if (length(domain) > 1) stop("This query only supports one domain")
    check_pgx_domain(domain, "CNV fraction data")
  
    pg.data <- list()
    if (!is.null(filters)){
        url <- paste0(domain,"/services/cnvstats/?filters=",transform_id(filters),"&limit=",limit,"&skip=",skip)

        if (codematches){
          suppressWarnings(all_biosample_id <- pgxmetaLoader(type = 'biosamples', 
                                                             biosample_id = NULL,individual_id = NULL,
                                                             filters = filters,codematches = TRUE,
                                                             skip = NULL,limit = 0,domain = domain,entry_point = "beacon",dataset = NULL))
          all_biosample_id <- all_biosample_id$biosample_id    
        }
        pg.data[[1]] <- read_cnvstats_json(url,codematches,all_biosample_id)
    } 
  
    if (!is.null(biosample_id)){
        url  <- paste0(domain,"/services/cnvstats/?biosampleIds=",transform_id(biosample_id))
        pg.data[[2]] <- read_cnvstats_json(url)
    }
    
    if (!is.null(individual_id)){
        url <- paste0(domain,"/services/cnvstats/?individualIds=",transform_id(individual_id))
        pg.data[[3]] <- read_cnvstats_json(url)
    }
    
    result <- list()
    result$arm_cnv_frac <- do.call(rbind,lapply(pg.data,function(x){return(x$arm_cnv_frac)}))
    result$chr_cnv_frac <- do.call(rbind,lapply(pg.data,function(x){return(x$chr_cnv_frac)}))
    result$genome_cnv_frac <- do.call(rbind,lapply(pg.data,function(x){return(x$genome_cnv_frac)}))
    
    return(result)
}

# function to query sample callset -----------------------------

pgxcallsetLoader <- function(biosample_id, individual_id, filters, limit, skip, codematches, domain){
    if (length(domain) > 1) stop("This query only supports one domain")
    check_pgx_domain(domain, "pgxmatrix data")
  
    pg.data <- list()

    if (!is.null(filters)){
        url  <- paste0(domain,"/services/samplematrix/?filters=",transform_id(filters),"&limit=",limit,"&skip=",skip)
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
        encoded_url <- URLencode(url)
        pg.data[[2]] <- read.table(encoded_url, header=TRUE, sep="\t",quote="")
    }  
    
    if (!is.null(individual_id)){
      url  <- paste0(domain,"/services/samplematrix/?individualIds=",transform_id(individual_id))
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

