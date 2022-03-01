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
    header <- gsub('\\.\\.','_',header)
    header <- gsub('___','_',header)
    header <- gsub('::','_',header)
    colnames(data) <- header
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
   stop_if(.x=length(filters) < 1, msg="\n at least one valid filter has to be provided \n")
    # check output format
    stop_if(.x=!(output %in% c('pgxseg','pgxmatrix')),msg="\n Output is invalid. Only support 'pgxseg' or 'pgxmatrix' \n")
    # check if filters exists
    idcheck <- pgidCheck(filters)

    if (!all(idcheck)){
        cat(" No results for id", filters[!idcheck], "in progenetix database.\n","\n Only query id:", filters[idcheck],'\n')
        filters <- filters[idcheck]
    }
    # start query
    cat("\n accessing", "IntervalFrequencies service","from Progenetix \n")

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

pgxSampleLoader <- function(biosample_id,filters,codematches){
    cat("\n accessing", "biosample information","from Progenetix \n\n")

    if (!(is.null(filters))){
        for (i in c(1:length(filters))) {
            if_next <- FALSE
            url <- paste0("http://progenetix.org/cgi/bycon/beaconServer/biosamples.py?filters=",
                              filters[i],"&output=table")
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

    # remove nonsense items
    if (any((res$id %in% c("DUO:0000004","" )))){
        res <- res[-(which(res$id %in% c("DUO:0000004","" ))),]
    }
    # remove duplicated items
    if (any(duplicated(res$id))){
        res <- res[-(which(duplicated(res$id))),]}

    if (codematches){
        idx <- res$histological_diagnosis__id %in% filters | res$sampled_tissue__id %in% filters | res$icdo_morphology__id %in% filters |
        res$icdo_topography__id %in% filters | res$external_references__id_PMID %in% filters | res$id %in% biosample_id | 
        res$external_references__id_geo.GSM %in% filters | res$external_references__id_geo.GSE %in% filters | 
        res$external_references__id_geo.GPL %in% filters | res$external_references__id_cellosaurus %in% filters | 
        res$external_references__id_arrayexpress %in% filters

        res <- res[idx,]
        if (dim(res)[1] == 0){
            cat("Warning: the option `codematches=TRUE` filters out all samples \n")
        }}
    
    colnames(res)[c(1,3)] <- c('biosample_id','callset_id')
    return(res)
}

pgxVariantLoader <- function(biosample_id, output, save_file,filename){
    
    if (save_file){
        stop_if(.x= !(output %in% c("pgxseg","seg")), msg = "output is invalid ('seg' or 'pgxseg')")
    }

    cat("accessing", "variant information","from Progenetix \n")

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
                var_meta <- x[c('position.assemblyId','variantInternalId','position.refseqId',
                       'position.start','position.end','variantType','referenceBases','alternateBases')]
                id_ind <-  which(names(x) =='caseLevelData.id')
                temp_row <- c()
                for (i in id_ind){
                    temp_row <- rbind(temp_row,c(x[i],x[i-1],x[i-2]))}
                var_meta <- matrix(rep(var_meta,dim(temp_row)[1]),nrow=dim(temp_row)[1],
                                   byrow = T)
                temp_row <- as.data.frame(cbind(temp_row,var_meta))
                colnames(temp_row) <- c("variant_id","biosample_id","analysis_id",
                                        "assemply","interval","chromosome","start",
                                        "end","variant_type","reference_bases","alternate_bases")
                return(temp_row)
            })
            temp <- Reduce(rbind,temp)
            temp$chromosome  <- gsub('chr','',temp$chromosome)
            temp$start <- as.numeric(temp$start)
            temp$end <- as.numeric(temp$end)
            temp <- temp[order(temp$start),]
            chr <- temp$chromosome
            chr[which(chr == 'X')] <- 23
            chr[which(chr == 'Y')] <- 24
            chr <- as.integer(chr)
            temp <- temp[order(chr),]
            temp <- temp[order(temp$biosample_id),]
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
  cat("\n accessing", "CNV coverage profiles","from Progenetix \n")
  pg.url <- "http://www.progenetix.org/beacon/callsets/?output=pgxmatrix"
  
  filter <- transform_id(filters)
  pg.url  <- paste(pg.url,'&filters=',filter,sep="")
  
  pg.url  <- ifelse(is.null(limit), pg.url, paste0(pg.url,"&limit=",limit))
  pg.url  <- ifelse(is.null(skip), pg.url, paste0(pg.url,"&skip=",skip))
  
  meta <- readLines(pg.url,n=7)
  meta <-  gsub("#meta=>\"","",meta[7])
  meta <-  gsub("\"","",meta)
  cat(paste("\n", meta, "\n"))
  pg.data  <- read.table(pg.url, header=T, sep="\t", na="NA")
  
  if (codematches){
    pg.data <- pg.data[pg.data$group_id %in% filters,]
    if (dim(pg.data)[1] == 0){
      cat("\n Warning: the option `codematches=TRUE` filters out all samples \n")
    }}
  
  return(pg.data)
}
