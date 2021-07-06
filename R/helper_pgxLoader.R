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
    names <- read.table(url, nrow = 1, stringsAsFactors = FALSE, sep = "\t")
    names[c(24:31,38:39)] <- c("UBERON_code","UBERON_label", "NCIT_code", "NCIT_label",
                               "icdom_code", "icdom_label","icdot_code", "icdot_label","PMID","PMID_label")
    data <- read.table(url, skip = 1, stringsAsFactors = FALSE, sep = "\t",fill=TRUE)
    data <- data[, c(1:31, 38:39)]
    names(data) <- names[c(1:31, 38:39)]
    return(data.frame(data))
}

read_variant_pgxseg <- function(url){
    result <- read.table(url, header = T, sep="\t")
    col <- c("start","end","log2")
    suppressWarnings(result[,col] <- sapply(result[,col], as.numeric))
    return(result)
}

read_variant_pgxseg_meta <- function(url){
    meta <- readLines(url)
    idx <- length(grep("#",meta))
    return(meta[c(1:idx)])
}

pgidCheck <- function(id){
    id_prefix <- unique(sub(":.*", "", id))
    total_url <- url(description=paste0("https://progenetix.org/services/collations?filters=",id_prefix[1],"&method=codematches&output=text"),
                     open='r')
    info <- read.table(total_url, header=F, sep="\t", na="NA",fill=TRUE,quote = "")
    close(total_url)
    if (length(id_prefix) > 1){
        for (i in c(2:length(id_prefix))){
            total_url <- url(description=paste0("https://progenetix.org/services/collations?filters=",id_prefix[i],"&method=codematches&output=text"),
                             open='r')
            temp <- read.table(total_url, header=F, na="NA",sep = '\t',fill=TRUE,quote = "")
            close(total_url)
            info <- rbind(info, temp)
        }
    }

    return(id %in% unlist(info[1]))
}

pgxFreqLoader <- function(output, codematches, group_id) {
    # check output format
    stop_if(.x=!(output %in% c('pgxseg','pgxmatrix')),msg="\n Output is invalid. Only support 'pgxseg' or 'pgxmatrix' \n")
    # check if group id exists
    idcheck <- pgidCheck(group_id)

    if (!all(idcheck)){
        cat(" No results for id", group_id[!idcheck], "in progenetix database.\n","\n Only query id:", group_id[idcheck],'\n')
        group_id <- group_id[idcheck]
    }
    # start query
    cat("\n accessing", "IntervalFrequencies service","from Progenetix \n")

    pg.url <- paste0("http://www.progenetix.org/services/intervalFrequencies/?output=",output)


    ## check group id again
    stop_if(.x=length(group_id) < 1, msg="\n at least one valid group id has to be provided \n")

    if (length(group_id)>1){
        filter <- transform_id(group_id)
        pg.url  <- paste(pg.url,'&filters=',filter,sep="")
        } else{
        pg.url  <- paste(pg.url, '&', 'id', '=', group_id,sep="")
    }

    if (codematches){
        pg.url  <- paste0(pg.url, '&method=codematches')
    }

    meta <- readLines(pg.url)[c(1:(length(group_id)+1))]
    meta_lst <- unlist(strsplit(meta,split = ';'))
    label <- meta_lst[grep("label",meta_lst)]
    label<- gsub('.*=','',label)
    count <- meta_lst[grep("sample_count",meta_lst)]
    count <- as.numeric(gsub('.*=','',count))
    meta <- data.frame(code = group_id, label=label, sample_count=count)
    pg.data  <- read.table(pg.url, header=T, sep="\t", na="NA")
    data_lst <- list()
    for (i in group_id){
        data_lst[[i]] <- pg.data[pg.data$group_id == i,]
    }
    data_lst[['total']] <- pg.data
    result <- list(meta = meta, data = data_lst)

    return(result)
}

pgxSampleLoader <- function(biosample_id,group_id,codematches){
    cat("\n accessing", "biosample information","from Progenetix \n\n")

    if (!(is.null(group_id))){
        for (i in c(1:length(group_id))) {
            if_next <- FALSE
            url <- paste0("http://progenetix.org/cgi/bycon/beaconServer/biosamples.py?filters=",
                              group_id[i],"&output=table")
            if (!(exists('res_1'))){
                try_catch(res_1 <- read_sample(url),.e= function(e){
                    cat(paste("No samples with the group id", group_id[i],"\n"))
                    if_next <<- TRUE})
                if (if_next){ next }
            }else {
                try_catch(temp <- read_sample(url),.e= function(e){
                    cat(paste("No samples with the group id", group_id[i],"\n"))
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
        idx <- res$NCIT_code %in% group_id | res$UBERON_code %in% group_id | res$icdom_code %in% group_id |
            res$icdot_code %in% group_id | res$PMID %in% group_id | res$id %in% biosample_id
        res <- res[idx,]}
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
            temp <- lapply(temp$resultSets[[1]]$results,unlist)
            temp <- as.data.frame(dplyr::bind_rows(temp))
            temp <- temp[,c("id","biosampleId","callsetId","digest","info.varLength","referenceName",
                            "start","end","variantType")]
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
