#' Count samples in one collation of a given filter
#'
#' This function returns the number of samples for every filter in 
#' Progenetix database.
#'
#' @param filters A single or a comma-concatenated list of identifiers such as c("NCIT:C7376","icdom-98353")
#' @return Count of samples in the given filter
#' @export

pgxCount <- function(
    filters=NULL) 
{
    if (is.null(filters)){
        stop("Please input filter")
    }
    if (length(filters) > 1){
        temp <- filters[1]
        for (i in c(2:length(filters))){
            temp <- paste(temp, filters[i], sep=",")
        }
        id_name <- temp
    } else{
        id_name <- filters
    }
    url <- paste0("https://progenetix.org/services/collations?filters=",id_name,"&method=counts&output=text")
    id_url <- url(description=url,
                  open='r')
    info <- read.table(id_url, header=F, sep="\t", na="NA")
    close(id_url)
    if (is.null(info)){
        return("No samples with the queried filter \n")
    }
  
    info <- info[unlist(info[1]) %in% filters,]
    rownames(info) <- seq(1:dim(info)[1])
    
    if (dim(info)[1] < length(filters)){
        cat("Attention: No results for some queried filters \n")
    }
    
    info <- info[,c(1,2,3,4)]
    colnames(info) <- c('filters','label','total','exact_match')
    return (info)
}

