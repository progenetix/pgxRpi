#' Count samples in one collation of a given id
#'
#' This function returns the number of samples for every group id in 
#' Progenetix database.
#'
#' @param group_id A single or a comma-concatenated list of identifiers such as c("NCIT:C7376","icdom-98353")
#' @return Count of samples in the given id
#' @export

pgxCount <- function(
    group_id=NULL) 
{
    if (is.null(group_id)){
        stop("Please input group id")
    }
    if (length(group_id) > 1){
        temp <- group_id[1]
        for (i in c(2:length(group_id))){
            temp <- paste(temp, group_id[i], sep=",")
        }
        id_name <- temp
    } else{
        id_name <- group_id
    }
    url <- paste0("https://progenetix.org/services/collations?filters=",id_name,"&method=counts&output=text")
    id_url <- url(description=url,
                  open='r')
    info <- read.table(id_url, header=F, sep="\t", na="NA")
    close(id_url)
    print(paste('Trying url:',url,"\n"))
    if (is.null(info)){
        return("No samples with queried group id")
    }
  
    info <- info[unlist(info[1]) %in% group_id,]
    rownames(info) <- seq(1:dim(info)[1])
    
    if (dim(info)[1] < length(group_id)){
        cat("\n Attention: No results for some queried id")
    }
    return (info)
}

