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
    
    filter=transform_id(filters)
    url <- paste0("https://progenetix.org/services/collations?filters=",filter)
    info <- rjson::fromJSON(file = url)
    res <- lapply(info$response$results, function(x){
      data.frame(filters=x$id,label=x$label,total_count=x$count,exact_match_count=x$codeMatches)
      })
    
    res <- Reduce(rbind,res)
    if (length(res) == 0){
        return("No samples with the queried filter \n")
    }
  
    return (res)
}

