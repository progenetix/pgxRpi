#' Count samples in one collation of a given filter
#'
#' This function returns the number of samples for every filter in 
#' Progenetix database.
#'
#' @param filters A single or a comma-concatenated list of identifiers such as c("NCIT:C7376","icdom-98353")
#' @param domain A string specifying the domain of database. Default is "http://progenetix.org".
#' @param dataset A string specifying the dataset to query. Default is "progenetix". Other available options are "cancercelllines".
#' @importFrom utils URLencode
#' @importFrom httr GET content
#' @return Count of samples in the given filter
#' @export
#' @examples
#' pgxCount(filters = "NCIT:C3512")

pgxCount <- function(filters=NULL,domain="http://progenetix.org",dataset="progenetix"){
    if (is.null(filters)){
        stop("Please input filter")
    }

    dataset <- match.arg(dataset, c("progenetix","cancercelllines","cellz","examplez"))
    # actual dataset name for cancercelllines
    if (dataset == "cancercelllines") dataset <- "cellz"
    
    filter <- transform_id(filters)
    url <- paste0(domain,"/services/collations?datasetIds=",dataset,"&filters=",filter)
    encoded_url <- URLencode(url)
    info <-  content(GET(url))
    res <- lapply(info$response$results, function(x){
        if (is.null(x$label)) x$label <- NA
        df <- data.frame(filters=x$id,label=x$label,total_count=x$count,exact_match_count=x$codeMatches)
    })
    
    res <- Reduce(rbind,res)
    if (length(res) == 0){
        return("No samples with the queried filter \n")
    }
  
    return (res)
}

