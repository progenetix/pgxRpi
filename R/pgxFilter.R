#' Query available filters
#'
#' This function retrieves available filters in the Progenetix database.
#'
#' @param prefix A string specifying the prefix of filters, such as 'NCIT' and 'PMID'. Default is NULL, which 
#' means that all available filters will be returned. When specified, it returns all filters with the specified prefix.
#' @param return_all_prefix A logical value determining whether to return all valid prefixes of filters used in Progenetix. 
#' If TRUE, the `prefix` parameter will be ignored. Default is FALSE.
#' @param domain A string specifying the domain of the Progenetix database. Default is "http://progenetix.org".
#' @param dataset A string specifying the dataset to query. Default is "Progenetix". Other available options are "cancercelllines".
#' @importFrom httr GET content
#' @return filter terms used in Progenetix.
#' @export
#' @examples
#' pgxFilter(prefix = "NCIT")

pgxFilter <- function(prefix=NULL, return_all_prefix=FALSE, domain="http://progenetix.org", dataset="Progenetix"){
    dataset <- match.arg(dataset, c("Progenetix","cancercelllines","cellz","examplez"))
    # actual dataset name for cancercelllines
    if (dataset == "cancercelllines") dataset <- "cellz"
    url <- paste0(domain,"/beacon/filtering_terms?datasetIds=",dataset)
    query  <- content(GET(url))
    data_lst <- lapply(query$response$filteringTerms,FUN = function(term){
      term <- lapply(term,function(x){if(is.null(x)){x <- NA} else{x}})
      return(as.data.frame(term))})
    
    data <- do.call(plyr::rbind.fill,data_lst)
    
    if (return_all_prefix){
        all_prefix <- data$id
        all_prefix <- c(data$id[!grepl("pgx:cohort",data$id)],"pgx:cohort")
        all_prefix <- gsub("\\d","",all_prefix)
        all_prefix[!grepl("pgx:",all_prefix)] <- gsub(":.*","",all_prefix[!grepl("pgx:",all_prefix)])
        all_prefix[grepl("pgx:icdo",all_prefix)] <-  gsub("\\-.*","",all_prefix[grepl("pgx:icdo",all_prefix)])
        all_prefix <- unique(all_prefix)
        return(all_prefix)
    }
    
    if (!is.null(prefix)){
      sel_data <- data[grep(prefix,data$id),]
      return(sel_data$id)
    }
    
    return(data$id)
}




