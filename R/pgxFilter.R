#' Query available filters
#'
#' This function retrieves available filters in the Progenetix database via the Beacon v2 API.
#'
#' @param prefix A string specifying the prefix of filters, such as 'NCIT'. Default is NULL, which 
#' means that all available filters will be returned. When specified, it returns all filters with the specified prefix.
#' @param return_all_prefix A logical value determining whether to return all valid prefixes of filters used in database. 
#' If TRUE, the `prefix` parameter will be ignored. Default is FALSE.
#' @param domain A string specifying the domain of the query data resource. Default is "http://progenetix.org".
#' @param entry_point A string specifying the entry point of the Beacon v2 API. Default is "beacon", resulting in the endpoint being "http://progenetix.org/beacon".
#' @param dataset A string specifying the dataset to query. When the parameter `domain` is "http://progenetix.org", available options are "progenetix" (by defualt) and "cancercelllines".
#' @importFrom httr GET content
#' @return Filter terms used in the data resource that you query.
#' @export
#' @examples
#' pgxFilter(prefix = "NCIT")

pgxFilter <- function(prefix=NULL, return_all_prefix=FALSE, domain="http://progenetix.org", entry_point = "beacon", dataset=NULL){

    url <- paste0(domain,"/",entry_point, "/filtering_terms")
    url <- ifelse(is.null(dataset), url, paste0(url,"?datasetIds=", transform_dataset_parameter(domain,dataset)))
    query  <- content(GET(url))
    
    data_lst <- lapply(query$response$filteringTerms,FUN = function(term){term$id})
    
    data <- do.call(c,data_lst)
    
    if (return_all_prefix){
        all_ids <- data
        # progenetix specific filter prefix 
        pgx_prefix <- NULL
        if (domain %in% c("http://progenetix.org","progenetix.org")){
          pgx_idx <- grepl("pgx:",all_ids)
          pgx_prefix <- all_ids[pgx_idx]
          pgx_prefix <- unique(gsub("\\-.*","",pgx_prefix))
          all_prefix <- all_ids[!pgx_idx]
        # general filter prefix
        } else{
          all_prefix <- all_ids
        }
        all_prefix <- unique(gsub(":.*","",all_prefix))
        all_prefix <- c(all_prefix,pgx_prefix)
        return(all_prefix)
    }
    
    if (!is.null(prefix)){
      sel_data <- data[grep(prefix,data)]
      return(sel_data)
    }
    
    return(data)
}




