#' Return object with selected filters
#'
#' This function helps select filters of interests from input object.
#'
#' @param data CNV frequency data returned by `pgxLoader` function
#' @param filters A single or a comma-concatenated list of identifiers such as c("NCIT:C7376","icdom-98353")
#' @return A list object with the same format of input but only containing selected filters
#' @export


pgxFreqselect  <- function(
    data,
    filters){
    meta_idx <- data$meta$code %in% filters
    res <- list(meta=data$meta[meta_idx,],data=NULL)
    mat <- data$data[[filters[1]]]
    res$data[[filters[1]]] <- data$data[[filters[1]]]
    if (length(filters) > 1){
        for (i in c(2:length(filters))){
            mat <- rbind(mat,data$data[[filters[i]]])
            res$data[[filters[i]]] <- data$data[[filters[i]]]
        }
    }
    rownames(mat) <- seq(1:dim(mat)[1])
    res$data[['total']]<- mat
    rownames(res$meta) <- seq(1:dim(res$meta)[1])
    return(res)
}
