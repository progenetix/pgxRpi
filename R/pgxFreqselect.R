#' Return object with selected ids
#'
#' This function helps select ids of interests from input object.
#'
#' @param data CNV frequency data returned by `pgxLoader` function
#' @param id A single or a comma-concatenated list of identifiers such as c("NCIT:C7376","icdom-98353")
#' @return A list object with the same format of input but only containing selected id
#' @export


pgxFreqselect  <- function(
    data,
    id){
    meta_idx <- data$meta$code %in% id
    res <- list(meta=data$meta[meta_idx,],data=NULL)
    mat <- data$data[[id[1]]]
    res$data[[id[1]]] <- data$data[[id[1]]]
    if (length(id) > 1){
        for (i in c(2:length(id))){
            mat <- rbind(mat,data$data[[id[i]]])
            res$data[[id[i]]] <- data$data[[id[i]]]
        }
    }
    rownames(mat) <- seq(1:dim(mat)[1])
    res$data[['total']]<- mat
    rownames(res$meta) <- seq(1:dim(res$meta)[1])
    return(res)
}
