#' Return object with selected filters
#'
#' This function helps select filters of interests from input object.
#'
#' @param data CNV frequency data returned by `pgxLoader` function
#' @param filters A single or a comma-concatenated list of identifiers such as c("NCIT:C7376","icdom-98353")
#' @return A list object with the same format of input but only containing selected filters
#' @export
#' @examples
#' total_freq <- pgxLoader(type="frequency", output ='pgxmatrix', filters=c("NCIT:C3512","NCIT:C3059","NCIT:C3716","NCIT:C4917"))
#' subset_freq <-  pgxFreqselect(data=total_freq, filters = c("NCIT:C3512","NCIT:C3059"))

pgxFreqselect  <- function(data, filters){
    meta_idx <- data$meta$code %in% filters
    new_total_count <- sum(data$meta[meta_idx,3])
    meta_idx[length(meta_idx)] <- TRUE
    res <- list(meta=data$meta[meta_idx,],data=NULL)
    res$meta[res$meta$code == 'total',3] <- new_total_count
    rownames(res$meta) <- seq_len(dim(res$meta)[1])
    
    mat <- data$data[[filters[1]]]
    res$data[[filters[1]]] <- data$data[[filters[1]]]
    if (length(filters) > 1){
        for (i in c(2:length(filters))){
            mat <- rbind(mat,data$data[[filters[i]]])
            res$data[[filters[i]]] <- data$data[[filters[i]]]
        }
    }
    rownames(mat) <- seq_len(dim(mat)[1])
    res$data[['total']]<- mat
    return(res)
}
