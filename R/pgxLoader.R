#' Load data from Progenetix database 
#'
#' This function loads various data from `Progenetix` database.   
#'
#' @param type A string specifying output data type. Available options are "biosample", 
#' "variant" or "frequency". 
#' @param output A string specifying output file format. When the parameter `type` is "variant",
#' available options are NULL, "pgxseg" ,"pgxmatrix", "coverage" or "seg" ; When the parameter `type` is "frequency",
#' available options are "pgxseg" or "pgxmatrix" .
#' @param filters A single or a comma-concatenated list of identifiers for cancer type,
#' literature, and cohorts such as c("NCIT:C7376","pgx:icdom-98353","PMID:22824167", "pgx:cohort-TCGAcancers")
#' @param codematches A logical value determining whether to exclude samples 
#' from child terms of specified filters. If FALSE, retrieved samples include child 
#' terms of specified filters. Default is FALSE.
#' @param biosample_id  A single or a comma-concatenated list of identifiers used 
#' in Progenetix database for identifying biosamples. 
#' @param limit Integer to specify the number of returned samples/profiles. 
#' @param skip Integer to specify the number of skipped samples/profiles. If skip = 2, limit=500, 
#' the first 2*500 =1000 profiles are skipped and the next 500 profiles are returned. 
#' @param save_file A logical value determining whether to save the variant data as file 
#' instead of direct return. Only used when the parameter `type` is "variant" and `output` is "pgxseg" or "seg". Default is FALSE.
#' @param filename A string specifying the path and name of the file to be saved. 
#' Only used if the parameter `save_file` is TRUE. Default is "variants.seg/pgxseg" 
#' in current work directory.
#' @return Data from Progenetix database
#' @export

pgxLoader <- function(
    type = NULL,
    output  = NULL, 
    filters= NULL,
    codematches = FALSE, 
    limit=NULL,
    skip=NULL,
    biosample_id = NULL,
    save_file=FALSE,
    filename=NULL){
    
    if (!(any(type %in% c("biosample", "variant","frequency")))){
        stop("The parameter 'type' is invalid (available: \"biosample\", \"variant\", or \"frequency\")")
    }
    

    if (!is.null(output)){
      if (type == 'variant'){
        if (output=='pgxmatrix'){
          type = 'callset'
        }  else if (output=='coverage'){
          type= 'coverage'
        } else if (!output %in% c('pgxseg','seg')){
          stop("The parameter 'output' is invalid (available: NULL, \"pgxseg\", \"pgxmatrix\", \"seg\", or \"coverage\")")
        }
      }}
     
  
  
      
    if (type %in% c('callset', 'coverage','frequency')){
      if (!is.null(biosample_id)){
        cat("\n WARNING: The parameter 'biosample_id' is not used in this query. Only 'filters' are accepted. \n")
      }
    }
      

  
    if (type == "variant"){
      if (is.null(biosample_id)){
        stop("The parameter biosample_id cannot be NULL")
      }
      if (!is.null(filters)){
        cat("\n WARNING: The parameter 'filters' is not used in this query. Only 'biosample_id' is accepted. \n")
      }
    }
  
  if (type == "biosample" & is.null(limit)){
    limit=0
  }
  
    switch(type,
           biosample = pgxSampleLoader(biosample_id = biosample_id, filters = filters,codematches = codematches,skip=skip,limit=limit),
           variant= pgxVariantLoader(biosample_id = biosample_id,output=output,save_file=save_file, filename = filename),
           frequency =pgxFreqLoader(output = output, codematches = codematches, filters=filters),
           callset=pgxcallsetLoader(filters = filters,limit=limit, skip=skip,codematches = codematches),
           coverage = pgxCovLoader(filters = filters, codematches = codematches,skip=skip,limit=limit))
} 

