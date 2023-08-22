#' Load data from Progenetix database 
#'
#' This function loads various data from `Progenetix` database.   
#'
#' @param type A string specifying output data type. Available options are "biosample", "individual",
#' "variant" or "frequency". The first two options return corresponding metadata, "variant" returns CNV variant data, and "frequency" returns precomputed CNV frequency based on data in Progenetix. 
#' @param output A string specifying output data format. When the parameter `type` is "variant",
#' available options are NULL, "pgxseg", "seg", "coverage", or "pgxmatrix"; When the parameter `type` is "frequency",
#' available options are "pgxfreq" or "pgxmatrix".
#' @param filters Identifiers for cancer type, literature, cohorts, and age such as c("NCIT:C7376", "pgx:icdom-98353", "PMID:22824167", "pgx:cohort-TCGAcancers", "age:>=P50Y"). 
#' @param codematches A logical value determining whether to exclude samples from child concepts of specified filters that belong to cancer type/tissue encoding system (NCIt, icdom/t, Uberon). 
#' If TRUE, retrieved samples only keep samples exactly encoded by specified filters. 
#' Do not use this parameter when `filters` include cancer-irrelevant filters such as PMID and cohort identifiers.
#' Default is FALSE.
#' @param filterLogic A string specifying logic for combining multiple filters when query metadata (the paramter `type` = "biosample" or "individual"). Available options are "AND" and "OR". Default is "AND".  An exception is filters associated with age that always use AND logic 
#' when combined with any other filter, even if filterLogic = "OR", which affects other filters.  Note that when `type` = "frequency", the combining logic is "OR", which is not changed by this parameter.
#' @param limit Integer to specify the number of returned biosample/individual/variant profiles for each filter. Default is 0 (return all). 
#' @param skip Integer to specify the number of skipped biosample/individual/variant profiles for each filter. E.g. if skip = 2, limit=500, 
#' the first 2*500 =1000 profiles are skipped and the next 500 profiles are returned. Default is NULL (no skip).
#' @param biosample_id Identifiers used in Progenetix database for identifying biosamples. 
#' @param individual_id  Identifiers used in Progenetix database for identifying individuals. 
#' @param save_file A logical value determining whether to save the segment variant data as file 
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
    filterLogic = "AND",
    limit=0,
    skip=NULL,
    biosample_id = NULL,
    individual_id=NULL,
    save_file=FALSE,
    filename=NULL){
    
    if (!(any(type %in% c("biosample", "variant","frequency","individual")))){
        stop("\n The parameter 'type' is invalid (available: \"biosample\", \"individual\", \"variant\", or \"frequency\")")
    }
    

    if (!is.null(output)){
      if (type == 'variant'){
        if (output=='pgxmatrix'){
          type = 'callset'
        }  else if (output=='coverage'){
          type= 'coverage'
        } else if (!output %in% c('pgxseg','seg')){
          stop("\n The parameter 'output' is invalid (available: NULL, \"pgxseg\", \"pgxmatrix\", \"seg\", or \"coverage\")")
        }
      }
      if (type == 'frequency'){
        if (!output %in% c('pgxfreq','pgxmatrix')){
          stop("\n The parameter 'output' is invalid (available: \"pgxfreq\" or \"pgxmatrix\")")
        }
      }
    }    
  
      
    if (type %in% c('callset', 'coverage','frequency')){
      if (!is.null(biosample_id)){
        cat("\n WARNING: The parameter 'biosample_id' is not used in this query. Only 'filters' are accepted. \n")
      }
      if (!is.null(individual_id)){
        cat("\n WARNING: The parameter 'individual_id' is not used in this query. Only 'filters' are accepted. \n")
      }
      if (length(filters) < 1){
        stop("\n The parameter 'filters' is missing. At least one valid filter has to be provided")
      }

      if (type %in% c('callset', 'coverage')){
        if (length(filters) > 1) stop("\n The parameter 'filters' is invalid. This query only supports one filter")
      }

      if (type == "frequency"){
        if (is.null(output)) stop("\n The parameter 'output' is missing.  It must be set as 'pgxfreq' or 'pgxmatrix'")  
      }
    }
      
    if (type == "individual"){
      if (!is.null(biosample_id)){
        cat("\n WARNING: The parameter 'biosample_id' is not used in this query. Only 'individual_id' or 'filters' are accepted. \n")
      }
    }
  
    if (type == "variant"){
      if (is.null(biosample_id)){
        stop("\n The parameter 'biosample_id' is missing. At least one valid biosample id has to be provided")
      }
      if (!is.null(filters)){
        cat("\n WARNING: The parameter 'filters' is not used in this query. Only 'biosample_id' is accepted. \n")
      }
      if (!is.null(individual_id)){
        cat("\n WARNING: The parameter 'individual_id' is not used in this query. Only 'biosample_id' is accepted. \n")
      }
    }

 
  
    switch(type,
           biosample = pgxmetaLoader(type=type,biosample_id= biosample_id,individual_id=individual_id,filters=filters,codematches=codematches,skip=skip,limit=limit,filterLogic=filterLogic),
           individual= pgxmetaLoader(type=type,biosample_id= biosample_id,individual_id=individual_id,filters=filters,codematches=codematches,skip=skip,limit=limit,filterLogic=filterLogic),
           variant= pgxVariantLoader(biosample_id = biosample_id,output=output,save_file=save_file, filename = filename),
           frequency =pgxFreqLoader(output = output, codematches = codematches, filters=filters),
           callset=pgxcallsetLoader(filters = filters,limit=limit, skip=skip,codematches = codematches),
           coverage = pgxCovLoader(filters = filters, codematches = codematches,skip=skip,limit=limit))
} 

