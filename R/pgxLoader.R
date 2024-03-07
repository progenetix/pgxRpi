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
#' @param domain A string specifying the domain of database. Default is "http://progenetix.org".
#' @param dataset A string specifying the dataset to query. Default is "progenetix". Other available options are "cancercelllines".
#' @importFrom utils URLencode modifyList read.table write.table
#' @importFrom httr GET content
#' @return Data from Progenetix database
#' @export
#' @examples
#' ## query metadata
#' biosamples <- pgxLoader(type="biosample", filters = "NCIT:C3512")
#' ## query segment variants
#' seg <- pgxLoader(type="variant", output = "pgxseg", biosample_id = "pgxbs-kftvgx4y")
#' ## query CNV frequency
#' freq <- pgxLoader(type="frequency", output ='pgxfreq', filters="NCIT:C3512")

pgxLoader <- function(
    type=NULL,
    output=NULL, 
    filters= NULL,
    codematches = FALSE, 
    filterLogic = "AND",
    limit=0,
    skip=NULL,
    biosample_id = NULL,
    individual_id=NULL,
    save_file=FALSE,
    filename=NULL,
    domain="http://progenetix.org",
    dataset="progenetix"){
    
    type <- match.arg(type, c("biosample", "individual","variant","frequency"))
    dataset <- match.arg(dataset, c("progenetix","cancercelllines","cellz","examplez"))
    
    # actual dataset name for cancercelllines
    if (dataset == "cancercelllines") dataset <- "cellz"
    
    if (is.null(output) & type %in% c("variant","frequency")){
        output <-  switch(type,
                          variant=NULL,
                          frequency=match.arg(output, c("pgxfreq" , "pgxmatrix")))
                          
    } else{
        output <-  switch(type,
                          variant=match.arg(output, c("pgxseg", "seg", "pgxmatrix", "coverage")),
                          frequency=match.arg(output, c("pgxfreq" , "pgxmatrix")),
                          biosample=NULL,
                          individual=NULL)
        # adapt to different variant formats 
        if (type == "variant"){
          type <- switch(output,
                         pgxmatrix = 'callset',
                         coverage = 'coverage',
                         seg="variant",
                         pgxseg="variant")
        }
    }

    
    if (type %in% c('callset', 'coverage')){
        if (length(filters) > 1) stop("\n The parameter 'filters' is invalid. This query only supports one filter")
    }
    
    if (type == "frequency"){
      checkMissingParameters(filters,"'filters'")
      checkUnusedParameters(biosample_id, "'biosample_id'", "'filters'")
      checkUnusedParameters(individual_id, "'individual_id'", "'filters'")
    }
    
    if (type == "individual"){
        checkMissingParameters(c(filters,individual_id,biosample_id), c("'filters'","'individual_id'","'biosample_id'"))
    }
    
    if (type == "biosample"){
      checkMissingParameters(c(filters,individual_id,biosample_id), c("'filters'","'individual_id'","'biosample_id'"))
    }
  
    if (type == "variant"){
        checkMissingParameters(biosample_id,"'biosample_id'")
        checkUnusedParameters(filters, "'filters'", "'biosample_id'")
        checkUnusedParameters(individual_id, "'individual_id'", "'biosample_id'")
    }
    
    options(timeout=500)
    switch(type,
           biosample = pgxmetaLoader(type=type,biosample_id= biosample_id,individual_id=individual_id,filters=filters,codematches=codematches,skip=skip,limit=limit,filterLogic=filterLogic,domain=domain,dataset=dataset),
           individual= pgxmetaLoader(type=type,biosample_id= biosample_id,individual_id=individual_id,filters=filters,codematches=codematches,skip=skip,limit=limit,filterLogic=filterLogic,domain=domain,dataset=dataset),
           variant = pgxVariantLoader(biosample_id = biosample_id,output=output,save_file=save_file, filename = filename,domain=domain,dataset=dataset),
           frequency = pgxFreqLoader(output = output, codematches = codematches, filters=filters,domain=domain,dataset=dataset),
           callset = pgxcallsetLoader(biosample_id = biosample_id, individual_id=individual_id, filters = filters,limit=limit, skip=skip,codematches = codematches,domain=domain,dataset=dataset),
           coverage = pgxCovLoader(biosample_id = biosample_id, individual_id=individual_id, filters = filters, codematches = codematches,skip=skip,limit=limit,domain=domain,dataset=dataset))
} 

