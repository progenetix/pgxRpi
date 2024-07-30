#' Load data from Progenetix database via the Beacon v2 API with some extensions
#'
#' This function loads various data from `Progenetix` database via the Beacon v2 API with some extensions (BeaconPlus).   
#'
#' @param type A string specifying output data type. Available options are "biosamples", "individuals", "analyses",
#' "g_variants" or "cnv_frequency". The first three options return corresponding metadata, "g_variants" returns variants data with the focus on CNV, 
#' and "cnv_frequency" returns precomputed CNV frequency based on data in Progenetix. 
#' @param output A string specifying output data format. When the parameter `type` is "g_variants",
#' available options are NULL, "pgxseg", "seg", "cnvfraction", or "pgxmatrix"; When the parameter `type` is "cnv_frequency", available options are "pgxfreq" or "pgxmatrix".
#' @param biosample_id Identifiers used in the query database for identifying biosamples. 
#' @param individual_id  Identifiers used in the query database for identifying individuals. 
#' @param filters Identifiers used in public repositories, bio-ontology terms, or custom terms such as c("NCIT:C7376", "PMID<!-- -->:22824167"). 
#' When multiple filters are used, they are combined using AND logic when the parameter `type` is "biosamples", "individuals", or "analyses"; OR logic when the parameter `type` is "cnv_frequency".
#' @param limit Integer to specify the number of returned profiles. Default is 0 (return all). 
#' @param skip Integer to specify the number of skipped profiles. E.g. if skip = 2, limit=500, the first 2*500 =1000 profiles are skipped and the next 500 profiles are returned. Default is NULL (no skip).
#' @param codematches A logical value determining whether to exclude samples from child concepts of specified filters in the ontology tree. 
#' If TRUE, only samples exactly matching the specified filters will be included.  Do not use this parameter when `filters` include ontology-irrelevant filters such as PMID and cohort identifiers.
#' Default is FALSE.
#' @param save_file A logical value determining whether to save variant data as a local file instead of direct return. Only used when the parameter `type` is "g_variants". Default is FALSE.
#' @param filename A string specifying the path and name of the file to be saved. Only used if the parameter `save_file` is TRUE. Default is "variants" in current work directory.
#' @param num_cores Integer to specify the number of cores used for the variant query. Only used when the parameter `type` is "g_variants". Default is 1.
#' @param dataset A string specifying the dataset to query. When the parameter `domain` is "http://progenetix.org", available options are "progenetix" (by defualt) and "cancercelllines".
#' @param domain A string specifying the domain of the query data resource. Default is "http://progenetix.org".
#' @param entry_point A string specifying the entry point of the Beacon v2 API. Default is "beacon", resulting in the endpoint being "http://progenetix.org/beacon".
#' @importFrom utils URLencode modifyList read.table write.table
#' @importFrom httr GET content
#' @return Data from Progenetix database
#' @export
#' @examples
#' ## query metadata
#' biosamples <- pgxLoader(type="biosamples", filters = "NCIT:C3512")
#' ## query segment variants
#' seg <- pgxLoader(type="g_variants", output = "pgxseg", biosample_id = "pgxbs-kftvgx4y")
#' ## query CNV frequency
#' freq <- pgxLoader(type="cnv_frequency", output ='pgxfreq', filters="NCIT:C3512")

pgxLoader <- function(
    type=NULL,
    output=NULL, 
    biosample_id = NULL,
    individual_id=NULL,
    filters= NULL,
    limit=0,
    skip=NULL,
    codematches = FALSE, 
    save_file=FALSE,
    filename="variant",
    num_cores=1,
    dataset=NULL,
    domain="http://progenetix.org",
    entry_point="beacon"){
    
    type <- match.arg(type, c("biosamples", "individuals","g_variants","cnv_frequency","analyses"))
       
    # specify output 
    if (is.null(output) & type %in% c("g_variants","cnv_frequency")){
        output <-  switch(type,
                          g_variants=NULL,
                          cnv_frequency=match.arg(output, c("pgxfreq" , "pgxmatrix")))
                          
    } else{
        output <-  switch(type,
                          g_variants=match.arg(output, c("pgxseg", "seg", "pgxmatrix", "cnvfraction")),
                          cnv_frequency=match.arg(output, c("pgxfreq" , "pgxmatrix")),
                          analyses=NULL,
                          biosamples=NULL,
                          individuals=NULL)

    # adapt to different variant formats 
        if (type == "g_variants"){
          type <- switch(output,
                         pgxmatrix = 'callset',
                         cnvfraction = 'cnvstats',
                         seg="g_variants",
                         pgxseg="g_variants")
        }
    }

    if (type %in% c("callset", "cnvstats")){
        if (length(filters) > 1) stop("This query only supports one filter")
    }
    
    if (type == "cnv_frequency"){
      check_missing_parameters(filters,"'filters'")
      check_unused_parameters(biosample_id, "'biosample_id'", "'filters'")
      check_unused_parameters(individual_id, "'individual_id'", "'filters'")      
    }

    
    if (type %in% c("individuals","biosamples","analyses","callset", "cnvstats")){
        check_missing_parameters(c(filters,individual_id,biosample_id), c("'filters'","'individual_id'","'biosample_id'"))
    }
  
    if (type == "g_variants"){
        check_missing_parameters(biosample_id,"'biosample_id'")
        check_unused_parameters(filters, "'filters'", "'biosample_id'")
        check_unused_parameters(individual_id, "'individual_id'", "'biosample_id'")
    }

    if (type %in% c("analyses","cnv_frequency")){
        if (codematches) warning("\n The parameter 'codematches' is not used in this query. \n")
    }

    
    options(timeout=500)
    switch(type,
           biosamples = pgxmetaLoader(type=type,biosample_id=biosample_id,individual_id=individual_id,filters=filters,codematches=codematches,skip=skip,limit=limit,domain=domain,entry_point=entry_point,dataset=dataset),
           individuals= pgxmetaLoader(type=type,biosample_id=biosample_id,individual_id=individual_id,filters=filters,codematches=codematches,skip=skip,limit=limit,domain=domain,entry_point=entry_point,dataset=dataset),
           analyses   = pgxmetaLoader(type=type,biosample_id=biosample_id,individual_id=individual_id,filters=filters,codematches=codematches,skip=skip,limit=limit,domain=domain,entry_point=entry_point,dataset=dataset),
           g_variants = pgxVariantLoader(biosample_id=biosample_id,output=output,save_file=save_file,filename=filename,domain=domain,entry_point=entry_point,dataset=dataset,num_cores=num_cores),
           cnv_frequency = pgxFreqLoader(output=output,filters=filters,domain=domain,dataset=dataset),
           callset = pgxcallsetLoader(biosample_id=biosample_id,individual_id=individual_id,filters=filters,limit=limit,skip=skip,codematches=codematches,domain=domain,dataset=dataset),
           cnvstats = pgxFracLoader(biosample_id=biosample_id,individual_id=individual_id,filters=filters,codematches=codematches,skip=skip,limit=limit,domain=domain,dataset=dataset))
} 

