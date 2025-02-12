#' Load data from Progenetix database via the Beacon v2 API with some extensions
#'
#' This function loads various data from `Progenetix` database via the Beacon v2 API with some extensions (BeaconPlus).   
#'
#' @param type A string specifying the type of output data. Available options include:
#'   - `"individuals"`: Returns information about individuals.
#'   - `"biosamples"`: Returns information about biosamples.
#'   - `"analyses"`: Returns information about analyses.
#'   - `"g_variants"`: Returns variants data.
#'   - `"filtering_terms"`: Returns all available filtering terms.
#'   - `"counts"`: Returns the count of results based on the specified filters.
#'   - `"cnv_frequency"`: Returns precomputed CNV frequency data from Progenetix.
#'   - `"cnv_fraction"`: Returns CNV fraction per sample based on Progenetix data.
#' @param output A string specifying the format of the output data. The available options depend on the value of the `type` parameter:
#'   - If `type` is `"g_variants"`, the available options are `NULL` (default), `"pgxseg"`, or `"seg"`.
#'   - If `type` is `"cnv_frequency"`, the available options are `"pgxfreq"` (default) or `"pgxmatrix"`.
#'   - If `type` is `"cnv_fraction"`, the available options are `NULL` (default) or `"pgxmatrix"`.
#' @param biosample_id Identifiers used in the query database for identifying biosamples. 
#' @param individual_id  Identifiers used in the query database for identifying individuals. 
#' @param filters Identifiers used in public repositories, bio-ontology terms, or custom terms such as `c("NCIT:C7376", "PMID:22824167")`. 
#' When multiple filters are used, they are combined using AND logic when the parameter `type` is `"individuals"`, `"biosamples"`, or `"analyses"`; OR logic when the parameter `type` is `"counts"` or `"cnv_frequency"`.
#' @param limit Integer to specify the number of returned profiles. Default is `0` (return all). 
#' @param skip An integer specifying the number of profiles to skip. For example, if `skip = 2` and `limit = 500`, the first `2 * 500 = 1000` profiles are skipped, 
#' and the next 500 profiles are returned. Default is `0`, meaning no profiles are skipped.
#' @param dataset Datasets to query from the Beacon response. Default is `NULL`, which includes results from all datasets.
#' @param codematches A logical value indicating whether to exclude samples from child concepts of the specified filters in the ontology tree. 
#' If `TRUE`, only samples that exactly match the specified filters will be included. This parameter should not be used when `filters` include ontology-irrelevant filters, such as PMID or cohort identifiers. 
#' Default is `FALSE`. This option is applicable only when querying data resources are Progenetix or cancercelllines.org.
#' @param save_file A logical value determining whether to save variant data as a local file instead of direct return. Only used when the parameter `type` is `"g_variants"`. Default is `FALSE`.
#' @param filename A string specifying the path and name of the file to be saved. This parameter is used only when `save_file` is set to `TRUE`. The default value is `"variants.tsv"`, saved in the current working directory.
#' @param domain The domain of the query data resource. Default is `"http://progenetix.org"`.
#' @param entry_point The entry point of the Beacon v2 API. Default is `"beacon"`, resulting in the default endpoint being "http://progenetix.org/beacon"
#' @param num_cores An integer specifying the number of cores to use for parallel processing during Beacon v2 phenotypic/meta-data queries from multiple domains or variant data queries from multiple biosamples. Default is `1`.
#' @importFrom utils URLencode modifyList read.table write.table
#' @importFrom httr GET content
#' @return Data from Progenetix database
#' @export
#' @examples
#' ## query metadata
#' biosamples <- pgxLoader(type="biosamples", filters = "NCIT:C3512")
#' ## query variants
#' seg <- pgxLoader(type="g_variants", biosample_id = "pgxbs-kftvgx4y")
#' ## query CNV frequency
#' freq <- pgxLoader(type="cnv_frequency", output ='pgxfreq', filters="NCIT:C3512")

pgxLoader <- function(
    type=NULL,
    output=NULL, 
    biosample_id = NULL,
    individual_id=NULL,
    filters= NULL,
    limit=0,
    skip=0,
    dataset=NULL,
    codematches = FALSE, 
    save_file=FALSE,
    filename="variants.tsv",
    domain="http://progenetix.org",
    entry_point="beacon",
    num_cores=1){
    
    type <- match.arg(type, c("biosamples", "individuals","g_variants","analyses","filtering_terms","cnv_frequency","cnv_fraction","counts"))
       
    # specify output 
    if (is.null(output) & type %in% c("g_variants","cnv_fraction")){
        output <-  switch(type,
                          g_variants=NULL,
                          cnv_fraction=NULL)                     
    } else{
        output <-  switch(type,
                          g_variants=match.arg(output, c("pgxseg", "seg")),
                          cnv_frequency=match.arg(output, c("pgxfreq" , "pgxmatrix")),
                          cnv_fraction=match.arg(output, "pgxmatrix"))

    }
      
    # parameter usage warnings     
    if (type %in% c("cnv_frequency","counts")){
        check_missing_parameters(filters,"'filters'")
        check_unused_parameters(biosample_id, "'biosample_id'", "'filters'")
        check_unused_parameters(individual_id, "'individual_id'", "'filters'")      
    }

    if (type=="g_variants"){
        if (any(domain %in% c("http://progenetix.org","progenetix.org","https://cancercelllines.org","cancercelllines.org"))){
            check_missing_parameters(biosample_id,"'biosample_id'")
        } 
        check_unused_parameters(individual_id, "'individual_id'", "'biosample_id'")
        check_unused_parameters(filters, "'filters'", "'biosample_id'")
    }

    if (type == "cnv_fraction"){
        if (length(filters) > 1) stop("This query only supports one filter")
    }

    if (type %in% c("analyses","filtering_terms","g_variants","cnv_frequency","counts")){
        if (codematches) warning("\n The parameter 'codematches' is not used in this query. \n")
    }


    # CNV fraction data are accessed by different endpoints
    if (type == "cnv_fraction" & !is.null(output)) type <- "samplematrix"

    options(timeout=500)
    switch(type,
           biosamples = pgxmetaLoader(type=type,biosample_id=biosample_id,individual_id=individual_id,filters=filters,codematches=codematches,skip=skip,limit=limit,domain=domain,entry_point=entry_point,dataset=dataset,num_cores=num_cores),
           individuals= pgxmetaLoader(type=type,biosample_id=biosample_id,individual_id=individual_id,filters=filters,codematches=codematches,skip=skip,limit=limit,domain=domain,entry_point=entry_point,dataset=dataset,num_cores=num_cores),
           analyses   = pgxmetaLoader(type=type,biosample_id=biosample_id,individual_id=individual_id,filters=filters,codematches=codematches,skip=skip,limit=limit,domain=domain,entry_point=entry_point,dataset=dataset,num_cores=num_cores),
           filtering_terms = pgxmetaLoader(type=type,biosample_id=NULL,individual_id=NULL,filters=NULL,codematches=FALSE,skip=NULL,limit=NULL,domain=domain,entry_point=entry_point,dataset=NULL,num_cores=num_cores),
           counts = pgxCount(filters,domain,entry_point,num_cores=num_cores),
           g_variants = pgxVariantLoader(biosample_id=biosample_id,output=output,save_file=save_file,filename=filename,domain=domain,entry_point=entry_point,dataset=dataset,num_cores=num_cores),
           cnv_frequency = pgxFreqLoader(output=output,filters=filters,domain=domain),
           samplematrix = pgxcallsetLoader(biosample_id=biosample_id,individual_id=individual_id,filters=filters,limit=limit,skip=skip,codematches=codematches,domain=domain),
           cnv_fraction = pgxFracLoader(biosample_id=biosample_id,individual_id=individual_id,filters=filters,codematches=codematches,skip=skip,limit=limit,domain=domain))     
} 

