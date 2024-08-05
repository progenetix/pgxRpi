#' Plot survival data of individuals
#'
#' This function provides the survival plot from individual metadata. 
#'
#' @param data The object returned by the `pgxLoader` function, which includes survival data about individuals.
#' @param group_id A string specifying which column is used for grouping in the Kaplan-Meier plot.
#' @param condition A string for splitting individuals into younger and older groups, following the ISO 8601 duration format. Only used if `group_id` is "age_iso".
#' @param return_data A logical value determining whether to return the metadata used for plotting. Default is FALSE.
#' @param ... Other parameters relevant to KM plot. These include `pval`, `pval.coord`, `pval.method`, `conf.int`, `linetype`, and `palette` (see ggsurvplot from survminer)
#' @return The KM plot from input data
#' @export
#' @examples
#' individuals <- pgxLoader(type="individuals",filters="NCIT:C3512")
#' pgxMetaplot(individuals, group_id="age_iso", condition="P65Y")
pgxMetaplot <- function(data, group_id, condition, return_data = FALSE,...){
    # remove samples without survival data
    meta.sel.df <- data[!is.na(data$followup_time) & !data$followup_state_id %in% c('EFO:0030039','',NA),]
    if(nrow(meta.sel.df) == 0) stop("\n No available survival data \n")
    
    # transform ISOdate interal to days
    plot_followup_time <- lubridate::time_length(meta.sel.df$followup_time,unit='day')
  
    # process group id
    ## for age related group id
    if (group_id == "age_iso"){
        plot_group_time <- meta.sel.df[,group_id]
        time_condition <- condition
        ### transform to days for comparison     
        plot_group_time <- lubridate::time_length(plot_group_time,unit='day')
        time_condition <- lubridate::time_length(time_condition,unit='day')
        if (is.numeric(condition) | is.na(time_condition)) stop("\n `condition` argument is invalid \n")
        
        ### split group
        plot_group_id <- rep(NA,length(plot_group_time))
        plot_group_id[plot_group_time > time_condition] <- paste0("Age>",condition)
        plot_group_id[plot_group_time <= time_condition] <- paste0("Age<=",condition)
    ## for other group ids
    } else{
        plot_group_id <- meta.sel.df[,group_id]
    }
  
  
    # Map EFO code to numerical representation for plotting
    plot_followup_state_id <- meta.sel.df$followup_state_id
    plot_followup_state_id[plot_followup_state_id == "EFO:0030041"] <- 0 # alive
    plot_followup_state_id[plot_followup_state_id == "EFO:0030049"] <- 1 # death
    plot_followup_state_id <- as.numeric(plot_followup_state_id)
    # fit 
    plot_data <- data.frame(followup_time=plot_followup_time, followup_state=plot_followup_state_id, group=plot_group_id)
    sfit <- survival::survfit(survival::Surv(followup_time,followup_state)~group,data=plot_data)
  
    splot.op <- list(
      fit=sfit,
      data=plot_data,
      palette=NULL,
      linetype=1,
      conf.int=FALSE,
      pval = FALSE,
      pval.method=FALSE,
      pval.coord = NULL,
      ggtheme=ggplot2::theme_light()
    )
    splot.op <- modifyList(splot.op,list(...))
  
    ggp <- do.call(survminer::ggsurvplot,splot.op)
    show(ggp$plot)
    if (return_data) return(meta.sel.df)
}
