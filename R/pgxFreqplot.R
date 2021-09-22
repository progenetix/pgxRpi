#' Plot CNV frequency data
#'
#' Thie function plots the frequency of deletions and amplifications
#'
#' @param data CNV frequency data returned by `pgxLoader` function. The frequency matrices
#' in `data` slot must be stored as `.pgxseg` format, which can be specified by `output` 
#' parameter of `pgxLoader` function. 
#' @param chrom A vector with chromosomes to be plotted. If NULL, return the plot
#' by genome. If specified the frequencies are plotted with one panel for each
#' chromosome. Default is NULL.
#' @param layout Number of columns and rows in plot. Only used in plot by chromosome.
#' Default is c(1,1).
#' @param filters String value to indicate which filter to be plotted, such as 1
#' (the first filters in `data` slot of object ) or 'NCIT:C4038' (specific filter name). The length of filters
#' is limited to one if the parameter `circos` is False. Default is 1.
#' @param circos A logical value to indicate if return a circos plot. If TRUE, it
#' can return a circos plot with multiple filters for display and comparison.
#' Default is FALSE.
#' @return The binned CNV frequency plot
#' @export


pgxFreqplot <- function(data,chrom=NULL,layout=c(1,1),filters=NULL,circos = FALSE){
    if (is.null(filters)){
        filters <- names(data$data)[1]
    }

    if (all(!(filters %in% names(data$data)))){
        stop(paste("The filter is not contained in data:", filters[!(filters %in% names(data$data))]))
    }

    if (circos){
        return(cplotpgxFreq(data$data[filters]))
    }

    if (length(filters) > 1){
        stop("The length of filters exceeds limit")
    }

    pos.unit="bp"
    type <- ifelse(is.null(chrom),"genome","bychrom")
    switch(type,
            genome = genomeFreq(data,pos.unit,filters),
            bychrom = chromosomeFreq(data,pos.unit,chrom,layout,filters))
}

genomeFreq <- function(data,pos.unit,id){
    data_pgxseg <- data$data[[id]]
    nr <- 1
    nc <- 1
    op <- getFreqPlotParameters(type="genome",nc=nc,nr=nr)
    data_pgxseg[,2] <- gsub("X", "23", data_pgxseg[,2])
    data_pgxseg[,2] <- gsub("Y", "24", data_pgxseg[,2])
    op$xlim <- getGlobal.xlim(op=op,pos.unit=pos.unit,chrom=as.numeric(unique(data_pgxseg[,2])))
    x <- adjustPos(position=data_pgxseg[,4],chromosomes=as.numeric(data_pgxseg[,2]),pos.unit=pos.unit,type="genome",op=op)
    xleft <- x$xleft
    xright <- x$xright
    if(dev.cur()<=1){       #to make Sweave work
        dev.new(width=op$plot.size[1],height=op$plot.size[2],record=TRUE)
    }


    #Initialize:
    row=1
    clm=1
    new = FALSE

    #Division of plotting window:
    frames <- framedim(1,1)

    fig <- c(frames$left[clm],frames$right[clm],frames$bot[row],frames$top[row])

    par(fig=fig,new=new,oma=c(0,0,1,0),mar=op$mar)

    op <- updateFreqParameters(data_pgxseg$loss_frequency,data_pgxseg$gain_frequency,op)
    plot(1,1,type="n",ylim=op$ylim,xlim=op$xlim,xaxs="i",main="",frame.plot=TRUE,yaxt="n",xaxt="n",ylab="",xlab="")
    chromPattern(pos.unit,op)

    id_name <- names(data$data[id])
    meta <- data$meta[data$meta$code == id_name,]
    op$main <- paste(id_name, ": ",meta[2], " (", meta[3], " samples)", sep="")


    title(main=op$main,line=op$main.line,cex.main=op$cex.main)

    #Add axes, labels and percentage lines:
    addToFreqPlot(op,type="genome")
    rect(xleft=xleft,ybottom=0,xright=xright,ytop=data_pgxseg$gain_frequency,col=op$col.gain,border=op$col.gain)
    rect(xleft=xleft,ybottom=0,xright=xright,ytop=-data_pgxseg$loss_frequency,col=op$col.loss,border=op$col.loss)
    abline(h=0,lty=1,col="grey82",lwd=1.5)

    op$chrom.lty = 1
    addChromlines(as.numeric(data_pgxseg[,2]),xaxis="pos",unit=pos.unit,cex=op$cex.chrom,op=op)
    addArmlines(as.numeric(data_pgxseg[,2]),xaxis="pos",unit=pos.unit,cex=op$cex.chrom,op=op)
}

chromosomeFreq <- function(data,pos.unit,chrom,layout,id){
    data_pgxseg <- data$data[[id]]
    nProbe <- nrow(data_pgxseg)
    nChrom <- length(chrom)

    #Grid layout
    nr <- layout[1]
    nc <- layout[2]

    #Get plot parameters:
    op <- getFreqPlotParameters(type="bychrom",nc=nc,nr=nr,chrom=chrom)

    #Margins for entire plot in window:
    if(all(op$title=="")){
        oma <- c(0,0,0,0)
    }else{
        oma <- c(0,0,1,0)
    }
    mar=c(0.2,0.2,0.3,0.2)

    data_pgxseg[,2] <- gsub("X", "23", data_pgxseg[,2])
    data_pgxseg[,2] <- gsub("Y", "24", data_pgxseg[,2])
    #Adjust positions to be plotted along xaxis; i.e. scale according to plot.unit, and get left and right pos for freq.rectangles to be plotted (either continuous
    #or 1 probe long):
    x <- adjustPos(position=data_pgxseg[,4],chromosomes=as.numeric(data_pgxseg[,2]),pos.unit=pos.unit,type="chromosome",op=op)
    xleft <- x$xleft
    xright <- x$xright

    #Divide the plotting window by the function "framedim":
    frames <- framedim(nr,nc)

    #make separate plots for each value of thres.gain/thres.loss


    if(dev.cur()<=1){       #to make Sweave work
        dev.new(width=op$plot.size[1],height=op$plot.size[2],record=TRUE)
    }

    #Initialize row and column index:
    row=1
    clm=1
    new = FALSE


    #Find default ylimits and at.y (tickmarks):
    use <- which(as.numeric(data_pgxseg[,2]) %in% chrom)
    op <- updateFreqParameters(data_pgxseg$loss_frequency[use],data_pgxseg$gain_frequency[use],op)

    #Make separate plots for each chromosome:

    for(c in 1:nChrom){

        #Frame dimensions for plot c:
        fig.c <- c(frames$left[clm],frames$right[clm],frames$bot[row],frames$top[row])
        par(fig=fig.c,new=new,oma=oma,mar=mar)

        #Make list with frame dimensions:
        frame.c <- list(left=frames$left[clm],right=frames$right[clm],bot=frames$bot[row],top=frames$top[row])

        #Select relevant chromosome number
        k <- chrom[c]

        #Pick out frequencies for this chromsome
        ind.c <- which(as.numeric(data_pgxseg[,2]) ==k)
        freqamp.c <- data_pgxseg$gain_frequency[ind.c]
        freqdel.c <- data_pgxseg$loss_frequency[ind.c]

        xlim <- c(0,max(xright[ind.c]))

        #Plot ideogram at below frequencies:
        if(op$plot.ideo){
            #Ideogram frame:
            ideo.frame <- frame.c
            ideo.frame$top <- frame.c$bot + (frame.c$top-frame.c$bot)*op$ideo.frac

            par(fig=unlist(ideo.frame),new=new,mar=op$mar.i)
            #Plot ideogram and get maximum probe position in ideogram:
            plotIdeogram(chrom=k,cyto.text=op$cyto.text,cyto.data=op$assembly,cex=op$cex.cytotext,unit=op$plot.unit)

            #Get maximum position for this chromosome:
            xmaxI <- chromMax(chrom=k,cyto.data=op$assembly,pos.unit=op$plot.unit)
            xlim <- c(0,xmaxI)

            new <- TRUE
        }

        #Freq.plot-dimensions:
        frame.c$bot <- frame.c$bot + (frame.c$top-frame.c$bot)*op$ideo.frac
        par(fig=unlist(frame.c),new=new,mar=op$mar)

        #Limits:
        if(!is.null(op$xlim)){
            xlim <- op$xlim
        }

        #Empty plot:
        plot(1,1,type="n",ylim=op$ylim,xlim=xlim,xaxs="i",main=op$main[c],frame.plot=FALSE,yaxt="n",xaxt="n",cex.main=op$cex.main,ylab="",xlab="")

        #Add axes, labels and percentageLines:
        chrom.op <- op
        chrom.op$xlim <- xlim

        addToFreqPlot(chrom.op,type="bychrom")

        #Plot frequencies as rectangles
        rect(xleft=xleft[ind.c],ybottom=0,xright=xright[ind.c],ytop=freqamp.c,col=op$col.gain,border=op$col.gain)
        rect(xleft=xleft[ind.c],ybottom=0,xright=xright[ind.c],ytop=-freqdel.c,col=op$col.loss,border=op$col.loss)


        #Add line at y=0 and x=0
        abline(h=0,lty=1,col="grey90")
        abline(v=0)


        id_name <- names(data$data[id])
        meta <- data$meta[data$meta$code == id_name,]
        op$title <- paste(id_name, ": ",meta[2], " (", meta[3], " samples)", sep="")


        #If page is full; start plotting on new page
        if(c%%(nr*nc)==0 && c!=nChrom){
            #Add main title to page:
            title(op$title,outer=TRUE,line=-0.6,cex.main=0.85)

            #Start new page when prompted by user:
            devAskNewPage(ask = TRUE)

            #Reset columns and row in layout:
            clm = 1
            row = 1
            new=FALSE

        }else{
            #Update column and row index:
            if(clm<nc){
                clm <- clm+1
            }else{
                clm <- 1
                row <- row+1
            }#endif
            new=TRUE
        }#endif

    }#endfor

    title(op$title,outer=TRUE,line=-0.6,cex.main=0.85)
}

cplotpgxFreq  <- function(data_pgxseg){
    cex_text <- ifelse(length(data_pgxseg) == 1, 0.8, 0.6)
    id_name <-  names(data_pgxseg)
    id_1 <- id_name[1]
    data_cir_1 <- data_pgxseg[[id_1]][,c(2,3,4,5,6)]
    data_cir_1[,5] <- -data_cir_1[,5]
    data_cir_1[,1] <- paste0('chr',data_cir_1$reference_name)
    circlize::circos.clear()
    circlize::circos.par(start.degree=90, gap.degree=c(rep(1,23),10))
    circlize::circos.initializeWithIdeogram(species = "hg38")
    circlize::circos.genomicTrack(data_cir_1, ylim = c(-100, 100),numeric.column = c("gain_frequency", "loss_frequency"),
                        panel.fun = function(region, value, ...) {
                            value1 = unlist(value[1])
                            value2 = unlist(value[2])
                            pos=unlist(region[1]+(region[2]-region[1])/2)
                            circlize::circos.barplot(value=value1,pos=pos, border = 'red',col = 'red')
                            circlize::circos.barplot(value=value2,pos=pos, border = 'blue',col = 'blue')
                            y = seq(-100,100,by=20)
                            circlize::circos.segments(0,y,circlize::CELL_META$xlim[2], y,col=ifelse(y>0,'#FF9999','#99CCFF'))
                        })
    text(0,0,id_1,cex = cex_text)
    if (length(id_name) > 1){
        for (i in c(2:length(id_name))){
            id <- id_name[i]
            data_cir <- data_pgxseg[[id]][,c(2,3,4,5,6)]
            data_cir[,5] <- -data_cir[,5]
            data_cir[,1] <- paste0('chr',data_cir$reference_name)
            circlize::circos.genomicTrack(data_cir, ylim = c(-100, 100),numeric.column = c("gain_frequency", "loss_frequency"),
                                          panel.fun = function(region, value, ...) {
                                              value1 = unlist(value[1])
                                              value2 = unlist(value[2])
                                              pos=unlist(region[1]+(region[2]-region[1])/2)
                                              circlize::circos.barplot(value=value1,pos=pos, border = 'red',col = 'red')
                                              circlize::circos.barplot(value=value2,pos=pos, border = 'blue',col = 'blue')
                                              y = seq(-100,100,by=20)
                                              circlize::circos.segments(0,y,circlize::CELL_META$xlim[2], y,col=ifelse(y>0,'#FF9999','#99CCFF'))
                                          })
            text(0,-0.1*(i-1),id,cex = cex_text)


        }
         }
}
