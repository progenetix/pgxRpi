#Input:
### type: plot type (genome or bychrom)
### nc: number of columns in plot
### nr: number of rows in plot
### thres.gain,thres.loss: aberration calling thresholds
### chrom: a vector giving the chromosomes to be plotted, only used for type=bychrom
### ... : other optional plot parameters specified by user

#Output:
### op: a list containing default and user modified plot parameters

###Required by: 
### plotFreq (genomeFreq and chromosomeFreq)

getFreqPlotParameters <- function(type,nc,nr,chrom=NULL,...){
    
    #Apply a scaling factor according to number of columns and rows in plot:
    #seems to work ok:
    cr <- nc*nr
    f <- 1-0.013*cr
    
    
    #Common default parameters for genome and bychrom:
    op <- list(ylab="% with gain or loss",
               plot.size=c(11.8,min(3*nr,8.2)),
               col.gain="#FFC633",
               col.loss="#33A0FF",
               plot.unit="mbp",
               percentLines=TRUE,
               continuous=TRUE,
               assembly="hg38",
               las=1,
               chrom.lty=5,
               chrom.side=1,
               chrom.col="darkgrey",
               cyto.text=FALSE,
               #Parameters that will be set later
               ylim=NULL,
               xlim=NULL,
               xlab=NULL,
               mgp=NULL,
               mar=NULL,
               at.x=NULL,
               at.y=NULL,
               ideo.frac=NA,
               #Parameters that depend on the grid layout of the plot
               f=f,
               main.line=0.6*f,
               cex.lab=0.9*f,
               cex.main=f,
               cex.axis=0.8*f,
               cex.cytotext=0.6*f,
               cex.chrom=0.8*f
    )
    
    #Defaults specific to plot type:
    #For genome plot:
    if(type=="genome"){
        op$main <- paste("CNV frequency")
        op$main.line <- 1.5*f
        op$xlab <- ""
        op$plot.ideo <- FALSE
        op$mar <- c(1.5*op$f,3*op$f,2.3*op$f,1*op$f)
    }
    
    #For chromosome plot:
    if(type=="bychrom"){
        op$main <- paste("Chromosome ",chrom,sep="")
        op$plot.ideo=TRUE
    }
    
    
    
    #Set/modify parameters more depending on user input:
    #Check for user modifications
    op <- modifyList(op,list(...))
    
    #Set assembly to refer to stored data instead of character string:
    op$assembly <- get(op$assembly)
    
    #Placement of labels and axis annotation:
    if(is.null(op$mgp)){
        op$mgp <- c(1.3,0.05,0)*f
        mgp.y <- c(2,0.5,0)*f
    }else{
        mgp.y <- op$mgp
    }
    op$mgp.y <- mgp.y
    
    #Xlabel
    if(is.null(op$xlab)){
        op$xlab <- paste("Position (",op$plot.unit,")",sep="")
    }
    
    #margins:
    if(is.null(op$mar)){
        op$mar <- if(op$plot.ideo) c(0.2*f,3*f,2.5*f,f) else c(1.5*f,3*f,2.5*f,f)    
    }
    
    #Set default ideo.frac and ideogram margins:
    if(op$plot.ideo){
        #ideogram margins:
        op$mar.i <- c(0.2*f,3*f,0,f)
        if(op$cyto.text){
            #Need to increase bottom margin:
            op$mar.i <- op$mar.i + c(2,0,0,0)
        }
        #Make sure left and right margins are equal for mar and mar.i:
        op$mar.i[c(2,4)] <- op$mar[c(2,4)] 
        if(is.na(op$ideo.frac)){
            #ideo.frac has not been defined by user:
            op$ideo.frac <- 0.05*sqrt(sqrt(cr))
            if(op$cyto.text){
                #Need larger space for ideogram:
                op$ideo.frac <- op$ideo.frac*2
            }
        }
    }else{
        op$ideo.frac <- 0
    } 
    
    #Check that we have a title for each plot:
    if(type=="genome"){
        op$main <- op$main[1]
    }
    if(type=="bychrom" && length(op$main)<length(chrom)){
        op$main <- rep(op$main[1],length(chrom))
    }

    return(op)
}

#Function that finds p-arm and chromosome stopping positions based on cytoband information

##Input:
### cyto.data: dataframe with cytoband information
### unit: the unit used to represent positions in data to be plotted (bp,kbp,mbp)

##Output:
### pstop: a vector giving the stopping positions for each p-arm (adjusted to match unit)
### chromstop: a vector giving the stopping position for each chromosome (adjusted to match unit)

##Required by:
### getArms
### addChromlines
### getGlobPos
### getGlobal.xlim

getArmandChromStop <- function(cyto.data, unit){
  
  #Sort cyto.data by chromosome number; let be represented by X=23 and Y=24:
  chrom <- cyto.data[,1]
  use.chrom <- gsub("chr","",chrom)  #Remove 'chr' from chrom-strings
  use.chrom[use.chrom=="X"] <- "23"	#Replace X by 23
  use.chrom[use.chrom=="Y"] <- "24"	#Replace Y by 24
  num.chrom <- as.numeric(use.chrom)	#Convert to numeric chromosomes
  
  #Order such that chromosomes are in increasing order from 1:24:
  ord.chrom <- order(num.chrom)
  cyto.data <- cyto.data[ord.chrom,,drop=FALSE] 	
  
  #Get chromosome stopping positions:
  chrom <- cyto.data[,1]
  chrom.stop <- which(chrom[1:length(chrom)-1]!=chrom[2:length(chrom)])
  chrom.stop <- c(chrom.stop,length(chrom))  #include last chromstop as well
  
  #Get p-arm stopping positions:
  arm.char <- substring(cyto.data[,4],1,1)   #Retrive first character in name which identifies p and q arms
  arm.stop <- which(arm.char[1:length(arm.char)-1]!=arm.char[2:length(arm.char)])
  p.stop <- arm.stop[-which(arm.stop%in%chrom.stop)]  #Remove qstops
  
  pos.chromstop <- cyto.data[chrom.stop,3]  #Local stopping position for each chromosome
  pos.pstop <- cyto.data[p.stop,3]		#Local stopping position for each p-arm
  
  #Factor used to convert positions into desired unit
  f <- switch(unit,
              bp = 1,
              kbp = 10^(-3),
              mbp = 10^(-6))
  
  return(list(pstop=pos.pstop*f,chromstop=pos.chromstop*f))
  
}

#Function to get a scaling factor such that positions may be converted from unit2 to unit1:

##Input:
### unit1: unit that we want to convert to
### unit2: unit that should be converted

##Output:
### factor: a scaling factor used in conversion of positions

##Required by:
### addChromlines
### chromMax
### getx
### plotIdeogram
### getGlobal.xlim

convert.unit <- function(unit1,unit2){
    
    factor <- NA
    #Allowed units:
    units <- c("bp","kbp","mbp")
    
    if(identical(unit1,unit2)){
        factor <- 1
    }else{
        if(identical(unit1,units[3])){
            if(identical(unit2,units[2])){
                factor <- 1/(10^3)
            }else{if(identical(unit2,units[1])){
                factor <- 1/(10^6)
            }}
        }else{
            if(identical(unit1,units[2])){
                if(identical(unit2,units[3])){
                    factor <- 10^3
                }else{if(identical(unit2,units[1])){
                    factor <- 1/(10^3)
                }}
            }else{
                if(identical(unit1,units[1])){
                    if(identical(unit2,units[3])){
                        factor <- 10^6
                    }else{if(identical(unit2,units[2])){
                        factor <- 10^3
                    }}
                }
            }
        }
    }
    
    if(is.na(factor)){
        if(all(units!=unit1)){
            stop("plot.unit must be one of 'kbp' and 'mbp'",call.=FALSE)
        }
        
    }
    
    return(factor)
}

#function that find global limits on xaxis (used in genome plots)

##Input: 
### op: list with plot parameters
### pos.unit: the unit used for positions in data
### chrom: a vector of unique chromosome numbers found in data

##Output:
### op: a list with updated plot parameters (xlim)

##Required by:
### plotFreq (genomeFreq)

##Requires:
### getArmandChromStop
### convert.unit

getGlobal.xlim <- function(op,pos.unit,chrom){  
    #Set xlim using chromosome information in cytoband; must transform to global information
    chromstop <- getArmandChromStop(op$assembly,pos.unit)$chromstop
    glob.chromstop <- cumsum(chromstop)   #Global stopping position for each chromosome
    scale.fac <- convert.unit(unit1=op$plot.unit,unit2=pos.unit)    #Scaling factor according to plot.unit
    #Not use chromosome X and Y if not in data
    if(!any(chrom==24)){
        glob.chromstop <- glob.chromstop[-24]
    }
    if(!any(chrom==23)){
        glob.chromstop <- glob.chromstop[-23]
    }
    xlim <- c(0,glob.chromstop[length(glob.chromstop)])*scale.fac
    
    return(xlim) 
}  

#Function to covert local posistions to global positions based on a defined set of local stop positions chromosomes given in cyto.data

##Input:
### chromosomes: vector with chromosome numbers
### position: vector with positions
### pos.unit: positon's unit
### cyto.data: data frame with cytoband information

##Output:
### glob.pos: vector with global positions

##Required by:
### adjustSeg
### getx 
### plotCircle

##Requires:
### getArmandChromStop

getGlobPos <- function(chromosomes, position, pos.unit, cyto.data){
    
    #Get local stopping posistions for each p-arm and each chromosome from cytoband data 
    l <- getArmandChromStop(cyto.data=cyto.data,unit=pos.unit)
    chromstop <- l$chromstop 
    
    #Need to make sure that none of the input positions are larger than chromosome stop postions:
    u.chrom <- unique(chromosomes)
    new.pos <- position
    for(j in 1:length(u.chrom)){
        probe.c <- chromosomes==u.chrom[j]
        #Check for positions that are larger than chromosome max position
        out <- which(position[probe.c] > chromstop[u.chrom[j]])
        #Replace these positions by max chrom position:
        new.pos[probe.c][out] <- chromstop[u.chrom[j]] 
    }
    glob.chromstop <- cumsum(chromstop)   #Global stopping position for each chromosome
    
    glob.pos <- new.pos + c(0,glob.chromstop[-length(glob.chromstop)])[chromosomes]  #Calculate global positions by adding global chromstop (for chrom > 1, for chrom=1 positions remain unchanged
    return(glob.pos)
    
}

# Function to return the values to be plotted along xaxis depending on the choice of xaxis (index or pos), and the type of plot
# If type=genome and xaxis=pos, global positions are calculated.
# Positions are scaled according to plotunit

### input:
### xaxis: index or pos
### type: plot type; genome, chromosome, sample or aspcf
### chromosomes: vector of same length as pos giving the corresponding chromosome numbers
### pos: vector giving probe positions 
### unit: unit used to represent positons   (bp,kbp, or mbp)
### op: a list containing other plot parameters

###Output:
### x: a vector containg the numbers to be plotted along xaxis of plot


##Required by:
### adjustPos

##Requires:
### getGlobPos
### convert.unit

getx <- function(xaxis,type,chromosomes,pos,unit,op){
    if(xaxis=="pos"){
        
        x <- pos
        if(type=="genome"){
            #Convert to global position:
            global.pos <- getGlobPos(chromosomes,pos,pos.unit=unit,cyto.data=op$assembly)   
            x <- global.pos
        }
        
        #Want to scale the x-axis to fit the desired unit given in plot.unit (default is mega base pairs)
        scale.fac <- convert.unit(unit1=op$plot.unit,unit2=unit)
        x <- x*scale.fac
        
    }else{
        #xaxis=="index"
        x <- 1:length(pos) 
    }#endif
    
    return(x)
    
}

#Function that returns the index where each chromosome starts (and the last chromosome ends)

##Input:
### v: a vector of chromosome numbers

## Output:
### cp: indeces for start of each chromosome and end of last chromosome

##Required by:
### addChromlines

separateChrom <- function(v){
    d <- diff(v)   #get difference between value (i+1) and value i in vector v
    cp <- which(d!=0)+1  #get changepoints
    
    #Add start of vector and (stop+1) of the whole vector
    cp <- c(1,cp,(length(v)+1))
    
    return(cp)
}

#Function that converts character arms to numeric

##Input:
### chrom : vector with chromosome numbers corresponding to each character arm
### char.arms : vector containing charcter arms; dentoed p or q

## Output:
### arms : vector with numeric arms calculated as chrom*2 - 1 or chrom*2

##Required by:
### adjustPos

numericArms <- function(chrom,char.arms){
    p.arm <- which(char.arms=="p")
    q.arm <- which(char.arms=="q")
    arms <- rep(NA,length(char.arms))
    arms[p.arm] <- chrom[p.arm]*2-1
    arms[q.arm] <- chrom[q.arm]*2
    
    return(arms)
}

#Function that checks if chrom is numeric, converts x/X and y/Y to 23 and 24 if not:

##Input:
### chrom: vector with chromosomes; numeric or character

## Output:
### chrom: numeric vector with chromosome numbers

##Required by:
### getArms

numericChrom <- function(chrom){ 
    if(!is.numeric(chrom)){
        if(is.factor(chrom)){
            #If chrom is factor; need to convert to character first
            chrom <- as.character(chrom)
        }
        #Replace X by 23:
        chrx <- c(which(chrom=="x"),which(chrom=="X"))
        chrom[chrx] <- 23
        #Replace Y by 24
        chry <- c(which(chrom=="y"),which(chrom=="Y"))
        chrom[chry] <- 24
        
        chrom <- as.numeric(chrom)
    }
    return(chrom)
}

#Get arm numbers from chromosomes, position and cytoband data:

##Input:
### chrom: vector of chromosome numbers
### pos: vector of positions
### pos.unit: unit for positions
### cyto.data: object specifying which genome assembly version should be used for cytoband data

##Output:
### arms: a character vector with arms; dentoed p and q 

##Required by:
### adjustPos
### multipcf
### pcf
### winsorize
### aspcf 


##Requires:
### getArmandChromStop
### numericChrom

getArms <- function(chrom, pos, pos.unit="bp", cyto.data){
    
    #Make sure chromosomes are numeric:
    chrom <- numericChrom(chrom)
    
    nProbe <- length(chrom)
    chrom.list <- unique(chrom)
    nChrom <- length(chrom.list)
    
    #Get local stopping posistions for each p-arm and each chromosome from cytoband data 
    l <- getArmandChromStop(cyto.data=cyto.data,unit=pos.unit)
    pStop <- l$pstop
    chromStop <- l$chromstop 
    
    #Intitialize
    arms <- rep(NA,nProbe)
    
    for(i in 1:nChrom){
        #Find corresponding arm numbers:
        c <- chrom.list[i]
        ind.c <- which(chrom==c)
        
        arms[ind.c] <- "q"
        p.arm <- ind.c[pos[ind.c]<=pStop[c]]
        arms[p.arm] <- "p"   #p-arm
        
    }
    
    return(arms)
}

#Function that scales positions according to plotunit, converts to global positions if type=genome, and finds start and stop positions (left and right) for recatangles to be plotted. Also makes sure that frequencies are shown as continuous if this is desired

##Input:
### position: the genomic postions to be plotted 
### chromosomes: the chromosomes corresponding to the positions
### pos.unit: the unit used for positions (bp,kbp,mbp)
### type: plot type (genome or chromosome)
###op: a list of other set plot parameters 

##Output:
### xleft: the left/start position of the plot rectangle
### xright: the right/stop position of the plot rectangle

##Required by:
### plotFreq (genomeFreq and chromosomeFreq)

##Requires:
### getx
### getArms
### numericArms
adjustPos <- function(position,chromosomes,pos.unit,type,op){
    
    if(type=="chromosome"){
        #Only need to scale positions first
        pos <- getx(xaxis="pos",type=type,chromosomes=NULL,pos=position,unit=pos.unit,op=op)
    }else if(type=="genome"){
        #Need to scale and convert to global pos:
        pos <- getx(xaxis="pos",type=type,chromosomes=chromosomes,pos=position,unit=pos.unit,op=op)
    }
    
    nPos <- length(position)
    #Define left-pos and right-pos for freqency-rectangle to be plotted:
    xleft <- pos
    xright <- pos
    
    #Should frequencies be plotted continously across probes?:
    if(op$continuous){
        #The rectangles should start and end halfway between two probes, except across different arms/chromosomes, fixing this below
        half <- (pos[2:nPos]-pos[1:(nPos-1)])/2
        xleft[2:nPos] <- xleft[2:nPos] - half
        xright[1:(nPos-1)] <- xright[1:(nPos-1)] + half   
    }else{
        #Let rectangle be one probe wide:
        xleft[2:nPos] <- xleft[2:nPos] - 0.5
        xright[1:(nPos-1)] <- xright[1:(nPos-1)] + 0.5
    }
    
    
    if(type!="genome"){
        #First find locations for change in arm number:
        char.arms <- getArms(chromosomes,position,pos.unit,op$assembly)
        arms <- numericArms(chromosomes,char.arms)
        #Where do arms start:
        n.arm <- length(unique(arms))
        
        if(n.arm>1){
            #Locate where arm starts (if more than one arm):
            sep.arm <- separateChrom(arms)
            sep.arm <- sep.arm[-c(1,length(sep.arm))]
            #Keep positions at arm-change 
            xleft[sep.arm] <- pos[sep.arm]
            xright[sep.arm-1] <- pos[sep.arm-1]
        }
    }else{
        #when type=genome: only want to prevent continous over chromosomes:
        n.chrom <- length(unique(chromosomes))
        if(n.chrom > 1){
            sep.chrom <- separateChrom(chromosomes)
            sep.chrom <- sep.chrom[-c(1,length(sep.chrom))]
            #keep positions at chrom change
            xleft[sep.chrom] <- pos[sep.chrom]
            xright[sep.chrom-1] <- pos[sep.chrom-1]   
        }
    }
    
    return(list(xleft=xleft,xright=xright))
    
}

#Function to find frame dimensions when a window is to be sectioned into nrow rows and ncol columns:

##Input:
## nrow: number of rows in plot
## ncol: number of columns in plot

##Output:
# a list giving the left, right, bottom and top dimensions for each of the nrow*ncol frames to be made in plot

##Required by:
### plotFreq (genomeFreq and chromosomeFreq)

framedim <- function(nrow,ncol){
    cl <- 0:(ncol-1)
    cr <- 1:ncol
    left <- rep(1/ncol,ncol)*cl
    right <- rep(1/ncol,ncol)*cr
    
    rt <- nrow:1
    rb <- (nrow-1):0
    top <- rep(1/nrow,nrow)*rt
    bot <- rep(1/nrow,nrow)*rb
    
    return(list(left=left,right=right,bot=bot,top=top))
}


#Function to set default ylim and at.y for plotFreq

##Input: 
### freq.del: vector with deletion frequencies
### freq.amp: vector with amplification frequencies
### op: list with plot parameters

##Output:
### op: list wiht updated plot parameters

##Required by: 
### plotFreq (genomeFreq and chromosomeFreq)

updateFreqParameters <- function(freq.del,freq.amp,op){
    #Y-limits; symmetric:
    max.freq <- max(c(freq.del,freq.amp))
    
    #Define tickmarks on y-axis
    if(max.freq>30){
        at.y <- seq(0,100,by=25)
    }else if(max.freq>10){
        at.y <- seq(0,100,by=10)
    }else{
        at.y <- seq(0,100,by=5)
    }
    if(is.null(op$at.y)){
        op$at.y <- c(at.y)
    }
    #Make sure ylim includes the first tickmark above max.freq
    q <- min(op$at.y[op$at.y>=max.freq])
    ylim <- c(-q,q)
    
    if(is.null(op$ylim)){
        op$ylim <- ylim
    }
    
    return(op)
}

# Function that separates chromosomes in genome-plots by color

## Required by:
## plotFreq (genomeFreq)

## Requires:
## getArmandChromStop
## convert.unit

chromPattern <- function(pos.unit,op) {
    #Use cytoband data information to get stopping points of chromosomes:
    chromstop <- getArmandChromStop(op$assembly,pos.unit)$chromstop
    scale.fac <- convert.unit(unit1=op$plot.unit,unit2=pos.unit)    #Scaling factor according to plot.unit
    chrom.mark <- c(1,cumsum(chromstop))*scale.fac
    #Drop chrom24 if no observations for this chrom in data/segments:
    #if(!any(segments[,2]==24)){
    # chrom.mark <- chrom.mark[-length(chrom.mark)]
    #} 
    #Let background be black to avoid white parts in arms without probes:
    #rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
    for (i in 1:(length(chrom.mark)-1)) {
        if(i%%2==0){
            rect(chrom.mark[i], par("usr")[3], chrom.mark[i+1], par("usr")[4], col = "grey95")#, border=NA)
        } 
    }
}

##Function that adds percentagelines, yaxis, xaxis and labels to frequency plots 
#Input:
### op: a list with plot parameters
### type: plot type (genome or bychrom)

##Required by:
### plotFreq (genomeFreq and chromosomeFreq)

##Requires:
### get.xticks
addToFreqPlot <- function(op,type){
    
    #Add y-lines if wanted
    if(is.logical(op$percentLines)){
        if(!op$percentLines){
            op$percentLines<- NULL
        }else{
            op$percentLines <- op$at.y
        }
    }
    if(!is.null(op$percentLines)){
        abline(h=c(-op$percentLines,op$percentLines),lty=3,col="grey82")
    }
    
    #Add yaxis and lab:
    #Make sure tickmarks are at percentLines
    if(!is.null(op$percentLines)){
        op$at.y <- op$percentLines
    }
    op$at.y <- c(-op$at.y,op$at.y)
    axis(side=2,cex.axis=op$cex.axis,at=op$at.y,mgp=op$mgp.y,las=op$las,tcl=-0.2,labels=abs(op$at.y))
    title(ylab=op$ylab,cex.lab=op$cex.lab,line=op$mgp.y[1])
    
    #Add xaxis:
    if(op$plot.ideo || type=="genome"){
        axis(side=1,labels=FALSE,tcl=0)
    }else{
        if(is.null(op$at.x)){
            op$at.x <- get.xticks(0,op$xlim[2],unit=op$plot.unit,ideal.n=6)
        }
        axis(side=1,tcl=-0.2,at=op$at.x,cex.axis=op$cex.axis,mgp=op$mgp)
        title(xlab=op$xlab,cex.lab=op$cex.lab,line=op$mgp[1])
    }
    
    
}

#Function used to separate chromosomes by stapled lines in genome plot
##Input:
### chromosomes: a vector of chromosome numbers
### xaxis: what should be plotted along xaxis; pos or index
### unit: the unit used for positions, bp, kbp or mbp
### ind: only used when called by plotSegments; vector with start position and last stop position for segments
### cex: size used to plot chromosome numbers
### op: other plot parameters

##Output:
### adds chromosome lines to existing genome plot

##Required by:
### plotFreq

##Requires:
### convert.unit
### getArmandChromStop
### separateChrom

addChromlines <- function(chromosomes,xaxis,unit,ind=NULL,cex,op){
    if(xaxis=="pos"){
        #Use cytoband data information to get stopping points of chromosomes:
        chromstop <- getArmandChromStop(op$assembly,unit)$chromstop
        scale.fac <- convert.unit(unit1=op$plot.unit,unit2=unit)    #Scaling factor according to plot.unit
        chrom.mark <- c(1,cumsum(chromstop))*scale.fac
        #Drop chrom24 if no observations for this chrom in data/segments:
        if(any(chromosomes==24)){
            chromosomes <- 1:24
        }else{
            chromosomes <- 1:23
            chrom.mark <- chrom.mark[-length(chrom.mark)]
        }
    }else{
        chrom.mark <- separateChrom(chromosomes)
        if(!is.null(ind)){
            #Segments are used to decide where to separate chromosomes; get index start where chromosome number changes and last index
            chrom.mark <- ind[chrom.mark]
        }
        chrom.mark <- chrom.mark - 0.5
    }
    #Separate chromosomes by vertical lines in existing plot:
    nChrom <- length(unique(chromosomes))
    arg <- list(chrom.lwd=1, chrom.lty=2, chrom.col="darkgrey",chrom.side=3, chrom.cex=cex,chrom.line=c(0,0.3))
    
    if(!is.null(op)){
        arg <- modifyList(arg,op)
    }
    abline(v=chrom.mark[2:(length(chrom.mark)-1)],col=arg$chrom.col,lwd=arg$chrom.lwd,lty=arg$chrom.lty)
    
    at <- (chrom.mark[1:nChrom]-1)+(chrom.mark[2:(nChrom+1)]-chrom.mark[1:nChrom])/2
    chrom.names <- unique(chromosomes)
    chrom.names <- gsub("23","X", chrom.names)
    chrom.names <- gsub("24","Y", chrom.names)
    #Plot half at bottom, half at top:
    bot <- seq(1,length(chrom.mark),2)
    top <- seq(2,length(chrom.mark),2)
    mtext(chrom.names[bot],side=1,line=arg$chrom.line[1],at=at[bot],cex=arg$chrom.cex)
    mtext(chrom.names[top],side=3,line=arg$chrom.line[2],at=at[top],cex=arg$chrom.cex)
    
}

# Function that separates chromosome arms in genome-plots by dashed lines

## Required by:
## plotFreq (genomeFreq)

## Requires:
## getArmandChromStop
## convert.unit

addArmlines <- function(chromosomes,xaxis,unit,ind=NULL,cex,op){
    if(xaxis=="pos"){
        #Use cytoband data information to get stopping points of chromosome arms:
        marks <- getArmandChromStop(op$assembly,unit)
        armstop <- c(marks$pstop[1],cumsum(marks$chromstop)[1:length(marks$chromstop)-1]+marks$pstop[2:length(marks$pstop)])
        scale.fac <- convert.unit(unit1=op$plot.unit,unit2=unit)    #Scaling factor according to plot.unit
        arm.mark <- armstop*scale.fac
        #Separate arms by vertical lines in existing plot:
        arg <- list(chrom.lwd=1, chrom.lty=2, chrom.col="darkgrey",chrom.side=3, chrom.cex=cex,chrom.line=c(0,0.3))
        
        if(!is.null(op)){
            arg <- modifyList(arg,op)
        }  
        abline(v=arm.mark[1:(length(arm.mark)-1)],col=arg$chrom.col,lwd=arg$chrom.lwd,lty=2)
    }
}


#Function that plots ideogram for given chromosome

#Input:
### chrom: number between 1-24 indicating which chromosome's ideogram is to be plotted
### cyto.text: should cytoband-names be plotted along ideogram?
### cex: the size used for cyto.text
### cyto.data: data frame with cytoband-information
### cyto.unit: the unit used to represent positons in cyto.data
### unit: the unit used for positions in the plot


##Required by:
### plotFreq (chromosomeFreq)

##Requires:
### convert.unit

plotIdeogram <- function(chrom,cyto.text=FALSE,cex=0.6,cyto.data,cyto.unit="bp",unit){
    
    if(chrom==23){
        chrom.cytoband <- cyto.data[cyto.data[,1]=="chrX",]
    }else{
        if(chrom==24){
            chrom.cytoband <- cyto.data[cyto.data[,1]=="chrY",]
        }else{
            chrom.cytoband <- cyto.data[cyto.data[,1]==paste("chr",chrom,sep=""),]
        }
    }
    
    cyto.start <- chrom.cytoband[,2]
    cyto.end <- chrom.cytoband[,3]
    scale <- convert.unit(unit1=unit,unit2=cyto.unit)
    
    xleft <- cyto.start*scale
    xright <- cyto.end*scale
    n <- length(xleft)
    chrom.length <- xright[n]-xleft[1]
    
    stain <- chrom.cytoband[,5]
    sep.stain <- c("gpos","gneg","acen","gvar","stalk")
    
    g <- sapply(sep.stain,grep,x=stain,fixed=TRUE)
    
    centromere <- g$acen
    stalk <- g$stalk
    col <- rep("",n)
    col[stain=="gneg"] <- "white"
    col[stain=="gpos100"] <- "black"
    col[stain=="gpos75"] <- "gray25"
    col[stain=="gpos50"] <- "gray50"
    col[stain=="gpos25"] <- "gray75"
    col[stain=="stalk"] <- "gray90"
    col[stain=="gvar"] <- "grey"
    col[stain=="acen"] <- "yellow"
    density <- rep(NA,n)
    angle <- rep(45,n)
    density[stain=="gvar"] <- 15
    
    
    ylow <- 0
    yhigh <- 1
    
    
    plot(x=c(0,max(xright)),y=c(ylow,yhigh),type="n",axes=FALSE,xlab="",ylab="",xlim=c(0,max(xright)),ylim=c(0,1),xaxs="i")
    
    #Rectangles:
    skip.rect <- c(1,centromere,n,stalk)
    rect(xleft[-skip.rect],rep(ylow,n-length(skip.rect)),xright[-skip.rect],rep(yhigh,n-length(skip.rect)),
         col=col[-skip.rect],border="black",density=density[-skip.rect],angle=angle[-skip.rect])
    
    #Round edges at ideogram start, stop and at centromere:
    draw.roundEdge(start=xleft[1],stop=xright[1],y0=ylow,y1=yhigh,col=col[1],bow="left",density=density[1],angle=angle[1],chrom.length=chrom.length)
    draw.roundEdge(start=xleft[centromere[1]],stop=xright[centromere[1]],y0=ylow,y1=yhigh,col=col[centromere[1]],bow="right",density=density[centromere[1]],
                   angle=angle[centromere[1]],lwd=1,chrom.length=chrom.length)
    draw.roundEdge(start=xleft[centromere[2]],stop=xright[centromere[2]],y0=ylow,y1=yhigh,col=col[centromere[2]],bow="left",density=density[centromere[2]],
                   angle=angle[centromere[2]],lwd=1,chrom.length=chrom.length)
    draw.roundEdge(start=xleft[n],stop=xright[n],y0=ylow,y1=yhigh,col=col[n],bow="right",density=density[n],angle=angle[n],chrom.length=chrom.length)
    
    #Draw stalk-segment:
    if(length(stalk)>0){
        for(i in 1:length(stalk)){
            drawStalk(xleft[stalk[i]],xright[stalk[i]],ylow,yhigh,col=col[stalk[i]])
        }
    }
    if(cyto.text){
        mtext(text=paste(chrom.cytoband[,4],"-",sep=" "),side=1,at=(xleft + (xright-xleft)/2),cex=cex,las=2,adj=1,xpd=NA)#,line=-1)#,outer=TRUE)
    }
    
}

# Function for plotIdeogram
draw.roundEdge <- function(start,stop,y0,y1,col,bow,density=NA,angle=45,lwd=1,chrom.length){
    #Y points in round edge:
    f <- rep(0,0)
    f[1] <- 0.001
    i=1
    half <- y0+(y1-y0)/2
    while(f[i]<half){
        f[i+1] <- f[i]*1.3
        i <- i+1
    }
    f <- f[-length(f)]
    
    Y <- c(y1,y1,y1-f,half,y0+rev(f),y0,y0)
    
    #X points in roundedge
    cyto.length <- stop-start
    
    share <- cyto.length/chrom.length
    if(share>0.2){
        #to create bow in end of chromosome 24
        share <- 0.2
    }
    
    if(bow=="left"){
        
        round.start <- start + cyto.length*(1-share)^20
        
        x <- seq(round.start,start,length.out=(length(f)+2))
        revx <- rev(x[-length(x)])
        x <- c(x,revx)
        X <- c(stop,x,stop)
    }else{
        if(bow=="right"){
            round.start <- stop - cyto.length*(1-share)^20
            x <- seq(round.start,stop,length.out=(length(f)+2))
            revx <- rev(x[-length(x)])
            x <- c(x,revx)
            X <- c(start,x,start)
            
        }
    }
    
    polygon(x=X,y=Y,col=col,border="black",density=density,angle=angle,lwd=lwd)
    
}

# Function for plotIdeogram
drawStalk <- function(start,stop,y0,y1,col){
    size <- stop-start
    x1 <- c(start,start+size/3,stop-size/3,stop)
    x2 <- rev(x1)
    x <- c(x1,x2)
    y_1 <- c(y0,y0+0.25,y0+0.25,y0)
    y_2 <- c(y1,y1-0.25,y1-0.25,y1)
    y <- c(y_1,y_2)
    polygon(x=x,y=y,col=col)
    
}

#Get the maximum position on a given chromosome

#Input:
### chrom: number between 1-24 indicating which chromosome's ideogram is to be plotted
### cyto.data: data frame with cytoband-information
### unit: the unit used for positions in the plot
### cyto.unit: the unit used to represent positons in cyto.data

# Output:
### max.pos: a scalar giving the maximum position on this chromosome (given the cytoband information)

##Required by:
### plotFreq (chromosomeFreq)

## Requires:
### convert.unit
chromMax <- function(chrom,cyto.data,pos.unit,cyto.unit="bp"){
    
    #Get scaling factor for positions
    s <- convert.unit(unit1=pos.unit,unit2=cyto.unit)
    
    #Get the rows in cytoband data that correspond to this chromosome
    if(chrom==23){
        txt <- "chrX"
    }else if(chrom==24){
        txt <- "chrY"
    }else{
        txt <- paste("chr",chrom,sep="") 
    }
    rows <- which(cyto.data[,1]==txt)
    
    ##Get max position for this chromosome and scale according to pos.unit
    max.pos <- max(cyto.data[rows,3])
    max.pos <- max.pos*s
    
    return(max.pos)
    
}
