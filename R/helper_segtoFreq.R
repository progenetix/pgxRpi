generate_genomic_intervals <- function(genome='hg38', bin.size= 1000000, soft.expansion = 0.1){
  data(list=c(paste0(genome,'_cytoband')))
  cb <- get(paste0(genome,'_cytoband'))
  arm.length <- numeric(48)
  chroms <- c(paste0("chr", 1:22), "chrX", "chrY")
  for (i in seq(24)){
    temp <- cb[cb$V1 == chroms[i],]
    j <- 2*(i-1)
    arm.length[j+1] <- max(temp[temp$V4 %in% c('p11.1','p11','p11.11'),'V3'])
    arm.length[j+2] <- max(temp$V3)
  }
  
  bin.size <- bin.size
  soft.expansion <- soft.expansion * bin.size
  count <- 1
  result <- list()
  for ( i in seq(24)){
    p_max <- arm.length[2*i-1]
    q_max <- arm.length[2*i]
    arm <- "p"
    start <- 0
    p_first <- p_max
    chr_result <- list()
    chr_count <- 1
    while (p_first >= bin.size + soft.expansion){
      p_first <- p_first-bin.size
    }
    end <- start + p_first
    while(start < q_max){
      interval.len <- bin.size
      if (end > q_max){
        end <- q_max
      }else if (q_max < end + soft.expansion){
        end <- q_max
        interval.len <- interval.len + soft.expansion
      }
      
      if (end >= p_max){
        arm <- "q"
      }
      
      size <- end-start
      chr_result[[chr_count]] <- data.frame(index= count, chromosome = gsub("chr","",chroms[i]), start=start, end=end, arm=arm, size= size )
      start <- end
      end <- end+interval.len
      count <- count+ 1
      chr_count <- chr_count+1
    }
    result[[chroms[i]]] <- Reduce(rbind, chr_result)
  }
  result <- Reduce(rbind,result)
  return(result)
}

extract.bin.feature <- function(data,genome,bin.size,overlap,soft.expansion){
    bins <- generate_genomic_intervals(genome = genome, bin.size = bin.size, soft.expansion = soft.expansion)
    exclude.sexchr <- !c('X','Y') %in% unique(data[,2])
    if (all(exclude.sexchr)){
      bins <- bins[!bins[,2] %in% c('X','Y'),]
    }
        
    total_dup <- list()
    total_del <- list()   
    
    for (sample_idx in seq_len(length(unique(data[,1])))){
        ind.data <- data[data[,1] %in% unique(data[,1])[sample_idx],]
        ind_dup <- rep(0,dim(bins)[1])
        ind_del<- rep(0,dim(bins)[1])
        for (j in seq_len(dim(ind.data)[1])){
            ind.seg.start <- ind.data[j,3]
            ind.seg.end <- ind.data[j,4]
            sel.bin <- which(bins$chromosome == ind.data[j,2] & bins[,4] > ind.seg.start & bins[,3] < ind.seg.end)
            overlap.dist <- sapply(sel.bin,function(x){min(bins[x,4],ind.seg.end)-max(bins[x,3],ind.seg.start)})
            sel.bin <- sel.bin[overlap.dist >= overlap]
            if (length(sel.bin) == 0){next}
            ind_dup[sel.bin] <- ind_dup[sel.bin] + as.numeric(ind.data[j,6] == 'DUP')
            ind_del[sel.bin] <- ind_del[sel.bin] + as.numeric(ind.data[j,6] == 'DEL')
        }
        ind_dup[ind_dup > 1] <- 1
        ind_del[ind_del > 1] <- 1
        total_dup[[sample_idx]] <- ind_dup
        total_del[[sample_idx]] <- ind_del
    }
    total_dup <- do.call(rbind,total_dup)
    total_del <- do.call(rbind,total_del)
    
    rownames(total_dup) <- unique(data[,1])
    rownames(total_del) <- unique(data[,1])
    
    feature.list <- list()
    feature.list[['dup']] <- total_dup
    feature.list[['del']] <- total_del
    feature.list[['bin']] <- bins
    return(feature.list)
}



