generate_genomic_intervals <- function(genome='hg38', bin.size= 1000000, soft.expansion = 100000){
  cb <- data.table::fread(paste0('http://hgdownload.cse.ucsc.edu/goldenpath/',genome,'/database/cytoBand.txt.gz'), header = F, data.table = FALSE)
  chroms <- c(paste0("chr", 1:22), "chrX", "chrY")
  cb <- cb[cb$V1 %in% chroms, ]
  rownames(cb) <- seq(1:nrow(cb))
  arm.length <- numeric(48)
  
  for (i in seq(24)){
    temp <- cb[cb$V1 == chroms[i],]
    j <- 2*(i-1)
    arm.length[j+1] <- max(temp[temp$V4 %in% c('p11.1','p11','p11.11'),'V3'])
    arm.length[j+2] <- max(temp$V3)
  }
  
  bin.size = bin.size
  soft.expansion = soft.expansion
  count = 1
  result <- list()
  for ( i in seq(24)){
    p_max = arm.length[2*i-1]
    q_max = arm.length[2*i]
    arm = "p"
    start = 0
    p_first = p_max
    chr_result = list()
    chr_count = 1
    while (p_first >= bin.size + soft.expansion){
      p_first = p_first-bin.size
    }
    end = start + p_first
    while(start < q_max){
      interval.len = bin.size
      if (end > q_max){
        end = q_max
      }else if (q_max < end + soft.expansion){
        end = q_max
        interval.len = interval.len + soft.expansion
      }
      
      if (end >= p_max){
        arm = "q"
      }
      
      size = end-start
      chr_result[[chr_count]] <- data.frame(index= count, chromosome = gsub("chr","",chroms[i]), start=start, end=end, arm=arm, size= size )
      start = end
      end = end+interval.len
      count = count+ 1
      chr_count = chr_count+1
    }
    result[[chroms[i]]] <- Reduce(rbind, chr_result)
  }
  result <- Reduce(rbind,result)
  return(result)
}

hg19_bins <- generate_genomic_intervals(genome = 'hg19')
hg38_bins <- generate_genomic_intervals(genome = 'hg38')
