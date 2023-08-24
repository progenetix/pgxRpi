# download cytoband data of GRCh38 
cytoband_hg38 <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz"
cb <- data.table::fread(cytoband_hg38, header = FALSE, data.table = FALSE)
cb <- cb[cb$V1 %in% c(paste0("chr", 1:22), "chrX", "chrY"), ]
rownames(cb) <- seq(1:nrow(cb))
cb$V1 <- factor(cb$V1)
cb$V4 <- factor(cb$V4)
cb$V5 <- factor(cb$V5)
hg38 <- cb

# download cytoband data of GRCh37 
cytoband_hg19 <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz"
cb <- data.table::fread(cytoband_hg19, header = FALSE, data.table = FALSE)
cb <- cb[cb$V1 %in% c(paste0("chr", 1:22), "chrX", "chrY"), ]
rownames(cb) <- seq(1:nrow(cb))
cb$V1 <- factor(cb$V1)
cb$V4 <- factor(cb$V4)
cb$V5 <- factor(cb$V5)
hg19 <- cb