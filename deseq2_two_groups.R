library (openxlsx)
library (DESeq2)
library (ggplot2)

system ("mkdir ./output")

## See Github RNA-Seq_mouse/gene_annotation.R
anno <- read.delim ("gencode.vM32.annotation.txt")

anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)

a <- read.delim ("./projects/star_results/subread.counts.txt", skip=1)
a <- a[ ,grep ("Gene|bam", colnames (a))]

a <- merge (a, anno, by.x="Geneid", by.y="gene_id", all.x=TRUE) 

a <- a[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", a$gene_type, invert=TRUE), ]
colnames (a) <- gsub ("Aligned.out.bam", "", colnames (a))
colnames (a) <- gsub ("star.", "", colnames (a))

  
write.xlsx (a, "./output/star_gene_raw_counts.xlsx", rowNames=F)

