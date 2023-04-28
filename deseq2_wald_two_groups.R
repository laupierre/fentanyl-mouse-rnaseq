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

write.xlsx (a, "./output/star_gene_raw_counts.xlsx", rowNames=F)
# excel file renamed as hippocampus_star_gene_raw_counts.xlsx


a <- a[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", a$gene_type, invert=TRUE), ]
colnames (a) <- gsub ("Aligned.out.bam", "", colnames (a))
colnames (a) <- gsub ("star.", "", colnames (a))


## DESeq2 

counts <- annot <- a

annot <- annot[ ,c("Geneid", "gene_name", "gene_type", "mgi_id", "external_gene_name", "description")]

row.names (counts) <- counts$Geneid
counts <- counts[ ,grep ("S", colnames (counts))]

samples <- data.frame (matrix (nrow=dim (counts)[2], ncol=3))
colnames (samples) <- c("sample", "condition", "sex")
samples$sample <- colnames (counts)
samples$sex <- as.factor (substr (gsub ("\\..*" , "", colnames (counts)), 1, 1))
samples$condition <- as.factor (substr (gsub ("\\..*" , "", colnames (counts)), 2, 2))

samples <- samples[!samples$sample %in% c("mF.01_S10", "fF.02_S11"), ]
samples

counts <- counts[ ,colnames (counts) %in% samples$sample]

idx <- match (samples$sample, colnames (counts))
samples <- samples[idx, ]
stopifnot (samples$sample == colnames (counts))


dds <- DESeqDataSetFromMatrix(countData = round (counts), colData = samples, design = ~ sex +condition)
                                 
# keep <- rowSums(counts(dds)) >= 10
keep <- rowSums(counts(dds) >= 10) >= dim (counts)[2]/2
dds <- dds[keep,]
dds

# first contrast of interest (F vs C)
ddsLRT <- DESeq(dds, test="LRT", full=~sex+condition, reduced=~sex)
resultsNames(ddsLRT)

res <- results(ddsLRT, contrast=list("condition_F_vs_C"), test="Wald")

res <- merge (data.frame (res), counts (dds), by="row.names")
#res <- merge (data.frame (res), round (counts (dds, normalized=TRUE)), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

write.xlsx (res, "hippocampus_deseq2_FentanylvsControl_differential_expression.xlsx", rowNames=F)


# second contrast of interest (W vs C)
ddsLRT <- DESeq(dds, test="LRT", full=~sex+condition, reduced=~sex)
resultsNames(ddsLRT)

res <- results(ddsLRT, contrast=list("condition_W_vs_C"), test="Wald")

res <- merge (data.frame (res), counts (dds), by="row.names")
#res <- merge (data.frame (res), round (counts (dds, normalized=TRUE)), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

write.xlsx (res, "hippocampus_deseq2_WithdrawalvsControl_differential_expression.xlsx", rowNames=F)


# third contrast of the two previous interaction terms, using the list() style of contrasts !!!
ddsLRT <- DESeq(dds, test="LRT", full=~sex+condition, reduced=~sex)
resultsNames(ddsLRT)

res <- results(ddsLRT, contrast=list("condition_W_vs_C", "condition_F_vs_C"), test="Wald")  ## This is equivalent to W/F

res <- merge (data.frame (res), counts (dds), by="row.names")
#res <- merge (data.frame (res), round (counts (dds, normalized=TRUE)), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

write.xlsx (res, "hippocampus_deseq2_WithdrawalvsFentanyl_differential_expression.xlsx", rowNames=F)


## PCA plot

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition, shape=sex)) +
  		geom_point(size=3) +
  		xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  		ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
		coord_fixed ()











