library(msigdbr)
library (openxlsx)
library (fgsea)

# The Gene Set Enrichment Analysis (GSEA) differs from Gene Ontology enrichment analysis in that it considers all genes  
# in contrast to taking only significantly differentially expressed genes. 



mygsea <- function (res, outfile) {

## BP from GO
#pathways <- msigdbr("mouse", category="C5", subcategory = "GO:BP")
#pathways <- split (as.character (pathways$ensembl_gene), pathways$gs_name)

#data.frame (msigdbr_collections())

# Pathways
pathways <- msigdbr("mouse", category="C2", subcategory = "CP:PID")
pathways <- split (as.character (pathways$ensembl_gene), pathways$gs_name)

# ranks are from lowest to highest
res <- res[order (res$log2FoldChange), ]
ranks <- res$log2FoldChange
names (ranks) <- gsub ("\\..*", "", res$Geneid)
  
fgseaRes <- fgsea(pathways = pathways, 
                  stats    = ranks,
                  minSize  = 10,
                  maxSize  = 500)
  
# see 20 pathways ordered by pval  
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam=0.5)


# Enrichment plot of a specific top pathway (up and down)
library (ggplot2)
plotEnrichment(pathways[[topPathwaysUp[[1]] ]], ranks) + ggtitle( topPathwaysUp[[1]] )

plotEnrichment(pathways[[topPathwaysDown[[1]] ]], ranks) + ggtitle( topPathwaysDown[[1]] )


# make leading edge more human-readable
library (org.Mm.eg.db)

fgseaRes[ , leadingEdge := mapIdsList(x=org.Mm.eg.db, 
                                      keys=leadingEdge,
                                      keytype="ENSEMBL", 
                                      column="SYMBOL")]


# save fgsea results in a file
fgseaRes <- fgseaRes[ ,c(1:6,8)]
fgseaRes <- data.frame (fgseaRes[order (fgseaRes$padj), ])
fgseaRes <- fgseaRes[order (fgseaRes$padj), ]
fgseaRes <- fgseaRes[ ,-4]
fgseaRes <- fgseaRes[ ,-4]
head (fgseaRes)

# In 8, we define the leading-edge subset to be those genes in the gene set S  
# that appear in the ranked list L at, or before, the point where the running sum reaches its maximum deviation from zero.

# Positive NES correspond to up-regulated gene sets whereas negative NES to down-regulated ones

write.xlsx (fgseaRes, outfile, rowNames=F)
}


# Fentanyl vs Control
res <- read.xlsx ("hippocampus_deseq2_FentanylvsControl_differential_expression.xlsx")
mygsea (res, "hippocampus_deseq2_FentanylvsControl_gsea_pathways.xlsx")

# Withdrawal vs Control
res <- read.xlsx ("hippocampus_deseq2_WithdrawalvsControl_differential_expression.xlsx")
mygsea (res, "hippocampus_deseq2_WithdrawalvsControl_gsea_pathways.xlsx")

# Withdrawal vs Fentanyl
res <- read.xlsx ("hippocampus_deseq2_WithdrawalvsFentanyl_differential_expression.xlsx")
mygsea (res, "hippocampus_deseq2_WithdrawalvsFentanyl_gsea_pathways.xlsx")








