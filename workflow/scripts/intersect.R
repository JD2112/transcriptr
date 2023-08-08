#!/usr/bin/env Rscript
library(data.table)
library(ggvenn)

edgeR <- fread("edgeR_results/COND_vs_CTR_filtered.tsv")
DESeq2 <- fread("deseq2_results/COND_vs_CTR_filtered.tsv")
limma <- fread("limma_results/COND_vs_CTR_filtered.tsv")

l <- list(edgeR, DESeq2, limma)
dtm = Reduce(function(...) merge(..., all = TRUE), l)


i <- list(
  edgeR = edgeR$id, 
  DESeq2 = DESeq2$id,
  limma = limma$id
)


v <- ggvenn(
  i, 
  fill_color = c("#0073C2FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)


ids <- intersect(intersect(edgeR$id, DESeq2$id), limma$id)
t <- dtm[id %in% ids ,]

colnames(t) <- c('id', 'name', 'chr', 'start', 'end', 'strand', 'biotype',
	'edgeR_logFC', 'edgeR_logCPM', 'edgeR_PValue', 'edgeR_FDR', 'edgeR_baseMean',	
	'DESeq2_log2FoldChange', 'DESeq2_lfcSE', 'DESeq2_stat', 'DESeq2_pvalue', 'DESeq2_padj',	
	'limma_logFC', 'limma_AveExpr', 'limma_t', 'limma_P.Value', 'limma_adj.P.Val', 'limma_B')

t <- t[order(edgeR_FDR, decreasing=FALSE),]

ggsave(filename="./venn.pdf", plot=v, device="pdf", width=29.7, height=21, dpi=300, units = "cm")
write.table(t, file="./intersection.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

