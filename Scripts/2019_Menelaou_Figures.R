#!/usr/bin/Rscript
rm(list=ls())


#-------------------------------------------------------------------------------------------------#
# Author  : Dr Malwina Prater
# Email   :     mn367@cam.ac.uk
#
# Project : CTR_edw23_0002
# PI      : Erica Watson (edw23@cam.ac.uk)
#
# R-Script to perform differential transcript analysis, produce figures and enrichment tables 
#
#-------------------------------------------------------------------------------------------------#

library("DESeq2")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("ggrepel")
library("dplyr")
library("cowplot")
library("pheatmap")
library("DOSE")
library("clusterProfiler")
library("org.Mm.eg.db")
library("karyoploteR")
library("TxDb.Mmusculus.UCSC.mm10.ensGene")
library("biomaRt")
library("ComplexHeatmap")
library("circlize")


message("+-------------------------------------------------------------------------------")
message("+ Set up some constants e.g. base directories")
message("+-------------------------------------------------------------------------------")


Project  <- "CTR_edw23_0002_TRANSFER_vs_NON_TRANSFER"
Base.dir <- "/Users/malwina/Documents/CTR-Groups/Erica_Watson/CTR_edw23_0002"
setwd(Base.dir)
list.files(Base.dir)
HTSeq.dir <- paste(Base.dir,"/combined_HTSEQ-counts", sep="")
list.files(HTSeq.dir)
Karyo.dir <- "/Users/malwina/Documents/CTR-Groups/Erica_Watson/CTR_edw23_0002/KARYOPLOTS/Karyoplots/"

elementTextSize <- 10

significance <- 0.05
l2fc         <- 1 

#Res.dir <- "/Users/malwina/Documents/CTR-Manuscripts/2019_Menelaou/2019-Menelaou"
#setwd(Res.dir)

message("+-------------------------------------------------------------------------------")
message("+ Set up the sample table")
message("+-------------------------------------------------------------------------------")


sample_table <-  read.csv(file.path("/Users/malwina/Documents/CTR-Manuscripts/2019_Menelaou/2019-Menelaou/Tables/SampleTable.csv"),  header = TRUE) # sep = "," or "\t"


message("+-------------------------------------------------------------------------------")
message("+ Retrieve ensEMBL annotations")
message("+-------------------------------------------------------------------------------")

ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description'), mart = ensembl)          
head(ensEMBL2id)
attributes_list <- listAttributes(ensembl)

ensEMBL2id2 <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', 'entrezgene'), mart = ensembl) 
ensEMBL2id2_bed <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position'), mart = ensembl) 
ensEMBL2id2_granges <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene','chromosome_name', 'start_position', 'end_position'), mart = ensembl) 

#write.csv(ensEMBL2id2, "ensEMBL2id2.csv")
#ensEMBL2id <- read.csv("ensEMBL2id.csv")

message("+-------------------------------------------------------------------------------")
message("+ Create ddsHTSeq object")
message("+-------------------------------------------------------------------------------")

ddsHTSeq<- DESeqDataSetFromHTSeqCount(sampleTable=sample_table, directory=HTSeq.dir, design= ~ Condition)
dds <- DESeq(ddsHTSeq)
colData(ddsHTSeq)

# shrinkage included:
res <- lfcShrink(dds, coef="Condition_transferred_vs_non.transferred", type="normal")

res<-res[order(res$padj),]
head(assay(dds))
colSums(assay(dds))
resultsNames(dds) # "Condition_transferred_vs_non.transferred"

dds <- estimateSizeFactors(dds)


message("+-------------------------------------------------------------------------------")
message("+ Create results object")
message("+-------------------------------------------------------------------------------")

results.df <- data.frame(results(dds, contrast= c("Condition", "transferred", "non-transferred")))
results.df$gene <- rownames(results.df)

counts <- counts(dds)
head(counts)
normCounts <- counts(dds, normalized=TRUE)
head(normCounts)



message("+-------------------------------------------------------------------------------")
message("+ Run transformations")
message("+-------------------------------------------------------------------------------")

rld <- rlogTransformation(dds, blind=T)


message("+-------------------------------------------------------------------------------")
message("+ Create PCA Plots")
message("+-------------------------------------------------------------------------------")

elementTextSize <- 12
topNum = 5000

pca = prcomp(t(assay(rld)))
rv = rowVars(assay(rld))
select = order(rv, decreasing = TRUE)[seq_len(min(topNum, length(rv)))]
pca = prcomp(t(assay(rld)[select, ]))


pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")

scores <- data.frame(sample_table$fileName, pca$x, sample_table$Condition)

#pdf(paste(Project, "_DESeq2_PCA_top", topNum, "MV.pdf", sep=""),width=8,height=7)
#par(bg=NA)
#ggplot(scores, aes(x = PC1, y = PC2, col = (factor(sample_table$Condition))) ) +
#  geom_point(size = 5 ) + 
#  xlab(pc1lab) + ylab(pc2lab) + 
#  scale_colour_manual(name="Pedigree", values = c("seagreen2","royalblue1",  "blue4", "#5f022f", "violetred", "olivedrab", "cornflowerblue", "blue")) +
#  scale_shape_discrete(name="Phenotype") +
#  theme(text = element_text(size=elementTextSize)) 
#dev.off()

Fig_PCA <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(sample_table$Condition))) ) +
  geom_point(size = 5 ) + 
  xlab(pc1lab) + ylab(pc2lab) + 
  scale_colour_manual(name="Pedigree", values = c("seagreen2","royalblue1",  "blue4", "#5f022f", "violetred", "olivedrab", "cornflowerblue", "blue")) +
  scale_shape_discrete(name="Phenotype") +
  theme(text = element_text(size=elementTextSize)) 


message("+-------------------------------------------------------------------------------")
message("+                                PCA explained                                  ")
message("+-------------------------------------------------------------------------------")

ensanno <- ensEMBL2id[,c(1:2)]
ensanno <- ensanno[!duplicated(ensanno),]
head(ensanno)
#ensanno <- ensanno[-1,]

loadings                         <- as.data.frame(pca$rotation)
loadings$ensembl_gene_id         <- rownames(loadings)
loadings                         <- merge(loadings, ensanno, by="ensembl_gene_id")
head(loadings)

pca.1         <-  loadings[ order(loadings$PC1,decreasing=TRUE), ]
pca.1.25      <-  pca.1[c(1:25),]
pca.1.25.plot <- ggplot(data=pca.1.25, aes(x=factor(pca.1.25$external_gene_name,levels=unique(pca.1.25$external_gene_name)), y=PC1)) + geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

pca.2         <-  loadings[ order(loadings$PC2,decreasing=TRUE), ]
pca.2.25      <-  pca.2[c(1:25),]
pca.2.25.plot <- ggplot(data=pca.2.25, aes(x=factor(pca.2.25$external_gene_name,levels=unique(pca.2.25$external_gene_name)), y=PC2)) + geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

pca.3         <-  loadings[ order(loadings$PC3,decreasing=TRUE), ]
pca.3.25      <-  pca.3[c(1:25),]
pca.3.25.plot <- ggplot(data=pca.3.25, aes(x=factor(pca.3.25$external_gene_name,levels=unique(pca.3.25$external_gene_name)), y=PC3)) + geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

pca.4         <-  loadings[ order(loadings$PC4,decreasing=TRUE), ]
pca.4.25      <-  pca.4[c(1:25),]
pca.4.25.plot <- ggplot(data=pca.4.25, aes(x=factor(pca.4.25$external_gene_name,levels=unique(pca.4.25$external_gene_name)), y=PC4)) + geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


Fig_PCA_loadings <- plot_grid(pca.1.25.plot, pca.2.25.plot, labels=c(" ", " "), ncol = 1, nrow = 2)


Fig_PCA_with_loadings <- plot_grid(Fig_PCA, Fig_PCA_loadings, labels=c("A", "B"), ncol = 1, nrow = 2, align = "hv", scale = c(1, 0.8))



message("+-------------------------------------------------------------------------------")
message("+                             DESeq2 MA plot                                    ")
message("+-------------------------------------------------------------------------------")


Title <- "MA plot"
MAplt <- ggplot(data = results.df, aes(x=baseMean, y=log2FoldChange)) + 
  geom_point(size=1, alpha=0.25, col="black") +
  geom_point(data=subset(results.df, (padj <= significance & log2FoldChange >= 1)),      size=1.25, alpha=0.5,  col="red") +
  geom_point(data=subset(results.df, ((padj <= significance) & (log2FoldChange <= -1))),     size=1, alpha=0.5,  col="blue") +
  scale_x_log10() +
  xlab("Normalised Read Count") + ylab("log2 Fold Change") + #ggtitle(paste(Project, Title, sep=" ")) +
  geom_abline(intercept = 1, slope = 0, colour='green', alpha=0.25) + 
  geom_abline(intercept = -1, slope = 0, colour='green', alpha=0.25) +
  scale_y_continuous(limits=c(-15,15), breaks=seq(-14,14,2), expand = c(0, 0))
MAplt

pdf(paste(Project, "_DESeq_MAplot", ".pdf", sep=""),width=7,height=5, onefile=FALSE)
par(bg=NA)
print({ MAplt })
dev.off()




message("+-------------------------------------------------------------------------------")
message("+                                   Heatmap                                     ")
message("+-------------------------------------------------------------------------------")

selected_gene_tbl <- read.csv("MANUSCRIPT____Bl6_TRANSFER_vs_NON-TRANSFER/Genes_for_heatmap.csv")
genes2plot <- selected_gene_tbl$ensembl_gene_id



mat          <- as.matrix(assay(rld[rownames(rld) %in% genes2plot,]))
mat.ann      <- unique(merge(mat, ensEMBL2id[,c(2,3)], by.x = "row.names", by.y= "ensembl_gene_id" ))
rownames(mat.ann) <- mat.ann$external_gene_name
matrix <- mat.ann[,-c(1,ncol(mat.ann))]
matrix          <- matrix - rowMeans(matrix)   #    MeanCentred

matrix <- matrix[match(selected_gene_tbl$external_gene_name.x, rownames(matrix)),]
matrix <- matrix[,c(5:8,1:4)]
matrix$split <- selected_gene_tbl$Gene_type
matrix$split <- factor(matrix$split, levels=c("Top up", "Top down", "Imprinted genes", "Placenta phenotype", "Prolactins", "Reg by DNA meth", "Reg histone mods"))
  
f1 = colorRamp2(seq(min(matrix[,c(1:8)]), max(matrix[,c(1:8)]), length = 3), c("green", "#000000", "red"), space = "RGB")
lgd1 = Legend(col_fun = f1, title = "Expression", at = c(seq(min(matrix[,c(1:8)]), max(matrix[,c(1:8)]), length = 3)))



ht1 = Heatmap(matrix[,c(1:8)], split = (matrix$split), col = f1, name = "non-transferred        transferred",  row_title = "", column_title = "non-transferred         transferred  ", show_row_names = TRUE, show_column_names = FALSE, heatmap_legend_param = list(title = "Expression", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = TRUE )


pdf(paste("Fig_2D",Project,  "gene_expression"  , "51", "ComplexHeatmap", ".pdf", sep="_"), onefile=FALSE, width=6, height=15) 
par(bg=NA)
draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),
     column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()







message("+-------------------------------------------------------------------------------")
message("+          Prepare results table for enrichment analysis                        ")
message("+-------------------------------------------------------------------------------")

results.df <- unique(results.df[order(results.df$padj),])
results.df$ensembl_gene_id <- results.df$gene
rownames(results.df) <- results.df$ensembl_gene_id
resdata <- merge(results.df, as.data.frame(counts(dds,normalized=TRUE)), by.x = "ensembl_gene_id", by.y='row.names', sort=FALSE)

# Annotate significant hits
resSig.df <- subset(results.df, padj < significance)
resSig.ann <- unique(merge(resSig.df, ensEMBL2id, by="ensembl_gene_id"))
resSig.ann$description <- gsub("..Source.*", "", resSig.ann$description)

resSig_l2fc1 <- subset(resSig.ann, abs(resSig.ann$log2FoldChange) > 1)
resSig_l2fc2 <- subset(resSig.ann, abs(resSig.ann$log2FoldChange) > 2)


#write.csv(resdata, paste(Project, "resdata", ".csv", sep = "_"))
#write.csv(resSig.ann, paste(Project, "resSig.ann", ".csv", sep = "_"))
#write.csv(resSig_l2fc1, paste(Project, "resSig_l2fc1", ".csv", sep = "_"))
#write.csv(resSig_l2fc2, paste(Project, "resSig_l2fc2", ".csv", sep = "_"))





message("+-------------------------------------------------------------------------------")
message("+                             KEGG pathways                                     ")
message("+-------------------------------------------------------------------------------")

Kegg_genes <- resSig_l2fc2[!is.na(resSig_l2fc2$entrezgene),]
head(Kegg_genes)
nrow(Kegg_genes)
Kegg_genes <- Kegg_genes[order(Kegg_genes$entrezgene, abs(Kegg_genes$log2FoldChange) ), ]        # sort by id and reverse of abs(value)
Kegg_genes <- Kegg_genes[!duplicated(Kegg_genes$entrezgene), ]                                   # take the first row within each id
colnames(Kegg_genes)
colnames(Kegg_genes)[colnames(Kegg_genes)=="entrezgene"] <- "Gene_ID" # exchange column name "entrez" to "Gene_ID"

Kegg_genes <- Kegg_genes[order(Kegg_genes$log2FoldChange, decreasing = T),]
foldchanges = Kegg_genes$log2FoldChange
names(foldchanges) = Kegg_genes$Gene_ID
head(foldchanges)

kk <- enrichKEGG(Kegg_genes$Gene_ID, organism="mmu", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05) #, readable=TRUE)
#head(summary(kk))
kk_results <- as.data.frame(kk)

write.csv(kk_results, paste(Project, "kk_results_enrichKegg_qval0.05_resSig_l2fc2_padj0.05", ".csv", sep = "_"))

keggresids = kk_results$ID
keggresids

# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu", new.signature=TRUE)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))

plot_pathway("mmu05332")




# KEGG analysis


# KEGG Module over-representation test
mkk <- enrichMKEGG(gene = Kegg_genes$Gene_ID, organism = 'mmu', keyType = "kegg")
head(mkk)
mkk2 <- setReadable(mkk, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
mkk2$geneID[[1]]



library(enrichplot)
cnetplot(mkk2, categorySize="pvalue", foldChange=foldchanges, colorEdge=F, node_label=T, layout = "kk")

mkk_df <- as.data.frame(mkk2)
mkk_df_entrez <- strsplit(mkk_df$geneID, "/")
write.table(as.data.frame(mkk_df), file = paste(Project, "mkk_df_kegg_modules", "_for_resSig_padj0.05_l2fc2.txt", sep = "_" ), sep = "\t", quote = F)






message("+-------------------------------------------------------------------------------")
message("+                        Cluster profiler ::: GO                                ")
message("+-------------------------------------------------------------------------------")

ego2_mf <- enrichGO(gene = unique(resSig_l2fc2$ensembl_gene_id), OrgDb = org.Mm.eg.db, keyType = 'ENSEMBL', ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = T)
ego2_bp <- enrichGO(gene = unique(resSig_l2fc2$ensembl_gene_id), OrgDb = org.Mm.eg.db,keyType = 'ENSEMBL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = T)
ego2_cc <- enrichGO(gene = unique(resSig_l2fc2$ensembl_gene_id), OrgDb = org.Mm.eg.db,keyType = 'ENSEMBL', ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = T)

ego_cc <- as.data.frame(ego2_cc)
ego_bp <- as.data.frame(ego2_bp)
ego_mf <- as.data.frame(ego2_mf)

ego_cc$GO <- "CellularComponent"
ego_bp$GO <- "BiologicalProcess"
ego_mf$GO <- "MolecularFuction"

ego2 <- rbind(ego_bp, ego_cc, ego_mf)


write.csv(ego2, paste(Project, "ego2_enrichGO_qval0.05_resSig_l2fc2_padj0.05", ".csv", sep = "_"))



#####  here i can plot tree of significant go terms (graph with nodes)! :::::

pdf(paste(Project, "ClusterProfiler", "plotGOgraph",  "l2fc2", "BP","GO_tree.pdf", sep="_"), width = 20, height = 15 )
par(bg=NA)
plotGOgraph(ego2_bp, useFullNames = TRUE ) #firstSigNodes = 10
dev.off()

pdf(paste(Project, "ClusterProfiler", "plotGOgraph", "l2fc2","MF","GO_tree.pdf", sep="_"), width = 20, height = 15 )
par(bg=NA)
plotGOgraph(ego2_mf, useFullNames = TRUE)
dev.off()

pdf(paste(Project, "ClusterProfiler", "plotGOgraph",  "l2fc2", "CC","GO_tree.pdf", sep="_"), width = 20, height = 15 )
par(bg=NA)
plotGOgraph(ego2_cc, useFullNames = TRUE)
dev.off()









#library(RMariaDB)
message("+-------------------------------------------------------------------------------")
message("+                                karyoplotes                                    ")
message("+-------------------------------------------------------------------------------")
ensembl   <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")

UP_l2fc1 <- subset(resSig_l2fc1, resSig_l2fc1$log2FoldChange >1)
DOWN_l2fc1 <- subset(resSig_l2fc1, resSig_l2fc1$log2FoldChange < -1)

txdb      <- TxDb.Mmusculus.UCSC.mm10.ensGene

all.genes <- genes(txdb)
all.transcripts <- transcripts(txdb)
#saveRDS(all.genes, "all.genes_GRanges.rds")
#all.genes <- readRDS("all.genes_GRanges.rds")


gene.UP.symbols <- UP_l2fc1$external_gene_name
genes.UP        <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'external_gene_name'),
                                   filters = 'external_gene_name', values =gene.UP.symbols, mart = ensembl))
seqlevelsStyle(genes.UP) <- "UCSC"

gene.DN.symbols <- DOWN_l2fc1$external_gene_name
genes.DN        <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'external_gene_name'),
                                   filters = 'external_gene_name', values =gene.DN.symbols, mart = ensembl))
seqlevelsStyle(genes.DN) <- "UCSC"
 
kp <- plotKaryotype(genome="mm10", plot.type=2)
kpDataBackground(kp, color = "#FFFFFF00", data.panel=2)
kp <- kpPlotDensity(kp, all.genes)
kpPlotRegions(kp, data=genes.UP, data.panel=2, col="red",  r0=-1.0, r1=0.5, lwd=2)
kpPlotRegions(kp, data=genes.DN, data.panel=2, col="blue", r0=-1.0, r1=0.5, lwd=2)


pdf(paste(Project, "_l2fc1_SigGenes_karyoplot.pdf", sep=""), width=8, height=10, onefile=FALSE)
par(bg=NA)
kp <- plotKaryotype(genome="mm10", plot.type=2)
kpDataBackground(kp, color = "#FFFFFF00", data.panel=2)
kp <- kpPlotDensity(kp, all.genes)
kpPlotRegions(kp, data=genes.UP, data.panel=2, col="red",  r0=-1.0, r1=0.5, lwd=2)
kpPlotRegions(kp, data=genes.DN, data.panel=2, col="blue", r0=-1.0, r1=0.5, lwd=2)
dev.off()






message("+-------------------------------------------------------------------------------")
message("+                     intersect resSig and EPUs for karyoplots v2               ")
message("+-------------------------------------------------------------------------------")

ensEMBL2id2_bed_mm10 <- getBM(attributes=c( 'chromosome_name', 'start_position', 'end_position','ensembl_gene_id'), mart = ensembl) 


#  546   135
BED_resSig_l2fc1_mm10 <- unique(ensEMBL2id2_bed_mm10[ensEMBL2id2_bed_mm10$ensembl_gene_id %in% resSig_l2fc1$ensembl_gene_id,])
BED_resSig_l2fc2_mm10 <- unique(ensEMBL2id2_bed_mm10[ensEMBL2id2_bed_mm10$ensembl_gene_id %in% resSig_l2fc2$ensembl_gene_id,])

#write.table(BED_resSig_l2fc1_mm10, "BED_resSig_l2fc1_mm10.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
#write.table(BED_resSig_l2fc2_mm10, "BED_resSig_l2fc2_mm10.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



resSig_l2fc0.6 <- subset(resSig.ann, abs(resSig.ann$log2FoldChange) > 0.6)
BED_resSig_l2fc0.6_mm10 <- unique(ensEMBL2id2_bed_mm10[ensEMBL2id2_bed_mm10$ensembl_gene_id %in% resSig_l2fc0.6$ensembl_gene_id,])
#write.table(BED_resSig_l2fc0.6_mm10, "BED_resSig_l2fc0.6_mm10.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


#   
#   ###   start with liftover of EPUs to mm10 :::
#   
#   Ren_supplementary_table5_chr -----> use website: https://genome.ucsc.edu/cgi-bin/hgLiftOver  ------>   hglft_to_mm10_EPUs_Ren.bed (48 EPUs failed to lift)
#   shorten the file to just 3 columns:
#   awk '{print $1"/t"$2"/t"$3}' hglft_to_mm10_EPUs_Ren.bed > hglft_to_mm10_EPUs_Ren_3COLS.bed
#   
#   ###  now intersect with DEGs :::
#   
#   bedtools intersect -wa -a BED_resSig_l2fc2_mm10_srt_chr.bed -b hglft_to_mm10_EPUs_Ren_3COLS.bed > Bedtools_intersect_mm10_l2fc2_EPU.bed
#   114 EPU
#   
#   bedtools intersect -wa -a hglft_to_mm10_EPUs_Ren_3COLS.bed -b BED_resSig_l2fc2_mm10_srt_chr.bed > Bedtools_intersect_EPU_mm10_l2fc2.bed
#   
#   bedtools intersect -c -a hglft_to_mm10_EPUs_Ren_3COLS.bed -b BED_resSig_l2fc2_mm10_srt_chr.bed > Bedtools_intersect_EPU_mm10_l2fc2_count.bed
#   bedtools intersect -c -a hglft_to_mm10_EPUs_Ren_3COLS.bed -b BED_resSig_l2fc1_mm10_srt_chr.bed > Bedtools_intersect_EPU_mm10_l2fc1_count.bed
#   
#   bedtools intersect -wao -a BED_resSig_l2fc2_mm10_srt_chr.bed -b hglft_to_mm10_EPUs_Ren_3COLS.bed > Bedtools_intersect_mm10_l2fc2_EPU_wao.bed
#   bedtools intersect -wao -a BED_resSig_l2fc1_mm10_srt_chr.bed -b hglft_to_mm10_EPUs_Ren_3COLS.bed > Bedtools_intersect_mm10_l2fc1_EPU_wao.bed
#   bedtools intersect -wao -a BED_resSig_l2fc0.6_mm10_srt_chr.bed -b hglft_to_mm10_EPUs_Ren_3COLS.bed > Bedtools_intersect_mm10_l2fc1_EPU_wao.bed
#   
#   




#EPUs_with_DEGs_l2fc2 <- read.table(file.path(paste0(Karyo.dir, "mm10/Bedtools_intersect_EPU_mm10_l2fc2_count.bed")),  header = FALSE, sep = "\t") 
#EPUs_with_DEGs_l2fc2 <- subset(EPUs_with_DEGs_l2fc2, EPUs_with_DEGs_l2fc2$V4 > 1)

EPUs_with_DEGs_l2fc1 <- read.table(file.path(paste0(Karyo.dir, "mm10/Bedtools_intersect_EPU_mm10_l2fc1_count.bed")),  header = FALSE, sep = "\t") 
EPUs_with_DEGs_l2fc1 <- subset(EPUs_with_DEGs_l2fc1, EPUs_with_DEGs_l2fc1$V4 > 1)


#DEGs_in_EPU_l2fc2 <- read.table(file.path(paste0(Karyo.dir, "mm10/Bedtools_intersect_mm10_l2fc2_EPU_wao.bed")),  header = FALSE, sep = "\t") 
DEGs_in_EPU_l2fc1 <- read.table(file.path(paste0(Karyo.dir, "mm10/Bedtools_intersect_mm10_l2fc1_EPU_wao.bed")),  header = FALSE, sep = "\t") 
#DEGs_in_EPU_l2fc2 <- subset(DEGs_in_EPU_l2fc2, DEGs_in_EPU_l2fc2$V5 != ".")
DEGs_in_EPU_l2fc1 <- subset(DEGs_in_EPU_l2fc1, DEGs_in_EPU_l2fc1$V5 != ".")

#DEGs_in_EPU_l2fc2 <- DEGs_in_EPU_l2fc2[DEGs_in_EPU_l2fc2$V6 %in% EPUs_with_DEGs_l2fc2$V2,]
DEGs_in_EPU_l2fc1 <- DEGs_in_EPU_l2fc1[DEGs_in_EPU_l2fc1$V6 %in% EPUs_with_DEGs_l2fc1$V2,]


colnames(DEGs_in_EPU_l2fc1) <- c("chr", "start", "end", "ensembl_gene_id", "chr2", "start2", "end2", "overlap_with_EPU")
#DEGs_in_EPU_l2fc2$EPU <- c("EPU_1", "EPU_1","EPU_1", "EPU_2","EPU_2","EPU_2","EPU_3","EPU_3", "EPU_4","EPU_4","EPU_5","EPU_5","EPU_6","EPU_6","EPU_6","EPU_6", "EPU_7","EPU_7")
DEGs_in_EPU_l2fc2$EPU <- c("EPU_1", "EPU_1",
                           "EPU_2", "EPU_2","EPU_2",
                           "EPU_3","EPU_3","EPU_3", "EPU_3",
                           "EPU_4","EPU_4",
                           "EPU_5","EPU_5",
                           "EPU_6","EPU_6",
                           "EPU_7", "EPU_7",
                           "EPU_8", "EPU_8","EPU_8","EPU_8",
                           "EPU_8","EPU_9",
                           "EPU_10","EPU_10")


DEGs_in_EPU_l2fc2_ann <- unique(merge(DEGs_in_EPU_l2fc2, resSig.ann[,c(1,3,7,8)], by= "ensembl_gene_id", all.x = TRUE))
#write.csv(DEGs_in_EPU_l2fc2_ann2, "DEGs_in_EPU_l2fc2_mm10.csv")


colnames(DEGs_in_EPU_l2fc1) <- c("chr", "start", "end", "ensembl_gene_id", "chr2", "start2", "end2", "overlap_with_EPU")
DEGs_in_EPU_l2fc1$EPU <- c("EPU_1", "EPU_1",
                           "EPU_2", "EPU_2","EPU_2",
                           "EPU_3","EPU_3","EPU_3", "EPU_3",
                           "EPU_4","EPU_4",
                           "EPU_5","EPU_5",
                           "EPU_6","EPU_6",
                           "EPU_7", "EPU_7",
                           "EPU_8", "EPU_8","EPU_8","EPU_8",
                           "EPU_8","EPU_9",
                           "EPU_10","EPU_10")


DEGs_in_EPU_l2fc1_ann <- unique(merge(DEGs_in_EPU_l2fc1, resSig.ann[,c(1,3,7,8)], by= "ensembl_gene_id", all.x = TRUE))
DEGs_in_EPU_l2fc1_ann <- DEGs_in_EPU_l2fc1_ann[order(DEGs_in_EPU_l2fc1_ann$start),]
DEGs_in_EPU_l2fc1_ann <- DEGs_in_EPU_l2fc1_ann[order(DEGs_in_EPU_l2fc1_ann$chr),]
#DEGs_in_EPU_l2fc1_ann <- unique(merge(DEGs_in_EPU_l2fc1_ann, ensEMBL2id2[,c(2,3)], by = "ensembl_gene_id", all.x = TRUE))
#write.csv(DEGs_in_EPU_l2fc1_ann, "DEGs_in_EPU_l2fc1_mm10.csv")
DEGs_in_EPU_l2fc1_ann <- read.csv("DEGs_in_EPU_l2fc1_mm10.csv")





#
#   Using the tools:  http://crossmap.sourceforge.net/   for a lift-over !!!
# 
# ```
# cd ~/Documents/CTR-Groups/Erica_Watson/CTR_edw23_0002/KARYOPLOTS/Karyoplots
# CrossMap.py bigwig mm10/mm9ToMm10.over.chain.gz mm9/wgEncodeLicrHistonePlacH3k04me1FAdult8wksC57bl6StdSig.bigWig mm10/mm10_wgEncode_Licr_HistonePlac_H3k04me1_FAdult8wks_C57bl6_StdSig
# CrossMap.py bigwig mm10/mm9ToMm10.over.chain.gz mm9/wgEncodeLicrHistonePlacH3k04me3FAdult8wksC57bl6StdSig.bigWig mm10/mm10_wgEncode_Licr_HistonePlac_H3k04me3_FAdult8wks_C57bl6_StdSig
# CrossMap.py bigwig mm10/mm9ToMm10.over.chain.gz mm9/wgEncodeLicrHistonePlacH3k27acFAdult8wksC57bl6StdSig.bigWig mm10/mm10_wgEncode_Licr_HistonePlac_H3k27ac_FAdult8wks_C57bl6_StdSig
# CrossMap.py bigwig mm10/mm9ToMm10.over.chain.gz mm9/wgEncodeLicrHistonePlacInputFAdult8wksC57bl6StdSig.bigWig mm10/mm10_wgEncode_Licr_HistonePlac_Input_FAdult8wks_C57bl6_StdSig
# ```
# 


message("+-------------------------------------------------------------------------------")
message("+                 Karyoplots 2: show regulatory regions                         ")
message("+-------------------------------------------------------------------------------")

# mm10 (?) The Ensembl Regulatory Build (Zerbino et al. 2015)
#    https://www.ensembl.org/info/genome/funcgen/regulatory_build.html

mart_funcgen <- biomaRt::useMart(biomart="ENSEMBL_MART_FUNCGEN", dataset = "mmusculus_regulatory_feature") # mm9
listDatasets(mart_funcgen)
listAttributes(mart_funcgen)
FUNCGEN_x <- getBM(attributes=c('chromosome_name', 'chromosome_start', 'chromosome_end', "regulatory_stable_id", "feature_type_name"), mart = mart_funcgen)          
#FUNCGEN <- getBM(attributes=c('chromosome_name', 'chromosome_start', 'chromosome_end', "activity"), mart = mart_funcgen)       
#write.table(FUNCGEN_x, "FUNCGEN_x.bed", sep = "/t", col.names = FALSE, row.names = FALSE) 




library("karyoploteR")
library("TxDb.Mmusculus.UCSC.mm10.ensGene")
library("biomaRt")
library("regioneR")
library("BSgenome")
library("BSgenome.Mmusculus.UCSC.mm10")
library("org.Mm.eg.db")
library("GenomicRanges")
library("GenomicFeatures")

#test.region    <- toGRanges("chr6",130035000,130245000, genome="mm9"); regionNum = 1; 
#test.region    <- toGRanges("chr3",92360000,92386000, genome="mm9"); regionNum = 2; 
#test.region    <- toGRanges("chr7",51040000,51080000,   genome="mm9"); regionNum = 2; gene.FCs   <- c(0.001,-4.67,2.03)
#test.region    <- toGRanges("chr1",175650000,175950000, genome="mm9"); regionNum = 3; gene.FCs   <- c(2.83,2.4)

#ensEMBL2id_mm10_entrez2symbol <- unique(ensEMBL2id2[,c(2,4)])
#head(ensEMBL2id_mm10_entrez2symbol)

head(resSig.ann)
head(ensEMBL2id)

base.url <- Karyo.dir


#histone.colours <- rev(c("purple", "red","green", "grey")) #  k4me3=red    K4me1=purp    K27ac=green    Inp=grey


test.region1    <- toGRanges("chr14",68000000,68600000, genome="mm10"); regionNum = 1; 
choose_chr <- "14"
histone.ymax   <- rev(c(10,5,5,5)) #
fc.ymax = 10
fc.ymin = -10


test.region2    <- toGRanges("chrX",169960000,170000000, genome="mm10"); regionNum = 1; 
choose_chr <- "X"
histone.ymax   <- rev(c(10,10,10,10))
fc.ymax = 5
fc.ymin = -5


test.region3    <- toGRanges("chr17",35960000,36230000, genome="mm10"); regionNum = 1; 
choose_chr <- "17"
histone.ymax   <- rev(c(20,10,10,10))
fc.ymax = 10
fc.ymin = -10


test.region4    <- toGRanges("chr3",92142000,92641000, genome="mm10"); regionNum = 1; 
choose_chr <- "3"
histone.ymax   <- rev(c(5,5,5,5))
fc.ymax = 5
fc.ymin = -5

test.region5    <- toGRanges("chr3",94587000,94821000, genome="mm10"); regionNum = 1; 
choose_chr <- "3"
histone.ymax   <- rev(c(10,10,10,10))
fc.ymax = 5
fc.ymin = -5

test.region6    <- toGRanges("chr4",120140000,120320000, genome="mm10"); regionNum = 1; 
test.region6b    <- toGRanges("chr4",119660000,120320000, genome="mm10"); regionNum = 1; 
choose_chr <- "4"
histone.ymax   <- rev(c(5,5,5,5))
fc.ymax = 5
fc.ymin = -5

	
test.region7    <- toGRanges("chr4",141360000,141450000, genome="mm10"); regionNum = 1; 
choose_chr <- "4"
histone.ymax   <- rev(c(10,10,10,10))
fc.ymax = 10
fc.ymin = -10

		
test.region8    <- toGRanges("chr6",129550000,131370000, genome="mm10"); regionNum = 1; 
choose_chr <- "6"
histone.ymax   <- rev(c(10,5,5,5))
fc.ymax = 5
fc.ymin = -5

	
test.region9    <- toGRanges("chr9",35550000,36710000, genome="mm10"); regionNum = 1; 
choose_chr <- "9"
histone.ymax   <- rev(c(10,5,5,5))
fc.ymax = 5
fc.ymin = -5

	
test.region10    <- toGRanges("chrX",169785000,169995000, genome="mm10"); regionNum = 1; 
choose_chr <- "X"
histone.ymax   <- rev(c(10,5,5,5))
fc.ymax = 5
fc.ymin = -5






test.region <- test.region6

tick.dist  <- c(100000, 10000, 100000)
mtick.dist <- c(20000, 2000, 20000)

pp                <- getDefaultPlotParams(plot.type=1)
pp$leftmargin     <- 0.15
pp$topmargin      <- 15
pp$bottommargin   <- 15
pp$ideogramheight <- 5
pp$data1inmargin  <- 10


##########
kp <- plotKaryotype(zoom = test.region, cex=2, plot.params = pp, genome="mm10")
kpAddBaseNumbers(kp, tick.dist = tick.dist[regionNum], minor.tick.dist = tick.dist[regionNum], add.units = TRUE, cex=1.3, digits = 6)
##########

genes.data                        <- makeGenesDataFromTxDb(txdb=txdb, karyoplot=kp,
                                                             plot.transcripts = TRUE, plot.transcripts.structure = TRUE)

genes.data$genes$external_gene_id <- sapply(genes.data$genes$gene_id, function(x) paste(ensEMBL2id[match(x, ensEMBL2id$ensembl_gene_id),]$external_gene_name))
#genes.data$genes$external_gene_id <- ensEMBL2id2[unique(ensEMBL2id2$ensembl_gene_id) %in% unique(genes.data$genes$gene_id),]$external_gene_name
gn2tn                             <- genes.data$genes$external_gene_id
names(gn2tn)                      <- genes.data$genes$gene_id


foldchanges_tmp <- resSig.ann[resSig.ann$ensembl_gene_id %in% genes.data$genes$gene_id,]$log2FoldChange
names(foldchanges_tmp) <- resSig.ann[resSig.ann$ensembl_gene_id %in% genes.data$genes$gene_id,]$ensembl_gene_id
genes.data$genes$log2FC <- sapply(genes.data$genes$gene_id, function(x) paste(foldchanges_tmp[match(x, names(foldchanges_tmp))]))
genes.data$genes$ensembl_gene_id <- genes.data$genes$gene_id
#genes.data$genes$entrezgene <- ensEMBL2id2[unique(ensEMBL2id2$ensembl_gene_id) %in% genes.data$genes$gene_id,]$entrezgene
#genes.data$genes$gene_id <- genes.data$genes$entrezgene


genes.data$genes$log2FC <- as.numeric(as.character(genes.data$genes$log2FC))
#genes.data$genes$log2FC[is.na(genes.data$genes$log2FC)] <- 0
gene.FCs <- as.numeric(genes.data$genes$log2FC)





##########
kpPlotGenes(kp, data=genes.data, gene.names=gn2tn, r0=0, r1=0.15, gene.name.cex = 1, plot.transcripts=TRUE, 
            add.transcript.names=FALSE, add.strand.marks=TRUE, mark.height=0.1, col="darkblue", 
            plot.transcripts.structure = FALSE, avoid.overlapping=TRUE)
##########



# CpGs

CpG.file <- paste0(base.url, "mm10/mm10_CpGIslands.bed")
CpGs     <- toGRanges(CpG.file, genome="mm10")
kpPlotRegions(kp, CpGs, border=NULL, col="darkgreen", r0=0.15, r1=0.2)
kpAddLabels(  kp, labels = "CpG Islands",             r0=0.15, r1=0.18, cex=1.4)

# EPUs
EPU.file <- paste0(base.url, "mm10/hglft_to_mm10_EPUs_Ren_3COLS.bed")
EPUs     <- toGRanges(EPU.file, genome="mm10")
kpPlotRegions(kp, EPUs, border="black", col="pink", r0=0.21, r1=0.3 )
kpAddLabels(  kp, labels = "EPUs",                  r0=0.21, r1=0.21, cex=1.4)

#logFC Plot
#fc.ymax = 5
#fc.ymin = -5
kpPoints(kp, data=genes.data$genes, y=gene.FCs, cex=1.5, ymax=fc.ymax, ymin=fc.ymin, r0=0.25, r1=0.45, col="red")
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r0=0.25, r1=0.45)
kpAddLabels(kp, labels = "log2 FC", srt=0, pos=2, r0=0.25, r1=0.45, cex=1.4, label.margin=0.035)

###reg build
#unique(FUNCGEN_x$feature_type_name)
###[1] "Promoter Flanking Region" "CTCF Binding Site"        "Enhancer"                 "Promoter"                 "Open chromatin"           "TF binding site"  
#Regulatory_build <- toGRanges(subset(FUNCGEN_x, FUNCGEN_x$chromosome_name == choose_chr), genome = "mm10")
#kpPlotRegions(kp, Regulatory_build, border="black", col=c("CTCF Binding Site" ="yellow", "Enhancer"="pink", "Promoter" = "red", "Promoter Flanking Region"= "orange", "Open chromatin" = "blue", "TF binding site" ="green"), r0=0.48, r1=0.57 )
#kpAddLabels(  kp, labels = "Reg. build Ensembl",                  r0=0.48, r1=0.50, cex=1.4)







#Histone marks

#histone.ymax   <- rev(c(10,10,10,10))
histone.colours <- rev(c("purple", "red","green", "grey")) #  k4me3 K4me1  K27ac  Inp
histone.marks <- rev(c(H3K4me1="mm10/mm10_wgEncode_H3K04me3.BigWig.bigwig",
                       H3K4me3="mm10/mm10_wgEncode_H3K04me1.BigWig.bigwig",
                       H3K27ac = "mm10/mm10_wgEncode_H3K27ac.BigWig.bigwig",
                       Input = "mm10/mm10_wgEncode_Input.BigWig.bigwig"))


total.tracks <- length(histone.marks)
out.at       <- autotrack(1:length(histone.marks), total.tracks, margin = 0.5, r0=0.50, r1=0.99)

for(i in seq_len(length(histone.marks))) {
  bigwig.file <- paste0(base.url, histone.marks[i])
  message(bigwig.file)
  at <- autotrack(i, length(histone.marks), r0=out.at$r0, r1=out.at$r1, margin = 0.3)
  kp <- kpPlotBigWig(kp, data=bigwig.file, ymax=histone.ymax[i], r0=at$r0, r1=at$r1, col = histone.colours[i], border=NULL, clipping = TRUE) #ymax="visible.region",
  
  message("calc ymax")
  computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
  message(computed.ymax)
  kpAxis(kp, ymin=0.2, ymax=computed.ymax, numticks = 2, r0=at$r0, r1=at$r1)
  kpAddLabels(kp, labels = names(histone.marks)[i], r0=at$r0, r1=at$r1, cex=1.4, label.margin = 0.035)
}






#    10  x 20
#    EPU9_kp.pdf





message("+-------------------------------------------------------------------------------")
message("+                 Plot sex-determining genes for each sample                    ")
message("+-------------------------------------------------------------------------------")

sex_genes <- c("ENSMUSG00000086503", "ENSMUSG00000069045", "ENSMUSG00000069044" ,"ENSMUSG00000069036")  
# XIST,RPS4Y1- not in mouse,DDX3Y,USP9Y,SRY

table <- resdata[resdata$ensembl_gene_id %in% sex_genes, -c(2:7)]
table <- unique(merge(table, ensEMBL2id[,c(1,2)], by = "ensembl_gene_id", all.x = TRUE))

table_t = setNames(data.frame(t(table[,-1])), table[,11])
table_t <- table_t[-c(1,nrow(table_t)),]
table_t$condition <- sample_table$Condition

table_t$condition <- as.factor(table_t$condition)
table_t$condition <- relevel(table_t$condition, "non-transferred")
table_t$sample <- rownames(table_t)
table_t <- table_t[order(table_t$condition),]

for (i in 1:4){
  table_t[[i]] <- as.numeric(as.character(table_t[[i]]))
}


table_t$Y_genes <- rowMeans(table_t[,c(1:3)])

plot_sex_XIST    <- ggplot(table_t, aes(y= Xist, x = sample, fill = condition)) + geom_col(alpha = 0.8) + scale_fill_manual(values=c( "firebrick2", "steelblue3"))  + ggtitle("XIST expression") +  ylab("Normalised gene expression") + xlab(" ") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ geom_hline(yintercept=200,  color = "grey", size=0.2)

plot_sex_Y_genes  <- ggplot(table_t, aes(y= Y_genes, x = sample, fill = condition)) + geom_col(alpha = 0.8) + scale_fill_manual(values=c( "firebrick2", "steelblue3"))  + ggtitle("Average expression of 3 Y-linked genes (Sry, Usp9y and Ddx3y)") +  ylab("Normalised gene expression") + xlab(" ") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=200,  color = "grey", size=0.2)


pdf(paste(Project, "XIST", "expression.pdf", sep="_"),width=8,height=4)
par(bg=NA)
plot_sex_XIST
dev.off()

pdf(paste(Project, "Y_genes", "expression.pdf", sep="_"),width=8,height=4)
par(bg=NA)
plot_sex_Y_genes
dev.off()

pdf(paste(Project, "plot_sex.pdf", sep="_"),width=10,height=10)
par(bg=NA)
plot_grid(plot_sex_XIST, plot_sex_Y_genes, labels=c("A", "B"), ncol = 1, nrow = 2)
dev.off()





message("+-------------------------------------------------------------------------------")
message("+                 Decon-rnaseq                         ")
message("+-------------------------------------------------------------------------------")

library(DESeq2)
library(ggplot2)
library(dplyr)
library(cowplot)
library(purrr)
library(DeconRNASeq)

NormCounts <- as.data.frame(counts(dds, normalized=T))


ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene', 'description'), mart = ensembl)          

Immune_cells_expr_matrix       <- read.csv("DeconRNAseq/srep40508-s1_immuCC_matrix_ensembl.csv")

Immune_cells_expr_matrix_anno  <- unique(merge(Immune_cells_expr_matrix, ensEMBL2id, by = "external_gene_name"))
head(Immune_cells_expr_matrix_anno)
#Immune_cells_expr_matrix_anno  <- Immune_cells_expr_matrix_anno[order(Immune_cells_expr_matrix_anno$entrezgene),]
Immune_cells_expr_matrix_anno  <- Immune_cells_expr_matrix_anno[!duplicated(Immune_cells_expr_matrix_anno$external_gene_name),]
rownames(Immune_cells_expr_matrix_anno) <- Immune_cells_expr_matrix_anno$ensembl_gene_id
Immune_cells_expr_matrix_anno  <- Immune_cells_expr_matrix_anno[,-c(1,27:29)]

NormCounts_Immune                 <- as.matrix(NormCounts[rownames(NormCounts) %in% rownames(Immune_cells_expr_matrix_anno),])
Immune_cells_expr_matrix_anno2    <- as.matrix(Immune_cells_expr_matrix_anno[rownames(Immune_cells_expr_matrix_anno) %in% rownames(NormCounts_Immune),])




library(DeconRNASeq)


Decon_results <- DeconRNASeq(datasets=as.data.frame(NormCounts_Immune), signatures=as.data.frame(Immune_cells_expr_matrix_anno2), checksig = F, known.prop = FALSE, use.scale = TRUE, fig = TRUE)

Decon_results_df <- as.data.frame(Decon_results$out.all)

rownames(Decon_results_df) <- paste0(sample_table$Condition, "_", sample_table$sampleName)
Decon_results_df$Condition <- sample_table$Condition

rowSums(Decon_results_df[,-c(26:28)])



CellTypesList <- c('Mast.Cells', 'Neutrophil.Cells', 'Eosinophil.Cells', 'B.Cells.Memory', 'B.Cells.Naive', 'Plasma.Cells', 'T.Cells.CD8.Actived', 'T.Cells.CD8.Naive', 'T.Cells.CD8.Memory', 'M0.Macrophage', 'M1.Macrophage', 'M2.Macrophage', 'Treg.Cells',  'T.Cells.CD4.Memory', 'T.Cells.CD4.Naive', 'T.Cells.CD4.Follicular', 'Th1.Cells', 'Th17.Cells', 'Th2.Cells', 'Monocyte', 'GammaDelta.T.Cells', 'NK.Resting', 'NK.Actived', 'DC.Actived', 'DC.Immature')




library(reshape)

#
# Note: I have some additional columns for annotating the plots, hence the column selection below
#
Decon_results_df.means      <- aggregate(Decon_results_df[, c(1:25)], list(Decon_results_df$Condition), mean)
Decon_results_df.means.melt <- melt(Decon_results_df.means)

Decon_results_df.vars      <- aggregate(Decon_results_df[, c(1:25)], list(Decon_results_df$Condition), var)
Decon_results_df.vars.melt <- melt(Decon_results_df.vars)

Decon_results_df.means.melt$vars     <- Decon_results_df.vars.melt$value

#
# Calculate a normalised variance
#
Decon_results_df.means.melt$normvars <- Decon_results_df.means.melt$vars/max(Decon_results_df.means.melt$vars)

#
# Tidy up the table a bit for plotting
#

colnames(Decon_results_df.means.melt) <- c("Group", "Immune.Cell", "mean", "variance", "normalisedvariance")
Decon_results_df.means.melt$Immune.Cell <- gsub("\\.", " ", Decon_results_df.means.melt$Immune.Cell)

#Decon_results_df.means.melt$Group <- as.factor(Decon_results_df.means.melt$Group)
#Decon_results_df.means.melt$Group <- factor(Decon_results_df.means.melt$Group, levels = c("normal_nerve_d_0", "implants_d_1", "implants_d_4","implants_d_7", "implants_d_14" ,"implants_d_28", "nerve_injury_d_1", "nerve_injury_d_4" , "nerve_injury_d_7",  "nerve_injury_d_14", "nerve_injury_d_28" ) )
#Decon_results_df.means.melt$Group <- gsub("d_", "day ", Decon_results_df.means.melt$Group)
#Decon_results_df.means.melt$Group <- gsub("_", " ", Decon_results_df.means.melt$Group)
#Decon_results_df.means.melt$Group <- gsub("implants", "Implants ", Decon_results_df.means.melt$Group)
#Decon_results_df.means.melt$Group <- gsub("nerve injury", "Nerve Injury ", Decon_results_df.means.melt$Group)

#
# Plot
#

#pdf(paste(Project, "_DESeq2_ImmuneCellTypeDeconvolution_bubble_grid2.pdf", sep=""), width=8,height=8, onefile=FALSE)
#png(paste(Project, "_DESeq2_ImmuneCellTypeDeconvolution_bubble_grid2.png", sep=""), width = 1500, height = 1500 )
#spar(bg=NA)
ggplot(Decon_results_df.means.melt, aes(y = Immune.Cell, x = Group)) +
  geom_point(aes(size = mean, colour = Group, alpha=(1-normalisedvariance))) + 
  xlab("") + ylab("") +
  #scale_alpha_manual(name="Test")  + 
  guides(alpha=guide_legend(title="1-(Normalised Variance)")) +
  guides(size=guide_legend(title="Mean Proportion")) +
  guides(colour=FALSE) +
  scale_color_manual(values = c("purple3", "royalblue2","royalblue2","royalblue2","royalblue2","royalblue2", "slateblue3","slateblue3","slateblue3","slateblue3","slateblue3","slateblue3", "olivedrab2", "olivedrab2","olivedrab2","olivedrab2","olivedrab2","#CCCC33" , "#99CCFF", "#339999", "olivedrab2", "#CCCC33" , "slategray") )  +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#dev.off()




sessionInfo() 

