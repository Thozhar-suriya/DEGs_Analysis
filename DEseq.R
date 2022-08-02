## Packages##
library("DESeq2")
library(dplyr)
library('biomaRt')
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

## Uploading Data ##

data <- read.csv("Count_data.csv", header = T, row.names = 1)
meta <- read.csv("Meta_data.csv", row.names = 1)
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = meta,
                              design = ~ Tissue)
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

 ### factor analysis ####

dds$Tissue <- factor(dds$Tissue, levels = c("Tumor","Normal"))
dds$Tissue <- relevel(dds$Tissue, ref = "Normal")
dds$Tissue <- droplevels(dds$Tissue)

### Deseq ###

dds <- DESeq(dds)
res <- results(dds)
res
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
df <- subset(res, res$pvalue <= 0.05 & (res$log2FoldChange > 2 | res$log2FoldChange < -2))
write.csv(df, "PRJNA435914_Korean_Male_UpandDown.csv")
df <- subset(res, res$pvalue <= 0.05 & res$log2FoldChange > 2)
write.csv(df, "PRJNA435914_Korean_Male_Upregulated.csv")
df <- subset(res, res$pvalue <= 0.05 & res$log2FoldChange < -2)
write.csv(df, "PRJNA435914_Korean_Male_Downregulated.csv")

### normalization ###

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 10)
write.csv(assay(vsd), "PRJNA435914_Korean_Male_vst_norm.csv")
write.csv(assay(rld), "PRJNA435914_Korean_Male_rlog_norm.csv")

## to extract count data ##

library(dplyr)
DGE.results.sorted <- read.csv("PRJNA435914_Korean_Male_UpandDown.csv", row.names = 1)
DGEgenes <- rownames(subset(DGE.results.sorted, padj < 0.05))
head(DGE.results.sorted)
normalizeddata <- read.csv("PRJNA435914_Korean_Male_vst_norm.csv", row.names = 1)
Visualization <- normalizeddata[DGEgenes,] %>% data.frame()
write.csv(Visualization, "PRJNA435914_Korean_Male_UpandDown_Visualization.csv")

### Finding Gene Symbol ##

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
df <- read.csv("PRJNA435914_Korean_Male_UpandDown.csv")
head(df)
genes <- df$ensembl_gene_id
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values=genes, mart= mart)
G_list
A <- merge(df, G_list, by = "ensembl_gene_id")
write.csv(A, "PRJNA435914_Korean_Male_UpandDown_withsymbol.csv")

### Volcano Plot ####
library(EnhancedVolcano)
keyvals <- ifelse(
res$log2FoldChange < -1.5, 'green',
ifelse(res$log2FoldChange > 1.5, 'red',
'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == 'red'] <- 'High'
names(keyvals)[keyvals == 'gray'] <- 'Normal'
names(keyvals)[keyvals == 'green'] <- 'Low'
EnhancedVolcano(res,
lab = rownames(res),
x = 'log2FoldChange',
y = 'pvalue',
selectLab = rownames(res)[which(names(keyvals) %in% c('high', 'low'))],
axisLabSize = 12,
ylim = c(0, 40),
xlim = c(-6, 6),
cutoffLineType = 'blank',
transcriptPointSize = 1.5,
transcriptLabSize = FALSE,
caption = NULL,
captionLabSize = 12,
title = 'Volcano plot',
subtitle = 'Differential Expression (PRJNA435914)',
colCustom = keyvals,
legend = c("NA", "Up-regulated", "Down-regulated"),
legendPosition = "top",
hline = NULL,
vline = NULL,
pCutoff = 0.05,
FCcutoff = 2,
gridlines.major = FALSE,
gridlines.minor = FALSE)

## Heatmap ##

library(pheatmap)
library(dplyr)
Counts <- read.csv("PRJNA435914_Korean_updown_Visualization.csv", header = TRUE, row.names = 1)
factors <- read.csv("PRJNA435914_Korean.csv", header = TRUE, row.names = 1)
factorsDS <- dplyr::select(factors, sex, Tissue)
factorsDS$sex <- factor(factorsDS$sex, levels = c("male", "female"))
SexCol <- c("darkorchid", "darkorange")
names(SexCol) <- levels(factorsDS$sex)
factorsDS$Tissue <- factor(factorsDS$Tissue, levels = c("Normal", "Tumour"))
TissueCol <- c("forestgreen", "blue1")
names(TissueCol) <- levels(factorsDS$Tissue)
AnnColour <- list(
sex = SexCol,
Tissue = TissueCol)
AnnColour
pheatmap(Counts, scale = "row", clustering_distance_rows = "manhattan", cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = 'average', annotation_names_row = FALSE, annotation_names_col = TRUE, annotation_colors = AnnColour, annotation_col = factorsDS, show_rownames = F, show_colnames = T, fontsize = 6)

## GO Analysis ##
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

##Data##
df = read.csv("TCGA_STAD_UpandDownwithsymbol", header=TRUE)
original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$hgnc_symbol
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

##GO##
gse1 <- gseGO(geneList=gene_list, ont ="ALL", keyType = "SYMBOL", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, OrgDb = organism, pAdjustMethod = "none")
require(DOSE)
dotplot(gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
dotplot(gse, showCategory=5, split=".sign") + scale_color_gradient(low = "#30F621", high = "#EC0712")
write.csv(gse1, "GO.csv")

##MF##
gse2 <- gseGO(geneList=gene_list, ont ="MF", keyType = "SYMBOL", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, OrgDb = organism, pAdjustMethod = "none")
require(DOSE)
dotplot(gse1, showCategory=5, split=".sign") + facet_grid(.~.sign)
dotplot(gse1, showCategory=5, split=".sign") + scale_color_gradient(low = "#30F621", high = "#EC0712")
write.csv(gse2, "MF.csv")

##CC##
gse3 <- gseGO(geneList=gene_list, ont ="CC", keyType = "SYMBOL", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, OrgDb = organism, pAdjustMethod = "none")
require(DOSE)
dotplot(gse2, showCategory=5, split=".sign") + facet_grid(.~.sign)
dotplot(gse2, showCategory=5, split=".sign") + scale_color_gradient(low = "#30F621", high = "#EC0712")
write.csv(gse3, "CC.csv")


