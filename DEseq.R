## Packages##
library("DESeq2")
library(dplyr)
library('biomaRt')
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(forcats)
library("ggrepel")

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
write.csv(df, "Data_UpandDown.csv")
df <- subset(res, res$pvalue <= 0.05 & res$log2FoldChange > 2)
write.csv(df, "Data_Upregulated.csv")
df <- subset(res, res$pvalue <= 0.05 & res$log2FoldChange < -2)
write.csv(df, "Data_Downregulated.csv")

### normalization ###

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 10)
write.csv(assay(vsd), "Data_vst_norm.csv")
write.csv(assay(rld), "Data_rlog_norm.csv")

## to extract count data ##

library(dplyr)
DGE.results.sorted <- read.csv("Data_UpandDown.csv", row.names = 1)
DGEgenes <- rownames(subset(DGE.results.sorted, padj < 0.05))
head(DGE.results.sorted)
normalizeddata <- read.csv("Data_vst_norm.csv", row.names = 1)
Visualization <- normalizeddata[DGEgenes,] %>% data.frame()
write.csv(Visualization, "Data_UpandDown_Visualization.csv")

### Finding Gene Symbol ##

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
df <- read.csv("Data_UpandDown.csv")
head(df)
genes <- df$ensembl_gene_id
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values=genes, mart= mart)
G_list
A <- merge(df, G_list, by = "ensembl_gene_id")
write.csv(A, "Data_UpandDown_withsymbol.csv")

### Volcano Plot ####
library(forcats)
library("ggrepel")
data <- read.csv("Data_UpandDown.csv", header = TRUE, row.names = 1)
data <- data %>%
  mutate(Expression = case_when(log2FoldChange >= 1 & pvalue <= 0.05 ~ "up",
                               log2FoldChange <= -1 & pvalue <= 0.05 ~ "down",
                               TRUE ~ "ns"))   
data %>%
  count(Expression) %>%
  knitr::kable()

data <- data %>%
  mutate(Expression = fct_relevel(Expression, "up", "down")) 

data %>%
  distinct(Expression) %>%
  pull()

sig_il_genes <- data %>%
  filter(symbol %in% c("CP",
                               "CHRNA3",
                               "TGM1",
                               "MMP11",
                               "RUNDC3A",
                               "CASQ2",
                               "CSTA",
                               "CNN1",
                               "H19",
                               "IDO1",
                               "ACTG2",
                               "IGF2",
                               "SYNPO2",
                               "DES"
  ))

cols <- c("up" = "#bb0c00", "down" = "green", "ns" = "grey") 
sizes <- c("up" = 3, "down" = 3, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

data %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(pvalue),
             fill = Expression,
             size = Expression,
             alpha = Expression)) + 
  geom_point(shape = 21, colour = "grey") + 
  #geom_hline(yintercept = -log10(0.05),
             #linetype = "dashed") + 
  #geom_vline(xintercept = c(log2(0.5), log2(2)),
             #linetype = "dashed") +
  geom_text_repel(data = sig_il_genes,  aes(label = symbol),
                  size=3, nudge_x = 0.01, fontface='bold') +
  scale_fill_manual(values = cols) + 
  scale_size_manual(values = sizes) +
  scale_alpha_manual(values = alphas) + 
  scale_x_continuous(breaks = c(seq(-5, 5, 1)),  
                     limits = c(-2, 2)) + labs(title = "Volcano plot") + theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),    
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.background = element_blank())

## Heatmap ##

library(pheatmap)
library(dplyr)
Counts <- read.csv("Data_updown_Visualization.csv", header = TRUE, row.names = 1)
factors <- read.csv("Data.csv", header = TRUE, row.names = 1)
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
df = read.csv("Data_UpandDownwithsymbol", header=TRUE)
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


