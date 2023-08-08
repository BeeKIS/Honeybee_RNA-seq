library(DESeq2)
library(tidyverse)
library(airway)
library(ggrepel)
library(EnhancedVolcano)
library(RColorBrewer)
library(pheatmap)
library(seriation)
library(dendextend)
library(DEGreport)
library(DOSE)
library(pathview)
library(clusterProfiler)
#Preparing cont data
##Read in count data
counts_data_0 <- read.csv("counts_data.csv")
counts_data <- counts_data_0[,-1]
rownames(counts_data) <- counts_data_0[,1]
##Read in sample info
colData_0 <- read.csv("sample_info.csv")
colData <- colData_0[,-1]
rownames(colData) <- colData_0[,1]
colData_0 <- colData
colData <- data.frame(colData[,2])
colnames(colData) <- "behaviour"
colData$behaviour <- as.factor(colData$behaviour)
rownames(colData) <- rownames(colData_0)
#Construct a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData=counts_data, 
                              colData=colData,
                              design = ~ behaviour)
#Pre-filtering
##Removing rows with low gene counts
##Keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]
#Run DESeq
dds <- DESeq(dds)
res <- results(dds)
#Exploring results
##Hyg VS NHyg
res_Hyg_vs_NHyg <- results(dds, contrast=c("behaviour", "hygienic", "non_hygienic"), alpha=0.05)
##Clo VS NHyg
res_Clo_vs_NHyg <- results(dds, contrast=c("behaviour", "closing", "non_hygienic"), alpha=0.05)
##Hyg VS Clo
res_Hyg_vs_Clo <- results(dds, contrast=c("behaviour", "hygienic", "closing"), alpha=0.05)
#PCA
rld <- rlog(dds, blind=TRUE)
output_file <- "PCA.pdf"
pdf(width = 8, height = 6, file = output_file)
plotPCA(rld, intgroup="behaviour")
dev.off()
#Correlation
##Etract the rlog matrix from the object
rld_mat <- assay(rld)   
##Compute pairwise correlation values
rld_cor <- cor(rld_mat) 
output_file <- "correlation.pdf"
pdf(width = 8, height = 6, file = output_file)
pheatmap(rld_cor, annotation = colData)
dev.off()
#Extracting significantly DEGs
##Hyg VS NHyg
res_Hyg_vs_NHyg <- res_Hyg_vs_NHyg%>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
sig_res_Hyg_vs_NHyg  <- res_Hyg_vs_NHyg  %>%
  filter(padj < padj.cutoff)
##Clo VS NHyg
padj.cutoff <- 0.05
res_Clo_vs_NHyg <- res_Clo_vs_NHyg %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
sig_res_Clo_vs_NHyg <- res_Clo_vs_NHyg %>%
  filter(padj < padj.cutoff)
##Hyg VS Clo
res_Hyg_vs_Clo <- res_Hyg_vs_Clo %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
sig_res_Hyg_vs_Clo  <- res_Hyg_vs_Clo  %>%
  filter(padj < padj.cutoff)
#Volcano plots
##Clo VS NHyg
res_Clo_vs_NHyg <- res_Clo_vs_NHyg %>% 
  mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 1.0)

output_file <- "Clo_vs_NHyg_vulcano_enhanced.pdf"
pdf(width = 8, height = 6, file = output_file)
keyvals <- ifelse(
  res_Clo_vs_NHyg_1$log2FoldChange < -1 & res_Clo_vs_NHyg_1$padj < 0.05, '#357EBDFF',
  ifelse(res_Clo_vs_NHyg_1$log2FoldChange > 1 & res_Clo_vs_NHyg_1$padj < 0.05, '#D43F3AFF',
         '#B8B8B8FF'))
keyvals[is.na(keyvals)] <- '#B8B8B8FF'
names(keyvals)[keyvals == '#D43F3AFF'] <- 'Up'
names(keyvals)[keyvals == '#B8B8B8FF'] <- 'Not Sig'
names(keyvals)[keyvals == '#357EBDFF'] <- 'Down'
EnhancedVolcano(res_Clo_vs_NHyg,
                  lab = res_Clo_vs_NHyg$gene,
                  x = 'log2FoldChange',
                  xlim = c(-6.5, 6.5),
                  y = 'padj',
                  ylim = c(0, 9),
                  title = 'Clo vs NHyg',
                  subtitle = NULL,
                  selectLab = c('', ''),
                  xlab = bquote(~Log[2]~ 'fold change'),
                  FCcutoff = 1.0,
                  pCutoff = 0.05,
                  pointSize = 2.0,
                  labSize = 4.0,
                  labCol = 'black',
                  labFace = 'bold',
                  boxedLabels = TRUE,
                  colCustom = keyvals,
                  colAlpha = 4/5,
                  legendPosition = 'none',
                  legendLabSize = 12,
                  legendIconSize = 3.0,
                  drawConnectors = TRUE,
                  widthConnectors = 1.0,
                  colConnectors = 'black',
                  caption = NULL
)
dev.off()
##Hyg VS NHyg
res_Hyg_vs_NHyg <- res_Hyg_vs_NHyg %>% 
  mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 1.0)

output_file <- "Hyg_vs_NHyg_vulcano_enhanced.pdf"
pdf(width = 8, height = 6, file = output_file)
EnhancedVolcano(res_Hyg_vs_NHyg_1,
                  lab = res_Hyg_vs_NHyg_1$gene,
                  x = 'log2FoldChange',
                  xlim = c(-6.5, 6.5),
                  y = 'padj',
                  ylim = c(0, 9),
                  title = 'Hyg vs NHyg',
                  subtitle=NULL,
                  selectLab = c('',''),
                  xlab = bquote(~Log[2]~ 'fold change'),
                  FCcutoff = 1.0,
                  pCutoff=0.05,
                  pointSize = 2.0,
                  labSize = 4.0,
                  labCol = 'black',
                  labFace = 'bold',
                  boxedLabels = TRUE,
                  colCustom = keyvals,
                  colAlpha = 4/5,
                  legendPosition = 'none',
                  legendLabSize = 12,
                  legendIconSize = 3.0,
                  drawConnectors = TRUE,
                  widthConnectors = 1.0,
                  colConnectors = 'black',
                  caption = NULL)
dev.off()
##Hyg vs Clo
res_Hyg_vs_Clo <- res_Hyg_vs_Clo %>% 
  mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 1.0)

output_file <- "Hyg_vs_Clo_vulcano_enhanced.pdf"
pdf(width = 8, height = 6, file = output_file)
EnhancedVolcano(res_Hyg_vs_Clo_1,
                  lab = res_Hyg_vs_Clo_1$gene,
                  x = 'log2FoldChange',
                  xlim = c(-6.5, 6.5),
                  y = 'padj',
                  ylim = c(0, 9),
                  title = 'Hyg vs Clo',
                  subtitle=NULL,
                  selectLab = c('',''),
                  xlab = bquote(~Log[2]~ 'fold change'),
                  FCcutoff = 1.0,
                  pCutoff=0.05,
                  pointSize = 2.0,
                  labSize = 4.0,
                  labCol = 'black',
                  labFace = 'bold',
                  boxedLabels = TRUE,
                  colCustom = keyvals,
                  colAlpha = 4/5,
                  legendPosition = 'right',
                  legendLabSize = 14,
                  legendIconSize = 3,
                  drawConnectors = TRUE,
                  widthConnectors = 1.0,
                  colConnectors = 'black',
                  caption = NULL)
dev.off()
#Heatmap
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts <- as.data.frame(normalized_counts)
normalized_counts <- tibble::rownames_to_column(normalized_counts, "gene")
normalized_counts <- normalized_counts %>%
  as_tibble()
norm_OEsig <- normalized_counts[,c(2:5, 6:10, 11:15)] %>% 
  filter(normalized_counts$gene %in% c(sig_res_Clo_vs_NHyg$gene, sig_res_Hyg_vs_Clo$gene, sig_res_Hyg_vs_NHyg$gene)) 

output_file <- "Heatmap.pdf"
pdf(width = 8, height = 6, file = output_file)
heat_colors <- rev(brewer.pal(6, "RdBu"))
ann_colors = list(behaviour = c(closing = "#EEA236FF", hygienic = "#5CB85CFF", non_hygienic = "#9632B8FF"))
pheatmap(norm_OEsig, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation_names_col  = F ,
         annotation = colData, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10,
         annotation_colors = ann_colors,
         height = 20)
dev.off()
#DGE analysis using LRT (likelihood ratio test)
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)
res_LRT <- results(dds_lrt)
res_LRT_tb <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
sigLRT_genes <- res_LRT_tb %>% 
  filter(padj < padj.cutoff)
clustering_sig_genes <- sigLRT_genes %>%
  arrange(padj) 
cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
##Show gene clusters across sample groups
clusters <- degPatterns(cluster_rlog, metadata = colData, time = "behaviour", col=NULL)

output_file <- "degPlotCluster.pdf"
pdf(width = 8, height = 6, file = output_file)
degPlotCluster(clusters$normalized , time="behaviour", col="colored")+
  theme_minimal()
dev.off()

groups <- clusters$df
