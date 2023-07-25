
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.16")

library(AnnotationHub)
library(clusterProfiler)
library(enrichplot)

## Merge the AnnotationHub dataframe with the results 
##res_ids_1 <- merge(x = res_LRT_tb, y = annot, by.x = "gene", by.y = "ensembl_gene_id")
##res_ids_1 <- subset(res_ids_1, select = -external_gene_name)
##res_ids_2 <- merge(x = res_LRT_tb, y = annot, by.x = "gene", by.y = "external_gene_name")
##res_ids_2 <- subset(res_ids_2, select = -ensembl_gene_id)
##res_ids <- rbind(res_ids_1, res_ids_2)
# Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
allOE_genes <- as.character(res_ids$entrezgene_id) 
# Extract significant results
sigOE <- dplyr::filter(res_ids, padj < 0.05)

sigOE_genes <- as.character(sigOE$entrezgene_id)

hub <- AnnotationHub()
honeybee_ens <- query(hub, "Apis mellifera")
apis_ens <- honeybee_ens[["AH109223"]]


allOE_genes <- readLines("allOE_genes.txt")
sigOE_genes <- readLines("sigOE_genes.txt")
##Run GO enrichment analysis
ego <- enrichGO(gene=sigOE_genes,
                universe=allOE_genes,
                keyType="ENTREZID",
                OrgDb = apis_ens,
                ont="BP",
                pAdjustMethod="BH",
                qvalueCutoff=0.05,
                readable=TRUE)


dotplot(ego, showCategory=50)

output_file <- "enrichGO_dotplot.pdf"
pdf(width=8, height=8, file=output_file)
dotplot(ego, showCategory=50)
dev.off()






