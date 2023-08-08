library(AnnotationHub)
library(ensembldb)
library(biomaRt)
#define gene ID of interest
gene_id <- rownames(counts_data)
mart <- useMart('metazoa_mart', host = 'https://metazoa.ensembl.org') 
mart <- useDataset('amellifera_eg_gene', mart)
#get annotation
annot <- getBM(
  mart = mart,
  attributes = c(
    'entrezgene_id',
    'ensembl_gene_id',
    'external_gene_name',
    'gene_biotype',
    'go_id',
    'name_1006', 
    'geneid', 
    'namespace_1003'),
  uniqueRows = TRUE, values = gene_id)

##combine annotations with DEG data
colnames(groups) <- c("ensembl_gene_id", "cluster")
groups_external <- groups
colnames(groups_external) <- c("external_gene_name", "cluster")
groups_external$external_gene_name <- gsub("\\.", "-", groups_external$external_gene_name)
annotations_groups_external <- merge(x=annot,y=groups_external,by="external_gene_name")
annotations_groups <- merge(x=annot,y=groups,by="ensembl_gene_id")
annotations_all <- rbind(annotations_groups, annotations_groups_external)




















