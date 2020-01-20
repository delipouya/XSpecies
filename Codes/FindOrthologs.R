source('Codes/Functions.R')
Initialize()

# Ortholog gene selection
# To compare transcription between species, we first created a gene ortholog table using the mouse genome as 
# the reference gene list. We performed gene homology search, using ensemble multiple species comparison tool 
# (http://www.ensembl.org/biomart/ martview/42ddd77f8b0f4aae7d9eefe32cc4518c). Each species was compared to 
# mouse and a high-quality ortholog genes list was extracted (gene order conservation score above 75, whole 
# genome alignment score above 75 and minimum sequence identity above 80%). To account for gene paralogs and gene-duplication events, 
# an aggregated table of ‘‘meta-genes’’ was created. Each meta-gene may include all gene symbols homologous to one mouse gene. 
# For each organism, read counts were combined across all manifestations of each meta-gene. 
# For example, if zebrafish’s actb1 had two reads, and actb2 three reads, the actb meta-gene received a total of five reads. 
# Missing genes in species were given an ‘‘A’’ value.


geneSymbolsMapped <- readRDS('objects/5.DEgeneSymbolsMappedEntrez.rds')
lapply(geneSymbolsMapped, head)

aGeneSetTable <- geneSymbolsMapped[[1]]
geneVector <- aGeneSetTable$ensembl_gene_id
geneVector <- geneVector[!is.na(geneVector)]

#### Query the genes using biomart
## Set up Biomart and choose the dataset
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset('rnorvegicus_gene_ensembl',mart=ensembl)

## checking the appropiate filter and attribute to use
listFilters(ensembl)[grep(listFilters(ensembl)[,1], pattern = 'symbol'),] #ensembl_gene_id
## finding the right output (attributes)
x = listAttributes(ensembl)[grep(listAttributes(ensembl)[,1], pattern = 'homolog'),] 
x[grep(x$description,pattern = 'Human') ,]

## check all the attributes -> hsapiens_homolog_dn(ds) > what do these mean? 
ortho_attr <- c('hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name', 
                'hsapiens_homolog_orthology_type', 'hsapiens_homolog_perc_id', 
                'hsapiens_homolog_perc_id_r1', 'hsapiens_homolog_orthology_confidence')

geneVectorMapped <- getBM(filters="ensembl_gene_id",
                          attributes= c('ensembl_gene_id',ortho_attr),
                          values=geneVector, mart= ensembl)

geneOrthologs <- merge(aGeneSetTable, geneVectorMapped, by.x='ensembl_gene_id', by.y='ensembl_gene_id', all.x=T)
View(head(geneOrthologs,20))





## this gives an error
getBM(filters="rgd_symbol",
      attributes= c('rgd_symbol', "ensembl_gene_id",'entrezgene_id', ortho_attr),values=test, mart= ensembl)





