#### in this script we'll find the pig's mitochondrial gene 

library(rtracklayer)
Sus_scrofa.Sscrofa11 <- rtracklayer::import('~/sus_scrofa/Sus_scrofa.Sscrofa11.1.99.gtf')
table(Sus_scrofa.Sscrofa11$gene_biotype)
pig_mit_genes_ensembl <- Sus_scrofa.Sscrofa11[Sus_scrofa.Sscrofa11@seqnames =='MT',]$gene_id
saveRDS(pig_mit_genes_ensembl, 'Data/pig_mit_genes_ensembl.rds')
