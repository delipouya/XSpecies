source('Codes/Functions.R')
Initialize()

library("SoupX")
library('multtest')
library("celda")
dataDirs = c("~/Desktop/SoupX_rat_DA/")
sc = load10X(dataDirs, keepDroplets = TRUE, cellIDs = rownames(pr_data_DA))
sc = estimateSoup(sc)


#get a rough sense of whether the expression of a 
#gene (or group of genes) in a set of cells is derived from the soup or not

processed_data <- readRDS('objects/rat_DA/2.seur_dimRed_rat_DA_mito_30_lib_1500.rds')
pr_data_DA <- as.data.frame(getEmb(processed_data, 'tsne'))
head(pr_data_DA)
length(sc$toc["Hbb",])
pr_data_DA$Hbb = sc$toc["Hbb",colnames(sc$toc) %in% rownames(pr_data_DA)]
gg = ggplot(pr_data_DA, aes(tSNE_1, tSNE_2, color=Hbb>0)) + geom_point(alpha=0.5)+theme_bw()
plot(gg)

## saveing the dimension reduction in the meta data
sc = setDR(sc, pr_data_DA)

#  checking if the expression of IGKC in these scattered cells is more than we would 
# expect by chance from the soup. To really answer this properly, we need to know how much 
# contamination is present in each cell, which will be the focus of the next sections.


# But we can get a rough idea just by calculating how many counts we would 
# expect for IGKC in each cell, by assuming that cell contained nothing but soup.
# The function soupMarkerMap allows you to visualise the ratio of observed counts 
# for a gene (or set of genes) to this expectation value. Let's try it out,

?soupMarkerMap
gg = plotMarkerMap(sc, "Hbb", pr_data_DA)
plot(gg)





