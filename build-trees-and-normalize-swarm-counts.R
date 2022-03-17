library("ggplot2")
library("ape")
library("vegan")

dat = read.table("~/Dropbox/PARASITES-BBW-JHV/JENNA-OCEANS/x_SWARMS-and-tax-for-anvio-counts.txt", header = TRUE, row.names = 1, sep = '\t')
reads = read.table("~/Dropbox/PARASITES-BBW-JHV/JENNA-OCEANS/x_Quality-and-primer-filtered-reads-per-sample.txt", header = TRUE, row.names = 1, sep = '\t')
tdat = t(dat)
propdat = prop.table(tdat)
tdist = vegdist(propdat, "bray")
tdclust = hclust(tdist)

dist = vegdist(dat, "bray")
dclust = hclust(dist)

write.tree(as.phylo(dclust),"x_SWARMS-and-tax-for-anvio-counts.tre")
write.tree(as.phylo(tdclust),"x_SWARMS-and-tax-for-anvio-counts-samples.tre")
