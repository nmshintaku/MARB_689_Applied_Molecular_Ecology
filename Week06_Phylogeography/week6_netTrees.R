.libPaths("/scratch/group/kitchen-group/Rlibs")
library(pegas)
library(phangorn)
#library(dplyr)
library(tidyverse)
library(vcfR)
library(adegenet)
library("RColorBrewer")
library(phytools)
library(psych)
library(vegan) 

# get color palette from paper
dk<- brewer.pal(n = 6, name = "Dark2")

# match assignments as the paper
dk2<-c("#7570B3","#666666","#D95F02","#E6AB02","#1B9E77","black")


######################
## Import COI  Data ##
######################

# import seq. alignment as a DNAbin object
mydna <- read.dna("/scratch/group/kitchen-group/MARB_689_Molecular_Ecology/06_phylogeo/EH_aligned_cox1.fa", format = "fasta")
mydna

# convert to read as an alignment
myalign <- as.alignment(mydna)

# view alignment
alview(mydna)

# identify nucleotide sites with polymorphisms from the alignment
seg.sites(mydna)

# read in metadata
pop<- read.table("/scratch/group/kitchen-group/MARB_689_Molecular_Ecology/06_phylogeo/EH_meta.txt", sep="\t", header=T)

##########################
# Haplotype Network- COI #
##########################

# haplotype assignment
hap <- haplotype(mydna)
hap

#calculate pairwise Hamming distance
d <- dist.dna(hap, "N")

# size of haplotype groups
sz <- summary(hap)

# extract individual's haplotype assignment
hapInfo <- utils::stack(setNames(attr(hap,"index"),rownames(hap)))
head(hapInfo)

# rename columns of table
names(hapInfo) <- c("index","haplotype")

# add names back to the table with group assignment from Bendif et al.
names <- as.data.frame(labels(mydna))
# split seq names into sample name and group assignment
meta<-names %>%  separate(`labels(mydna)`, into = c("sample", "grp"), sep="_(?=[^_]+$)") %>%
  separate(sample, into = c("accession", "strain"), sep="_(?=[^_]+$)")
head(meta)
# make into one table
merged <- data.frame(cbind(hapInfo,meta[hapInfo$index,]))
head(merged)

# let's make the minimum spanning network haplotype network 
net <- msn(d)
plot(net,size=sz)

# make table input for pie-chart
pie <- table(merged$haplotype,merged$grp)
head(pie)

# plot haplotype network with colored pie-chart
#plot(net,size=attr(net,"freq"),pie=pie)
plot(net,pie=pie, size=sz, bg=dk2,  fg="white")
legend("bottomleft", colnames(pie), col=dk2, pch=19, ncol=2)

###########
# NJ tree #
###########
# plot both NJ trees together
par(mfrow = c(1, 1),oma = c(3,3,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

# function to estimate NJ tree
estimate_tr <- function(m) nj(dist(m, method="manhattan"))

# tree estimate
tre <- estimate_tr(dist(mydna,method="manhattan"))

# 1000 bootstrap replicates
bstrees <- boot.phylo(tre, mydna, estimate_tr, trees = TRUE,B=1000)$trees

# bootstrap support/node labels
clad <- (prop.clades(tre, bstrees, rooted = FALSE)/1000)*100

# consensus tree
con <- consensus(bstrees, p=0.5)
plot(con)

# plot tree with node labels
tre$tip.label <- meta$strain
plot(ladderize(tre),  cex=1.5, main="cox1")
nodelabels(clad, bg="white", frame="none", col="blue",cex=1,adj = 1.2)

#######################
## Import tufA  Data ##
#######################

# import seq. alignment as a DNAbin object
mytuf <- read.dna("/scratch/group/kitchen-group/MARB_689_Molecular_Ecology/06_phylogeo/EH_aligned_tufA.fa", format = "fasta")
mytuf

# convert to read as an alignment
myalign <- as.alignment(mytuf)

# view alignment
alview(mytuf)

# identify nucleotide sites with polymorphisms from the alignment
seg.sites(mytuf)

###########################
# Haplotype Network- tufA #
###########################

# haplotype assignment
haptuf <- haplotype(mytuf)
haptuf

#calculate pairwise Hamming distance
dtuf <- dist.dna(haptuf, "N")

# size of haplotype groups
sztuf <- summary(haptuf)

# extract individual's haplotype assignment
haptufInfo <- stack(setNames(attr(haptuf,"index"),rownames(haptuf)))
head(haptufInfo)

# rename columns of table
names(haptufInfo) <- c("index","haplotype")

# make into one table
mergedtuf <- data.frame(cbind(haptufInfo,meta[haptufInfo$index,]))
head(mergedtuf)

# let's make the haplotype network 
nettuf <- msn(dtuf)
plot(nettuf,size=sztuf)

# make table input for pie-chart
pietuf <- table(mergedtuf$haplotype,mergedtuf$grp)
head(pietuf)

# plot haplotype network with colored pie-chart
plot(nettuf,pie=pie, size=sztuf, bg=dk2,  fg="white")
legend("topright", colnames(pietuf), col=dk2, pch=19, ncol=2)

###########
# NJ tree #
###########

# tree estimate
tretuf <- estimate_tr(dist(mytuf,method="manhattan"))

# 1000 bootstrap replicates
bstreestuf <- boot.phylo(tretuf, mytuf, estimate_tr, trees = TRUE,B=1000)$trees

# bootstrap support/node labels
cladtuf <- (prop.clades(tretuf, bstreestuf, rooted = FALSE)/1000)*100

# consensus tree
contuf <- consensus(bstreestuf, p=0.5)
plot(contuf)

# plot tree with node labels
tretuf$tip.label <- meta$strain
plot(ladderize(tretuf),  cex=1.5, main="tufA")
nodelabels(cladtuf, bg="white", frame="none", col="blue",cex=1,adj = 1.2)


# plot together
obj<-cophylo(ladderize(tre), ladderize(tretuf), fsize=2)
plot(obj)
nodelabels.cophylo(clad,frame="none", col="blue",cex=1,adj = 1.2)
nodelabels.cophylo(cladtuf,which="right",frame="none", col="blue",cex=1,adj = -0.2)
