.libPaths("/scratch/group/kitchen-group/Rlibs")

# install new packages
#library("BiocManager")
#BiocManager::install("phyloseq", lib="/scratch/group/kitchen-group/Rlibs")
#BiocManager::install("decontam", lib="/scratch/group/kitchen-group/Rlibs")

library("phyloseq")
#library("decontam")
library(reshape2)
library(vegan)
library(ggplot2)
library(microshades)

# set your working directory, change "kitchens" to your specific directory name
setwd("/scratch/group/kitchen-group/MARB_689_Molecular_Ecology/class_working_directories/kitchens/12_amplicon")

# import your files
table <- read.table(file = "feature-table.txt", sep = "\t", header = T, row.names = 1, skip = 1, comment.char = "")
OTU = otu_table(table,taxa_are_rows=TRUE)
metadata <- read.table("../../../12_amplicons/00-ANALYSIS/metadata.txt",sep = "\t", header = T,row.names=1)
SAMPLE <- sample_data(metadata)
taxonomy <- read.table(file = "../../../12_amplicons/00-ANALYSIS/taxonomy.tab", sep = "\t", header = TRUE ,row.names = 1)
TAX = tax_table(as.matrix(taxonomy))
TREE <- read_tree("tree.nwk")

# create phyloseq object
physeq <- phyloseq(OTU, TAX, SAMPLE,TREE)
physeq

# plot alpha diversity
plot_richness(physeq, x="Location", measures=c("Observed", "Shannon","Simpson")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

# create an ordination plot
ord <- ordinate(physeq,"PCoA","unifrac")
plot_ordination(physeq,ord,type="samples",color="Location")

# taxa bar plot
# Use microshades function prep_mdf to agglomerate, normalize, and melt the phyloseq object
mdf_prep <- prep_mdf(physeq, subgroup_level = "Class")

# Create a color object for the specified data
p<-plot_bar(physeq, fill = "Class", facet_grid="~Location")
p=p+facet_grid(~Location, scale="free_x", drop=TRUE)
p

#color_obj <- create_color_dfs(mdf_prep, group_level = "Phylum", subgroup_level = "Class", cvd = F)
#mdf_v1 <- color_obj$mdf
#cdf_v1 <- color_obj$cdf
#plot_1 <- plot_microshades(mdf_v1, cdf_v1, group_label = "Phylum Class")