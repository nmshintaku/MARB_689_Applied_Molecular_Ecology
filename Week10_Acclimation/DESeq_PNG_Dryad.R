## Install and load in required packages
# install DESeq2
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if (!require("DESeq2", quietly = TRUE)){
  BiocManager::install("DESeq2")
  library(DESeq2)
}

if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('gplots')) install.packages('gplots'); library('gplots')
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')
if (!require('genefilter')) BiocManager::install("genefilter"); library('genefilter')
if (!require('plotrix')) install.packages('plotrix'); library('plotrix')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('ggrepel')) install.packages('ggrepel'); library('ggrepel')
if (!require('pheatmap')) install.packages('pheatmap'); library('pheatmap')
if (!require('factoextra')) install.packages('factoextra'); library('factoextra')
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')
if (!require('tibble')) install.packages('tibble'); library('tibble')

# set working directory on your personal computer
setwd("G:/My Drive/Applied Molecular Ecology/2025/W10_Acclimation/RscriptInputFiles")

# read in the coral host counts table.
countsHost<-read.table("allcounts_HostOnly.txt",header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample - you must open R in the same folder or change the working directory
head(countsHost) 

# How many transcripts are in the coral counts table?
length(countsHost[,1]) 

# split by site
countUpa<-countsHost[,30:59]

countDobu<-countsHost[,1:29]

#######################
# Creating table of conditions for your experiment, origin and transplant site

origin=treatment<-c(1:length(names(countsHost)))
origin[grep("D",names(countsHost))]="Dobu"
origin[grep("I",names(countsHost))]="Upa"
treatment[grep("C",names(countsHost))]="control"
treatment[grep("B",names(countsHost))]="seep"

conditions=data.frame(id=colnames(countsHost),cbind(origin,treatment))
head(conditions)

conditions$origin <- as.factor(conditions$origin)
conditions$treatment <- factor(conditions$treatment, levels=c("control", "seep"))

# split by site
conUpa<-conditions[30:59,]
conDobu<-conditions[1:29,]
############################################
# load in data and QC
# read in the count and metadata table into a deseq object
dds <- DESeqDataSetFromMatrix(countData = countsHost,
                              colData = conditions,
                              design= ~ origin * treatment)

dds <- DESeqDataSetFromMatrix(countData = countUpa,
                              colData = conUpa,
                              design= ~ treatment)

dds <- DESeqDataSetFromMatrix(countData = countDobu,
                              colData = conDobu,
                              design= ~ treatment)
# How many genes are present?
dds

# let's remove low read count genes
rs=rowSds(counts(dds)) #using standard deviation as quality filtering metric based on analyses above
theta=0.6 #peak rej now occur at theta of 0.6 using SD for all..
use=(rs>quantile(rs,probs=theta)) ###
table(use) 

dds<-dds[use,]

# how many genes are left?
dds

# extract variance stabilized transformed counts
vsd <- vst(dds)

# pearson correlation between samples
vsd_mat <- assay(vsd) 
vsd_cor <- cor(vsd_mat)
head(vsd_cor) 

pheatmap(vsd_cor)

#all genes PCA
pca <- prcomp(t(assay(vsd)),center = TRUE,
              scale. = TRUE)
percentVar <- round(100 * pca$sdev^2/sum(pca$sdev^2),1)

intgroup.df <- as.data.frame(colData(dds)[, c("origin","treatment"), drop = FALSE])
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4],
                group = group, 
                intgroup.df, names = colnames(dds))

# PC 1 & 2
a<-ggplot(d, aes(PC1, PC2, color=treatment, shape=origin, fill=treatment, label=names, linetype=origin)) + 
  theme_bw()+
  geom_point(size=6) +
  stat_ellipse(type = "t")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  scale_color_manual(values=c(control="blue",seep="red"))+
  scale_fill_manual(values=c(control="white",seep="white")) +
  scale_shape_manual(values=c(19,21)) +
  geom_text_repel(point.padding	=0.3, size=4,min.segment.length = Inf, )

# PC2& 3
b<-ggplot(d, aes(PC2, PC3, color=treatment, shape=origin, fill=treatment, label=names, linetype=origin)) + 
  theme_bw()+
  geom_point(size=6) +
  stat_ellipse(type = "t")+
  xlab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance"))+
  scale_color_manual(values=c(control="blue",seep="red"))+
  scale_fill_manual(values=c(control="white",seep="white")) +
  scale_shape_manual(values=c(19,21)) +
  geom_text_repel(point.padding	=0.3, size=4,min.segment.length = Inf, )

# PC3 & 4
c<-ggplot(d, aes(PC3, PC4, color=treatment, shape=origin, fill=treatment, label=names, linetype=origin)) + 
  theme_bw()+
  geom_point(size=6) +
  stat_ellipse(type = "t")+
  xlab(paste0("PC3: ",percentVar[3],"% variance")) +
  ylab(paste0("PC4: ",percentVar[4],"% variance"))+
  scale_color_manual(values=c(control="blue",seep="red"))+
  scale_fill_manual(values=c(control="white",seep="white")) +
  scale_shape_manual(values=c(19,21)) +
  geom_text_repel(point.padding	=0.3, size=4,min.segment.length = Inf, )

ggarrange(a, b, c, labels="auto")

# put the outlier samples here that you want to remove
outliers<-c("A","B","C")

outliers<-c("DC2","IB5","IC1","IC4", "DC1","DB15","IB3","DB9","IB15")

# filter the dds object to remove those outliers before proceeding
dds<-dds[,!colnames(dds) %in% outliers]
dds

# re-run the PCA with the filtered dds object


##################################################
# Differential expression analysis on all host data in DESeq2
# calculate the differential expression
dds <- DESeq(dds)

# outlier detection of genes
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

# detect outliers
m <- nrow(attr(dds, "dispModelMatrix"))
p <- ncol(attr(dds, "dispModelMatrix"))
cooksCutoff <- qf(0.99, p, m - p)
cooksOutlier <- mcols(dds)$maxCooks > cooksCutoff

w <- which(cooksOutlier)

# that is the outliers, four as reported in the summary(res)
rownames(dds)[w]

# name of the different test
resultsNames(dds) # lists the coefficients

# results by origin
resO <- results(dds, name="origin_Upa_vs_Dobu", alpha=0.1)
summary(resO)

# results by treatment
resT <- results(dds, name="treatment_seep_vs_control", alpha=0.1)
summary(resT)

# results by interaction
resI <- results(dds, name="originUpa.treatmentseep", alpha=0.05)
summary(resI)

# PCA based on top 50 genes in coral host

# top 50 genes for PCA
topnum=50

ori<-data.frame(head(resO[order(resO$padj),],topnum)) %>% filter(padj < 0.05) %>% mutate(grp="ori")
treat<-data.frame(head(resT[order(resT$padj),],topnum)) %>% filter(padj < 0.05) %>% mutate(grp="treat")
int<-data.frame(head(resI[order(resI$padj),],topnum)) %>% filter(padj < 0.05) %>% mutate(grp="int")

vsd <- vst(dds)
rv <- assay(vsd)[rownames(ori), ]
pca <- prcomp(t(rv),center = TRUE,
              scale. = TRUE)
percentVar <- round(100 * pca$sdev^2/sum(pca$sdev^2),1)

intgroup.df <- as.data.frame(colData(dds)[, c("origin","treatment"), drop = FALSE])
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4],
                group = group, 
                intgroup.df, names = colnames(dds))

# PC 1 & 2
ggplot(d, aes(PC1, PC2, color=treatment, shape=origin, fill=treatment, label=names, linetype=origin)) + 
  theme_bw()+
  geom_point(size=6) +
  stat_ellipse(type = "t")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  scale_color_manual(values=c(control="blue",seep="red"))+
  scale_fill_manual(values=c(control="white",seep="white")) +
  scale_shape_manual(values=c(19,21)) +
  geom_text_repel(point.padding	=0.3, size=4,min.segment.length = Inf, )

##############################
#plot individual genes

l <- list()

genelist<-rownames(treat)

for(i in genelist){
  print(i)
  
  d <- plotCounts(dds, gene=i, intgroup=c("origin","treatment"), 
                  returnData=TRUE)
  d$origin <- factor(d$origin, levels = c("Dobu", "Upa"))
  d$treatment <- factor(d$treatment, levels = c("control", "seep"))
  
  p<-ggplot(d, aes(x=origin, y=count,color=treatment,fill=treatment,shape=origin)) + 
    geom_point(position=position_jitter(w=0.1,h=0),size=3) +
    theme_bw()+
    ggtitle(label=paste0(i,"; p.adj=",round(resT$padj[match(i,rownames(resT))],3))) +
    theme(text = element_text(size=10), axis.text.x = element_text(face = "italic"),
          legend.position = "none",
          plot.title = element_text(size=8, face="bold"),
          plot.subtitle = element_text(size=6))+
    xlab("")+
    ylab("Normalized Read Counts") +
    scale_shape_manual(values=c(19,21)) +
    scale_color_manual(values=c(control="blue",seep="red"))+
    scale_fill_manual(values=c(control="white",seep="white")) 
  
  name <- paste(i)
  l[[name]] <- p
  
}

# print to new windows
marrangeGrob(grobs=l,nrow=3,ncol=3)

########################## counting, venn diagram:

tab<-as.data.frame(resO) %>% rownames_to_column() %>% 
  left_join(as.data.frame(resT) %>% rownames_to_column() %>% select(rowname,log2FoldChange,padj), by=c("rowname")) %>%
  mutate(padj.ori=padj.x, padj.treat=padj.y, lFC.ori=log2FoldChange.x, lFC.treat=log2FoldChange.y) %>%
  left_join(as.data.frame(resI) %>% rownames_to_column() %>% select(rowname,log2FoldChange,padj), by=c("rowname")) %>%
  mutate(padj.int=padj, lFC.int=log2FoldChange) %>% select(-padj.x, -padj.y, -padj,-log2FoldChange.x, -log2FoldChange.y, -log2FoldChange)

inter<-tab %>% filter (tab$padj.int<=0.1)
ori<-tab %>% filter (tab$padj.ori<=0.1)
treat<-tab %>% filter (tab$padj.treat<=0.1)


candidates=list("seep"=treat$rowname,"origin"=ori$rowname,"interaction"=inter$rowname)
venn(candidates)


######################################

# Heatmaps
df <- as.data.frame(colData(dds)[,c("origin","treatment")])

# by treatment
pheatmap(assay(vsd)[treat$rowname,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df,scale="row")

# by origin
pheatmap(assay(vsd)[ori$rowname,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df,scale="row")

pheatmap(assay(vsd)[inter$rowname,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df,scale="row")

