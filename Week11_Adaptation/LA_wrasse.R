.libPaths("/scratch/group/kitchen-group/Rlibs")
library(tidyverse)
library(vcfR)
library(adegenet)
library(poppr)
library(reshape2)
library(ggplot2)
library(vegan) 
library(sp)
library(raster)
library(qvalue)

if(!requireNamespace("lfmm", quietly = TRUE)) {  
  remotes::install_github("bcm-uga/lfmm", lib="/scratch/group/kitchen-group/Rlibs")
}
library(lfmm)
#install.packages("pcadapt", lib="/scratch/group/kitchen-group/Rlibs/")
library(pcadapt)

#if(!requireNamespace("qvalue", quietly = TRUE)) {  
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#  BiocManager::install(version = "3.16", lib="/scratch/group/kitchen-group/Rlibs")
#  BiocManager::install("qvalue", lib="/scratch/group/kitchen-group/Rlibs")
#}

#################
## Import Data ##
#################

#Read in genome VCF file 
vcf_gen <- read.vcfR("/scratch/group/kitchen-group/MARB_689_Molecular_Ecology/07_hybrids/west.filt.maf0.01.recode.vcf")

# convert vcf to genind file format
gind <- vcfR2genind(vcf_gen)
gind

# loci names
loci<-locNames(gind)

#add population or species information to the genind pop slot
poptab_gen<-read.table("/scratch/group/kitchen-group/MARB_689_Molecular_Ecology/07_hybrids/Sampleinfo_metadata.txt",
                       check.names=FALSE, header=T, na.strings = c("", "NA"),
                       stringsAsFactors = FALSE, sep="\t")

head(poptab_gen)

# add region to the genlight strata
gind@pop <- as.factor(poptab_gen$Sample_ID)
strata(gind) <- data.frame(poptab_gen$`age(years)`,poptab_gen$Sex,poptab_gen$`total_length(mm)`,poptab_gen$Sample_ID)

#########################
# download climate data #
#########################

# download world climate bio data to your directory
file.copy("/scratch/group/kitchen-group/MARB_689_Molecular_Ecology/07_hybrids/wc10", ".",recursive=TRUE)
r <- getData("worldclim",var="bio",res=10, download=F)

# extract out 1=Annual Mean Temperature, 5=Max Temperature of Warmest Month
# 12=Annual Precipitation, and 15=Precipitation Seasonality
r <- r[[c(1,12,15)]]
names(r) <- c("Temp","Prec","Season")

# extract longitude and latitude from the pop-gen file
#Note - the order is LONGITUDE then LATITUDE
coords <- data.frame(poptab_gen[,7:6])

# create points to extract
points <- SpatialPoints(coords, proj4string = r@crs)

values <- extract(r,points)

# data frame with coordinates and env. variables
df <- cbind.data.frame(ID=poptab_gen$Individual_ID,coordinates(points),values)

na_DF <- df[!is.na(df$Temp),]
dim(na_DF)
str(na_DF)

#write.table(na_DF,"climateData.txt",sep="\t", row.names=F)

###########
# PCAdapt #
###########
# load in file
path_to_file <- "/scratch/group/kitchen-group/MARB_689_Molecular_Ecology/07_hybrids/west.filt.maf0.01.recode.vcf"
filename <- read.pcadapt(path_to_file, type = "vcf")

# number of PC to retain to describe the genetic structure of the data
x <- pcadapt(input = filename, K = 20)

# scree plot
plot(x, option = "screeplot", K=5)

# PC 1 & 2
plot(x, option = "scores", pop = poptab_gen$Sample_ID)

# PC 2 & 3
plot(x, option = "scores", i = 2, j = 3, pop = poptab_gen$Sample_ID)

# run pcadapt on chosen number of retained PCs
x <- pcadapt(filename, K = 3)
summary(x)

# SNPs with q-values less than α (expected FDR) are considered outliers.
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliersPC <- which(qval < alpha)
length(outliersPC)

# extract names of loci
outliersName <-loci[][outliersPC]

# which PCs are correlated with the outlier SNPs?
snp_pc <- get.pc(x,outliers)
head(snp_pc, 10)


#######################
# Redundancy Analysis #
#######################
# impute missing data
gen.imp <- apply(gind@tab[rownames(gind@tab) %in% na_DF$ID,], 2, 
                 function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(gen.imp)) # No NAs
dim(gen.imp)

# reduce climate table to just env predictors
pred <- na_DF[,c(4:6)]

# redundancy analysis
wrasse.rda <- rda(gen.imp ~ ., data=pred, scale=T)
wrasse.rda

# adjusted r2 
RsquareAdj(wrasse.rda)

# variance explained by each canonical axis
summary(eigenvals(wrasse.rda, model = "constrained"))

# most variation is in first CA
screeplot(wrasse.rda)

# Variance Inflation Factors for the predictor variables used in the model
vif.cca(wrasse.rda)

# shows are predictors are highly correlated
# set up pretty mapping factors and colors
levels(poptab_gen$Sample_ID) <- c("Flatanger","Austevoll","Stavanger","Kungsbacka","Strömstad","Kristiansand")
eco <- as.factor(poptab_gen$Sample_ID)
bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#33a02c","#e31a1c")

# axes 1 & 2
plot(wrasse.rda, type="n", scaling=3)
points(wrasse.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(wrasse.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) # the fish
text(wrasse.rda, scaling=3, display="bp", col="red", cex=1)                           # the predictors
legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)


# identify RDA candidates
load.rda <- summary(wrasse.rda)$species[,1:3]

hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 

# function to ID SNPs
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x) ## find loadings +/- z SD from mean loading     
  x[x < lims[1] | x > lims[2]]           # locus names in these tails
}

# find candidates
cand1 <- outliers(load.rda[,1], 3) 
cand2 <- outliers(load.rda[,2], 3) 
cand3 <- outliers(load.rda[,3], 3) 

# combine candidates into one object
wrasse.rda.cand <- c(names(cand1), names(cand2), names(cand3)) ## just the names of the candidates

# are there duplicates across axes?
length(wrasse.rda.cand[duplicated(wrasse.rda.cand)]) ## duplicate detection (detected on multiple RDA axes)

# write out unique candidates only
wrasse.rda.cand <- wrasse.rda.cand[!duplicated(wrasse.rda.cand)] ## unique candidates 

# Set up the color scheme for plotting:
bgcol  <- ifelse(colnames(gen.imp) %in% wrasse.rda.cand, 'gray32', '#00000000')
snpcol <- ifelse(colnames(gen.imp) %in% wrasse.rda.cand, 'red', '#00000000')

## a.es 1 & 2 - zooming in to just the SNPs here...
plot(wrasse.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), main="Wrasse RDA, axes 1 and 2")
points(wrasse.rda, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3)
points(wrasse.rda, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3)
text(wrasse.rda, scaling=3, display="bp", col="#0868ac", cex=1)

## a.es 2 & 3
plot(wrasse.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3), main="Wrasse RDA, axes 2 and 3")
points(wrasse.rda, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3, choices=c(2,3))
points(wrasse.rda, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3, choices=c(2,3))
text(wrasse.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(2,3))

# correlation of axes with predictors
intersetcor(wrasse.rda)[,1:3]

#############################
# Latent Factor Mixed Model #
#############################

# PCA on environmental predictors
pred.pca <- rda(pred,scale=T)
summary(pred.pca)$cont

# screeplot
screeplot(pred.pca, main = "Screeplot of Env Predictor Variables with Broken Stick", bstick=TRUE, type="barplot")

## correlations between the PC axis and predictors:
round(scores(pred.pca, choices=1:4, display="species", scaling=0), digits=3)

# store PC1 scores in new variable
pred.PC1 <- scores(pred.pca, choices=1, display="sites", scaling=0)

# PCA on genotypes
gen.pca <- rda(gen.imp, scale=T)
screeplot(gen.pca, main = "Screeplot of Genetic Data with Broken Stick", bstick=TRUE, type="barplot")

# set K value
K <- 3

# run LFMM
wrasse.lfmm <- lfmm_ridge(Y=gen.imp, X=pred.PC1, K=K) ## change K as you see fit

# calculate test statistics for the predictors
wrasse.pv <- lfmm_test(Y=gen.imp, X=pred.PC1, lfmm=wrasse.lfmm, calibrate="gif")

names(wrasse.pv) # this object includes raw z-scores and p-values, as well as GIF-calibrated scores and p-values

# genomic inflation factor
wrasse.pv$gif

# histograms of p-values
hist(wrasse.pv$pvalue[,1], main="Unadjusted p-values")        
hist(wrasse.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

# adjust p-values to q-values
wrasse.qv <- qvalue(wrasse.pv$calibrated.pvalue)$qvalues
length(which(wrasse.qv < 0.1)) ## h.w many SNPs have an FDR < 10%?

# which loci are outliers?
(wrasse.FDR.1 <- colnames(gen.imp)[which(wrasse.qv < 0.1)]) ## identify which SNPs these are

# fix loci names to match PCAdapt
wrasse.FDR.1<-gsub("\\.0", "",wrasse.FDR.1)
wrasse.FDR.1<-gsub("\\.1", "",wrasse.FDR.1)
wrasse.rda.cand<-gsub("\\.0", "",wrasse.rda.cand)
wrasse.rda.cand<-gsub("\\.1", "",wrasse.rda.cand)

# overlap between RDA and LFMM
intersect(wrasse.FDR.1, wrasse.rda.cand) ## found by both LFMM and RDA

# LFMM and PCAdapt
intersect(wrasse.FDR.1, outliersName)

# LFMM and PCAdapt
intersect(wrasse.rda.cand, outliersName)
