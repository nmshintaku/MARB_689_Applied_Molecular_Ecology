.libPaths("/scratch/group/kitchen-group/Rlibs")
library("vcfR")
library(adegenet)
library(poppr)
library(reshape2)
library("ggplot2")
library(tidyverse)
library(LEA)

#################
## Import Data ##
#################

#Read in genome VCF file 
vcf_gen <- read.vcfR("/scratch/group/kitchen-group/07_hybrids/west.filt.maf0.01.recode.vcf")

# subset to 200 SNPs used later
loc<-read.table("/scratch/group/kitchen-group/07_hybrids/loci.txt",check.names=FALSE, header=T, na.strings = c("", "NA"),
                stringsAsFactors = FALSE, sep="\t")

gt<-vcf_gen[vcf_gen@fix[,1] %in% loc$LocusNames]

# convert vcf to genind file format
gind <- vcfR2genind(vcf_gen)
gind

#add population or species information to the genind pop slot
poptab_gen<-read.table("/scratch/group/kitchen-group/07_hybrids/Sampleinfo_metadata.txt",
                       check.names=FALSE, header=T, na.strings = c("", "NA"),
                       stringsAsFactors = FALSE, sep="\t")

head(poptab_gen)

# add region to the genlight strata
gind@pop <- as.factor(poptab_gen$Sample_ID)
strata(gind) <- data.frame(poptab_gen$`age(years)`,poptab_gen$Sex,poptab_gen$`total_length(mm)`,poptab_gen$Sample_ID)

#####################
## Data filtering ###
#####################
# filter the SNPs based on missing data 
loci <- missingno(gind, type = "loci", cutoff = 0.05, quiet = FALSE, freq = FALSE)
loci

loci.individuals <- missingno(loci, type = "geno", cutoff = 0.05, quiet = FALSE, freq = FALSE)
loci.individuals

# filter the SNPs based on the MAF
loci.individuals.maf <- informloci(loci.individuals, cutoff = 2/nInd(loci.individuals), MAF = 0.05, quiet = FALSE)
loci.individuals.maf


##################################
## Principal Component Analysis ##
##################################
# perform PCA on filtered SNPs
x.pca <- dudi.pca(tab(loci.individuals.maf, NA.method = "mean"), scannf = FALSE, scale = FALSE)

# how much variation(inertia) is covered by PC1 or PC2?
summary(x.pca)

# plot the PCA results
s.class(x.pca$li, pop(loci.individuals.maf), col = funky(6),cstar=TRUE, label="")


##########
## DAPC ##
##########
# perform DAPC analysis on filtered SNPs
d.dapc <- dapc(loci.individuals.maf, n.pca = 20, n.da = 2)

# see how the samples were assigned to the regions
summary(d.dapc)

# let's plot the results
scatter(d.dapc,1,2, col= funky(6), bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4,cstar=TRUE, label="")

# look at just LD1
scatter(d.dapc,1,1, col= funky(6), bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)

############
### sNMF ###
############
# load in all SNPs
file.copy("/scratch/group/kitchen-group/07_hybrids/west.filt.maf0.01.recode.vcf",".")
wgeno<-vcf2geno("west.filt.maf0.01.recode.vcf")

set.seed(10032023)

# run sNMF, a STRUCTURE-like analysis
w_snmf <- snmf(input.file = wgeno,
               K = 1:6,
               entropy = TRUE,
               repetitions = 5,
               project = "new",
               alpha = 100
)

# what is the best K?
plot(w_snmf, cex = 1.2, col = "lightblue", pch = 19)

# which run (K) has the best fit?
ce <-  cross.entropy(w_snmf, K = 3)
ce
best_run <- which.min(ce)

# extract probability assignments of individuals into a table
q_mat <- LEA::Q(w_snmf, K = 3, run = best_run) 
colnames(q_mat) <- paste0("P", 1:3)
head(q_mat)

# convert in dataframe with information from the population table
q_df <- q_mat %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(individual = poptab_gen$Individual_ID,
         region = poptab_gen$Sample_ID)
q_df

# more data manipulation for plotting
q_df_long <- q_df %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q") 
q_df_long

# setting colors
myCol_2<-c("#B33B77","#FBB116","#6F81BC","white", "black",'#B0D09F',"#6FAA50")

# let's plot it!
q_df_long %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_manual(values = myCol_2, labels = c("P1", "P2", "P3")) +
  labs(fill = "Region") +
  theme_minimal() +
  # some formatting details to make it pretty
  theme(panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text( size=6, angle=90,hjust=2),
        strip.background = element_rect(fill = "transparent", color = "black"),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()
  )

## Now, repeat the steps above but with a subset of 200 SNPs
# convert subset of SNPs into geno format required by sNMF
w_tidy <- gt %>% 
  extract_gt_tidy() %>% 
  select(-gt_DP, -gt_PS, -gt_AD, -gt_GT_alleles, -gt_GQ) %>% 
  mutate(gt_GT= str_replace(gt_GT, "\\|+", "/"),
         gt1 = str_split_fixed(gt_GT, "/", n = 2)[,1],
         gt2 = str_split_fixed(gt_GT, "/", n = 2)[,2],
         geno_code = case_when(
           # homozygous for reference allele = 0
           gt1 == 0 & gt2 == 0 ~ 0,
           # heterozygous = 1
           gt1 == 0 & gt2 == 1 ~ 1,
           gt1 == 1 & gt2 == 0 ~ 1,
           # homozygous for alternate allele = 2
           gt1 == 1 & gt2 == 1 ~ 2,
           # missing data = 9
           gt1 == "" | gt2 == ""  ~ 9
         )) %>% 
  select(-gt_GT, -gt1, -gt2)

# more reformatting
wgeno_200 <- w_tidy %>% 
  pivot_wider(names_from = Indiv, values_from = geno_code) %>% 
  select(-Key)

# write geno file to your working directory
write.table(wgeno_200, 
            "geno.geno",
            col.names = FALSE,
            row.names = FALSE,
            sep = "")

# load in geno file you just made (sNMF is particular about this step)
w_snmf_200 <- snmf(input.file = "geno.geno",
                   K = 1:6,
                   entropy = TRUE,
                   repetitions = 5,
                   project = "new",
                   alpha = 100
)

plot(w_snmf_200, cex = 1.2, col = "lightblue", pch = 19)

# which run has the best fit?
ce <-  cross.entropy(w_snmf_200, K = 3)
ce
best_run <- which.min(ce)

q_mat <- LEA::Q(w_snmf_200, K = 3, run = best_run) 

colnames(q_mat) <- paste0("P", 1:3)

head(q_mat)

q_df <- q_mat %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(individual = poptab_gen$Individual_ID,
         region = poptab_gen$Sample_ID)
q_df

q_df_long <- q_df %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q") 

q_df_long

q_df_long %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_manual(values = myCol_2, labels = c("P1", "P2", "P3")) +
  labs(fill = "Region") +
  theme_minimal() +
  # some formatting details to make it pretty
  theme(panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text( size=6, angle=90,hjust=2),
        strip.background = element_rect(fill = "transparent", color = "black"),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()
  )


###############
## SNAPCLUST ##
###############

### SNAPCLUST - hybrid detection with subset of SNPs
gt_gind<-vcfR2genind(gt)
# add region to the genlight strata
gt_gind@pop <- as.factor(poptab_gen$Sample_ID)
strata(gt_gind) <- data.frame(poptab_gen$`age(years)`,poptab_gen$Sex,poptab_gen$`total_length(mm)`,poptab_gen$Sample_ID)
gt_gind

## optimal K, look for k = 1:8
a.aic <- snapclust.choose.k(8, gt_gind)
#pdf ("optimalK_aic.pdf", width=10, height=7)
plot(a.aic, type = "b", cex = 2, xlab = "k", ylab = "AIC")
points(which.min(a.aic), min(a.aic), col = "blue", pch = 20, cex = 2)
#dev.off()

## same data, using BIC
a.bic <- snapclust.choose.k(8, gt_gind, IC = BIC)
#pdf ("optimalK_bic.pdf", width=10, height=7)
plot(a.bic, type = "b", cex = 2, xlab = "k", ylab = "BIC")
points(which.min(a.bic), min(a.bic), col = "red", pch = 20, cex = 2)
#dev.off()

# Cluster assignment based on optimal K above
res.default <- snapclust(loci.individuals.maf[order(loci.individuals.maf$strata$poptab_gen.Sample_ID)], k = 2,hybrids = FALSE, pop.ini = "kmeans")
compoplot(res.default,show.lab = TRUE,posi="topright",cex.names = 0.65)

# Two samples were assiged to the southern population. Now we will test for different hybrid classes in the data to see if they are possible hybrids.

#F1 hybrids
res.hyb <- snapclust(gt_gind[order(gt_gind$strata$poptab_gen.Sample_ID)], k = 2, hybrids = TRUE)
#pdf ("K2_F1_SNPs.pdf", width=10, height=5)

compoplot(res.hyb,  show.lab = TRUE, col =myCol_2,
          posi="topright",cex.names = 0.5)
#dev.off()

# log-likelihood estimate of fit
res.hyb$ll

#F1 and 1st generation backcrosses
res2.back <- snapclust(gt_gind[order(gt_gind$strata$poptab_gen.Sample_ID)], k=2, hybrids = TRUE, n.starts= 15, max.iter=10, 
                      hybrid.coef = c(0.25, 0.5),cex.names = 0.65)

#pdf ("K2_F1andBC_SNPs.pdf", width=10, height=5)
compoplot(res2.back,show.lab = TRUE,col = myCol_2,posi="topright",cex.names = 0.5)
#dev.off()

# log-likelihood estimate of fit
res2.back$ll

#F1, 1st and 2nd generation backcrosses
res3.back <- snapclust(gt_gind[order(gt_gind$strata$poptab_gen.Sample_ID)], k=2, hybrids = TRUE, hybrid.coef = c(0.25,0.5,0.125))
#pdf ("K2_laterG_SNPs.pdf", width=10, height=5)
compoplot(res3.back,show.lab = TRUE,col = myCol_2,posi="topright",cex.names = 0.5)
#dev.off()

# log-likelihood estimate of fit
res3.back$ll

###############################
## D-statistic and F4 Ration ##
###############################

# This analysis was run on the HPRC command line using the code below in the 07_hybrids directory
#./Dsuite/Build/Dsuite Dtrios -c -n wrasse west.filt.maf0.01.recode.vcf SETs.txt

# read in the Dstat table
Dstat<-read.table("/scratch/group/kitchen-group/07_hybrids/SETs_wrasse_BBAA.txt", header=T)
Dstat

# let's assign pop names
pops <- c("Flatanger", "Stavanger", "Kristiansand", "Stromstad","Kungsbacka")

# for the analysis I randomly choose Kungsbacka as the outgroup

# plot heatmap of D-stat and F4 ratios
ggplot(Dstat %>% filter(p.value < 0.05),aes(P2, P3, alpha =p.value)) + geom_tile(aes(fill = Dstatistic))+
  scale_fill_distiller(type = "seq",direction = 1,palette = "Blues") + scale_alpha(range = c(1, 0.2))

ggplot(Dstat %>% filter(p.value < 0.05),aes(P2, P3, alpha =p.value)) + geom_tile(aes(fill = f4.ratio))+
  scale_fill_distiller(type = "seq",direction = 1,palette = "Blues")+ scale_alpha(range = c(1, 0.2))
