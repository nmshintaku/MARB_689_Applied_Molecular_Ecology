### Variant Effect Predictor
Another popular variant annotation tool is the Ensembl Variant Effect Predictor (VEP). This package is used commonly in model systems but has only limited support for new genome assemblies that not present in the Ensembl database.

Here is an example of how to run the program using the subset of coral SNPs.
load module on FASTER:
```
module load GCC/11.3.0
module load VEP/107
```

Compress and index the reference GFF file. This is required for VEP to run.
```
grep -v "#" Kandelia_obovata.gff3 | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > Kandelia_obovata.gff.gz
tabix -p gff Kandelia_obovata.gff.gz
```

Run VEP:
```
vep -i hqfilter_snp_Kandelia_obovata.vcf --format vcf \
--gff /scratch/group/kitchen-group/MARB_689_Molecular_Ecology/map_reference/Kandelia_obovata/Kandelia_obovata.gff.gz \
--fasta /scratch/group/kitchen-group/MARB_689_Molecular_Ecology/map_reference/Kandelia_obovata/Kandelia_obovata_genome.fa \
--tab -o Ko_VEP.txt
```

If you want to find where your SNPs intersect genes, you can use bedtools:
```
# load modules
module load GCC/13.2.0
module load BEDTools/2.31.1

# cut chromosome and position
cut -f 1,2  hqfilter_rd2_snp_Kandelia_obovata.recode.vcf > Ko_SNP.bed

# make bed file format
awk  -v OFS='\t'  '{print $1, $2, $2}' hqfilter_rd2_snp_Kandelia_obovata.recode.vcf > Ko_SNP.bed
sed '/^#/d' Ko_SNP.bed > Ko_SNP2.bed

# intersect the position with the gff file
bedtools intersect -wb -a Ko_SNP2.bed -b /scratch/group/kitchen-group/MARB_689_Molecular_Ecology/map_reference/Kandelia_obovata/Kandelia_obovata.gff3 | grep "mRNA" | cut -f 1,2,7,8,12 | sed "s/;/\t/g" | cut -f 1-5 | sed "s/ID=//" > Ko_SNP_annotated.tsv

```

To then know what genes these are, you will need to join this table with the annotation table in `map_reference` directories.