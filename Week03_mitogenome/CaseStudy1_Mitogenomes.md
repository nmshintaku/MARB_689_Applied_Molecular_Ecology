## MARB 689 Case Study 1: What are mitochondrial markers/genomes good for?
January 27, 2025

#### NAME:

### Identify a unknown specimen based on mtDNA sequences
Last week your colleague sent back assembled mitochondrial genomes from a collaborative whole genome sequencing project you started last year. The problem is your colleague mixed up the sequencing barcodes and now its unclear what specimen each mitogenome belongs to. They hope you can help sort this out.

1. So, where to begin? What mitochondrial markers could be helpful to determine whose who? Why?
<br/>
<br/>

2. Pick one of the three mitochondrial genomes in `02_mitogenome` directory. Circle:   A   B   C  

To test the marker(s) you identified above, let's annotate your chosen mitogenome. In your `class_working_directory`, copy over the `annotate.sh` batch script from the `02_mitogenome` in the group space. Open that script using the text editor `nano` and change the defined variables at the top of the script.

```
# Change the information below
USER=kitchens
SPECIMEN=A
```

Once you've updated the variables, execute the script. (**hint**- sbatch \<scriptName>)

3. a) What type of files did this job produce?
   b) What directory are they stored in?
   (**hint**- look up the tool used to see what files it produces)
<br/>
<br/>
<br/>
<br/>

Pull out one of the mitochondrial marker genes you identified above using `grep`. (**hint**- navigate to `/scratch/group/kitchen-group/MARB_689_Molecular_Ecology/class_working_directories/[netID]/[A/B/C]_unknown/[A/B/C]_unknown_MitoFinder_mitfi_Final_Results` to find the gene predictions).See prior notes for examples on how to use `grep` to match strings in the output. Create a new `.fa` file to paste the sequence(s) in.

Let's see if this marker gene shares sequence homology with sequences deposited in the NCBI nt database. Copy the `blast_nt.sh` script from the `02_mitogenome` directory into your working directory and modify the header with the requested information. Then, execute the script. This script takes a very LONG time to run, but submit it anyways to see how long it actually takes.

BLAST ( **B**asic **L**ocal **A**lignment **S**earch **T**ool) is a tool to identify similar regions between sequences by comparing the query (your sequence(s)) against a database of known sequences (nt= nucleotide database, nr= protein database). It allows for variation in sequences including gaps, insertions and deletions. There are different types of BLAST depending on the data type of the query and the database searched. In this case, you used `blastn` to search a nucleotide sequence against the nucleotide (nt) database.

Other BLAST commands:
_blastx_ - nucleotide query against protein database
_blastp_ - protein query against protein database
_tblastx_ - nucleotide query against translated nucleotide database
_tblastn_ - protein query against translated nucleotide database

Copy your sequence and submit it to the online version of BLAST can be found here: https://blast.ncbi.nlm.nih.gov/Blast.cgi

4. What type of output did this job produce?
<br/>
<br/>

5. a) Based on the BLAST report online, what is the closest organism that your sequence matched in the nt database?
<br/>

b) What is the Max Score of the top hit?
<br/>

c) What is Query Coverage of the top hit?
<br/>

d) What is the Percent Identity of the top hit?
<br/>

e) What is the E-value? Is this considered a good E-value?
<br/>


6. Try another mitochondrial gene. Do the results converge on the same top hit? If not, why do you think this could happen?
<br/>
<br/>

7. Which marker are you more confident for your taxonomic assignment? Why?
<br/>
<br/>

8. Does it appear that the mitochondrial genome of this specimen is present in the nt database?
<br/>
<br/>

### Pros and Cons of mitochondrial markers and mitogenomic analyses
9. Now that you've worked with mitochondrial markers, list 3 pros and 3 cons of using mitochondrial markers in the space below. For each, put a literature reference to support your point.

Pro 1:
<br/>
<br/>
<br/>
<br/>

Pro 2:
<br/>
<br/>
<br/>
<br/>

Pro 3:
<br/>
<br/>
<br/>
<br/>

Con 1:
<br/>
<br/>
<br/>
<br/>

Con 2:
<br/>
<br/>
<br/>
<br/>

Con 3:
<br/>
<br/>
<br/>
<br/>

### How many reads does it take to assemble a mitogenome?
The raw reads used to assemble the three specimens were also sent by your colleague. You wonder how many of the reads are needed to get a high quality mitogenome assembly and whether you can sequence at a lower depth (=save money!) in the future.

To test this, copy the `downsample.sh` script into your working directory from the `02_mitogenome` group directory. Like before, change the variables at the top of the script (under the header). In this you will need to define the number of reads you'd like to randomly sample. Shoot for somewhere between 1 to 10 million. Run the `downsample.sh` script.

Now, check the parameters for the assembler NovoPlasty on GitHub: https://github.com/ndierckx/NOVOPlasty. Copy the `config_specimen<?>.txt` from the `02_mitogenome` group directory into your directory. Also, copy the `Assemble.sh` into your working directory.

Once you've set up the script and configuration file with your variables/parameters, execute the `Assemble.sh` script.

10. How does the sampled assembled mitogenome you produced compare to the assembly from all reads?
<br/>
<br/>
<br/>
<br/>
11. How does it compare to the reference mitogenome? What are the number of variants (mutations)?
<br/>
<br/>
<br/>
<br/>

12. What species is the reference mitogenome? What is the lowest taxonomic group (Phylum, Class, Order, Family, Genus, Species) that your specimen and the reference specimen share?
<br/>
<br/>
<br/>
<br/>
