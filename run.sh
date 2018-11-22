## Author: Ariadna Cilleros Portet
## Date: 13/10/2018
## eQTL Hands-On

#Task 2: Subset the input VCF file, keeping only biallelic SNPs and indels with MAF ≥ 0.05 and those sample IDs for which we have both expression and genotype data (for duplicated IDs keep only the first). To do that, we will use bcftools. Besides, discard those variants with less than 10 individuals per genotype group.
eqtl@ubuntu:~$ cd eQTL
eqtl@ubuntu:~/eQTL$ sudo setxkbmap -layout es
eqtl@ubuntu:~/eQTL$ sudo docker run -v $PWD:$PWD -w $PWD -it dgarrimar/eqtlmapping
# Get GEUVADIS samples from the metadata
root@0f637262ab49:/home/eqtl/eQTL# cut -f1 input/unprocessed/geuvadis/geuvadis.metadata.txt | sed '1d' | sort | uniq > tmp/geuvadis.samples.txt 

# Subset the VCF (common samples, biallelic SNPs and indels, MAF >= 0.05, no duplicates)
root@0f637262ab49:/home/eqtl/eQTL# bcftools view -v snps,indels -m 2 -M 2 -q 0.05:minor -S tmp/geuvadis.samples.txt -Ob input/unprocessed/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | bcftools norm -d all -Oz -o tmp/genotypes.chr22.vcf.gz
Lines   total/split/realigned/skipped:	105298/0/0/0

# Subset the VCF so that there are at least 10 individuals per genotype group and compress it (for indexing we require 'bgzip' compression)
root@46f46e2b44aa:/home/eqtl/eQTL# filter.genotype.py -t 10 -g <(zcat tmp/genotypes.chr22.vcf.gz) | bgzip > input/processed/genotypes.chr22.vcf.gz

#Index the VCF
root@46f46e2b44aa:/home/eqtl/eQTL# tabix -p vcf input/processed/genotypes.chr22.vcf.gz

#Q1: What do the bcftools options employed mean? 
The following options are common to many bcftools commands. See usage for specific commands to see if they apply.

FILE
    Files can be both VCF or BCF, uncompressed or BGZF-compressed. The file "-" is interpreted as standard input. Some tools may require tabix- or CSI-indexed files. 
-c, --collapse snps|indels|both|all|some|none|id

    Controls how to treat records with duplicate positions and defines compatible records across multiple input files. Here by "compatible" we mean records which should be considered as identical by the tools. For example, when performing line intersections, the desire may be to consider as identical all sites with matching positions (bcftools isec -c all), or only sites with matching variant type (bcftools isec -c snps  -c indels), or only sites with all alleles identical (bcftools isec -c none).

    none
        only records with identical REF and ALT alleles are compatible 
    some
        only records where some subset of ALT alleles match are compatible 
    all
        all records are compatible, regardless of whether the ALT alleles match or not. In the case of records with the same position, only the first will be considered and appear on output. 
    snps
        any SNP records are compatible, regardless of whether the ALT alleles match or not. For duplicate positions, only the first SNP record will be considered and appear on output. 
    indels
        all indel records are compatible, regardless of whether the REF and ALT alleles match or not. For duplicate positions, only the first indel record will be considered and appear on output. 
    both
        abbreviation of "-c indels  -c snps" 
    id
        only records with identical ID column are compatible. Supported by bcftools merge only. 

-f, --apply-filters LIST
    Skip sites where FILTER column does not contain any of the strings listed in LIST. For example, to include only sites which have no filters set, use -f .,PASS. 
--no-version
    Do not append version and command line information to the output VCF header. 
-o, --output FILE
    When output consists of a single stream, write it to FILE rather than to standard output, where it is written by default. 
-O, --output-type b|u|z|v
    Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when piping between bcftools subcommands to speed up performance by removing unnecessary compression/decompression and VCF←→BCF conversion. 
-r, --regions chr|chr:pos|chr:from-to|chr:from-[,…]
    Comma-separated list of regions, see also -R, --regions-file. Note that -r cannot be used in combination with -R. 
-R, --regions-file FILE
    Regions can be specified either on command line or in a VCF, BED, or tab-delimited file (the default). The columns of the tab-delimited file are: CHROM, POS, and, optionally, POS_TO, where positions are 1-based and inclusive. The columns of the tab-delimited BED file are also CHROM, POS and POS_TO (trailing columns are ignored), but coordinates are 0-based, half-open. To indicate that a file be treated as BED rather than the 1-based tab-delimited file, the file must have the ".bed" or ".bed.gz" suffix (case-insensitive). Uncompressed files are stored in memory, while bgzip-compressed and tabix-indexed region files are streamed. Note that sequence names must match exactly, "chr20" is not the same as "20". Also note that chromosome ordering in FILE will be respected, the VCF will be processed in the order in which chromosomes first appear in FILE. However, within chromosomes, the VCF will always be processed in ascending genomic coordinate order no matter what order they appear in FILE. Note that overlapping regions in FILE can result in duplicated out of order positions in the output. This option requires indexed VCF/BCF files. Note that -R cannot be used in combination with -r. 
-s, --samples [^]LIST
    Comma-separated list of samples to include or exclude if prefixed with "^". The sample order is updated to reflect that given on the command line. Note that in general tags such as INFO/AC, INFO/AN, etc are not updated to correspond to the subset samples. bcftools view is the exception where some tags will be updated (unless the -I, --no-update option is used; see bcftools view documentation). To use updated tags for the subset in another command one can pipe from view into that command. 
 -S, --samples-file FILE
    File of sample names to include or exclude if prefixed with "^". One sample per line. See also the note above for the -s, --samples option. The sample order is updated to reflect that given in the input file. The command bcftools call accepts an optional second column indicating ploidy (0, 1 or 2) or sex (as defined by --ploidy, for example "F" or "M")
 -t, --targets [^]chr|chr:pos|chr:from-to|chr:from-[,…]
    Similar as -r, --regions, but the next position is accessed by streaming the whole VCF/BCF rather than using the tbi/csi index. Both -r and -t options can be applied simultaneously: -r uses the index to jump to a region and -t discards positions which are not in the targets. Unlike -r, targets can be prefixed with "^" to request logical complement. For example, "^X,Y,MT" indicates that sequences X, Y and MT should be skipped. Yet another difference between the two is that -r checks both start and end positions of indels, whereas -t checks start positions only. Note that -t cannot be used in combination with -T. 
-T, --targets-file [^]FILE
    Same -t, --targets, but reads regions from a file. Note that -T cannot be used in combination with -t. 
    With the call -C alleles command, third column of the targets file must be comma-separated list of alleles, starting with the reference allele. Note that the file must be compressed and index.
 --threads INT
    Number of output compression threads to use in addition to main thread. Only used when --output-type is b or z. Default: 0.

#Q2: How many variants do you get in input/processed/genotypes.chr22.vcf.gz? 
root@46f46e2b44aa:/home/eqtl/eQTL# zcat input/processed/genotypes.chr22.vcf.gz | grep -v "#" | wc -l 
#74656

#Q3: How many samples do you have before and after subsetting? Hint: The genotype file contains information from many samples, but only a subset of them has gene expression data in the GEUVADIS project. Note that the GEUVADIS metadata contains duplicated sample IDs, as gene expression in GEUVADIS is calculated from paired-end RNA-Seq data. You can use bcftools to print the sample IDs and count them.
root@46f46e2b44aa:/home/eqtl/eQTL# bcftools stats input/unprocessed/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > tmp/stats.before
root@46f46e2b44aa:/home/eqtl/eQTL# less tmp/stats.before
#2504
bcftools stats input/processed/genotypes.chr22.vcf.gz > tmp/stats.after
root@46f46e2b44aa:/home/eqtl/eQTL# less tmp/stats.before
#445

#Task 3: Use the GENCODE human gene annotation to obtain the information above. Q1: Which version of GENCODE is GEUVADIS using? v12
#Q2: To which genome assembly does this annotation correspond? GRCh37
#Q3: How many protein coding genes are annotated in the last version? 19940
#Hint: Check GEUVADIS paper. To run the scripts in bin without the need of bash script.sh or ./bin/script.sh add them to the PATH environment variable. 
#Q4: Which command do you use to do this?
PATH=$PATH:~/eQTL/bin

#----
# Set the variable 'release' to the version of GENCODE used in GEUVADIS (e.g. `release=99`) and download the corresponding GENCODE annotation
root@46f46e2b44aa:/home/eqtl/eQTL# release=12
root@46f46e2b44aa:/home/eqtl/eQTL# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$release/gencode.v$release.annotation.gtf.gz
--2018-10-15 14:41:11--  ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_12/gencode.v12.annotation.gtf.gz
           => 'gencode.v12.annotation.gtf.gz'
Resolving ftp.ebi.ac.uk (ftp.ebi.ac.uk)... 193.62.192.4
Connecting to ftp.ebi.ac.uk (ftp.ebi.ac.uk)|193.62.192.4|:21... connected.
root@362679a2dbfe:/home/eqtl/eQTL# zcat input/unprocessed/gencode/gencode.annotation.gtf.gz | grep "gene_type \"protein_coding\"\|gene_type \"lincRNA\"" | gtf2bed.sh > tmp/gencode.annotation.bed

Logging in as anonymous ... Logged in!
==> SYST ... done.    ==> PWD ... done.
==> TYPE I ... done.  ==> CWD (1) /pub/databases/gencode/Gencode_human/release_12 ... done.
==> SIZE gencode.v12.annotation.gtf.gz ... 26249874
==> PASV ... done.    ==> RETR gencode.v12.annotation.gtf.gz ... done.
Length: 26249874 (25M) (unauthoritative)

gencode.v12.annotation.gtf.gz       100%[=================================================================>]  25.03M  7.43MB/s    in 3.7s    

2018-10-15 14:41:16 (6.76 MB/s) - 'gencode.v12.annotation.gtf.gz' saved [26249874]
root@362679a2dbfe:/home/eqtl/eQTL# PATH=$PATH:$PWD/bin 
root@362679a2dbfe:/home/eqtl/eQTL# ls  
bin  input  result  tmp
root@362679a2dbfe:/home/eqtl/eQTL# echo $PWD
/home/eqtl/eQTL
root@362679a2dbfe:/home/eqtl/eQTL# zcat input/unprocessed/gencode/gencode.annotation.gtf.gz | grep "gene_type \"protein_coding\"\|gene_type \"lincRNA\"" | gtf2bed.sh > tmp/gencode.annotation.bed
root@362679a2dbfe:/home/eqtl/eQTL# 

#Note that we implicitly have the information that we wanted. Q5: But how to get the TSS positions and the gene lengths from it? With 'awk' command.
#Q6: to which BED coordinates would correspond the GTF coordinates chr1 10 20? Why? chr1 9 20, because BED is 0 based and GTF is 1 based.
#Q7: Why do we need to use tmpfile below? To edit a file, in which you wanna work at the same time. 
#Hint Let's consider the TSS as the first base of the (first transcript of the) gene (5' -> 3'). Note BED positions are 0-indexed and GTF positions are 1-indexed (some extra help here).

# Compute gene lengths 
root@362679a2dbfe:/home/eqtl/eQTL# awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$3-$2,$6}' tmp/gencode.annotation.bed > tmpfile; mv tmpfile tmp/gencode.annotation.bed
# Compute TSS positions. Note that for genes in the '+' strand, the TSS is the start position, and for genes in the '-' strand it is the end position!
root@362679a2dbfe:/home/eqtl/eQTL# awk 'BEGIN{OFS="\t"}{if($6=="+"){print $1,$2,$2+1,$4,$5,$6}else{print $1,$3-1,$3,$4,$5,$6}}' tmp/gencode.annotation.bed > tmpfile; mv tmpfile tmp/gencode.annotation.bed
# Remove 'chr' from chromosome names (-i option to modify the file 'in place')
root@362679a2dbfe:/home/eqtl/eQTL# sed -i "s/^chr//" tmp/gencode.annotation.bed

#Task 4: Obtain the gene expression file in the format required by QTLtools. Hint Remember that the gene expression file contains data for all genes in the annotation, while the BED file that you just created only contains PC and lincRNA.

# Join the bed file with the expression file
# (Both files should be row-ordered by gene ID. Column order and header are lost in the output file)
root@ab402860b250:/home/eqtl/eQTL# join -1 4 -2 1 -t $'\t' <(sort -k4,4 tmp/gencode.annotation.bed) <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | sort -k1,1) > tmp/joint.tsv

# Subset chr22 (same as the VCF file)
root@ab402860b250:/home/eqtl/eQTL# awk '$2==22' tmp/joint.tsv > tmp/joint.chr22.tsv

# Recover the column order, sort rows by chr and start position (WARNING: this command may not work within the docker container for WSL users)
root@ab402860b250:/home/eqtl/eQTL# paste <(awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$5,$6}' tmp/joint.chr22.tsv) <(cut -f1-6 --complement tmp/joint.chr22.tsv) | sort -k1,1V -k2,2n > tmp/joint.chr22.bed

# Recover the header
root@ab402860b250:/home/eqtl/eQTL# cat <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | head -1 | sed "s/TargetID/#chr\tstart\tend\tgene\tlength\tstrand/") tmp/joint.chr22.bed > tmp/genes.chr22.rpkm.bed

#The last step in gene expression pre-processing is gene expression normalization. We will also remove lowly expressed genes accross all the samples. Q1: Of all genes considered, which have lower expression levels, protein-coding or lincRNA? Protein?coding genes hae less expression.
#Q2: Why do we need gene expression to be normal? Porqué sinó cuando hagamos la normalización por quartiles, los valores se veran arrastrados por unos erroneos. 
#Q3: How would you check that quantile normalization worked? Because if you plot he distribution of the genes and the samples, you will see that they will be really similar. 
#Q4: and that gene expression of a gene follows a normal distribution? Doing a plot you will see a bell curve (Gauss Distribution).
 
# Run gene expression normalization (quantile normalization + gene expression to normal distribution)
# Filter out genes with less than 0.1 RPKM in 50% of the samples
root@93979a29ccf4:/home/eqtl/eQTL# normalize.R -i tmp/genes.chr22.rpkm.bed -o tmp/genes.chr22.norm.bed

# Compress and index the final gene expression file
root@93979a29ccf4:/home/eqtl/eQTL# bgzip tmp/genes.chr22.norm.bed
root@93979a29ccf4:/home/eqtl/eQTL# tabix -p bed tmp/genes.chr22.norm.bed.gz
root@93979a29ccf4:/home/eqtl/eQTL# mv tmp/genes.chr22.norm.bed.gz* input/processed

#Task 5: Check that normalization worked! Hint: Plot the distribution of the values in rows (genes) and columns (samples). Compare the plots before and after normalization. Q1: What can you see? We can clearly see that after the normalitzation there's a normal distribution of the data, and before in the boxplot it's messy. 
 
#Before: 
root@a8f11dad9b87:/home/eqtl/eQTL# check.norm.R -i tmp/genes.chr22.rpkm.bed -o result/plots/check.norm.2.pdf

#After: 
root@a8f11dad9b87:/home/eqtl/eQTL# check.norm.R -i input/processed/genes.chr22.norm.bed.gz -o result/plots/check.norm.pdf

#Task 6: Review the genotype and gene expression metadata to identify potential covariates for your analysis. Q1: Which ones would you select? Hint: Maybe you find useful this small oneliner:
root@a8f11dad9b87:/home/eqtl/eQTL# head -1 input/unprocessed/geuvadis/geuvadis.metadata.txt | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++) {print i, $i}}'

#Source Name, Comment [ENA_SAMPLE], Characteristics [Organism], Term Source REF, Term Accession Number, Characteristics [Strain], Characteristics [population], Protocol REF, Technology Type

#Task 7: Let's use the pca tool in QTLtools to obtain PCs of our expression and genotype data. Q1: What do the parameters employed mean? Q2: Which information do the output files contain?
#This is going to produce two files: one that contains the individual coordinates on the Principal Components (PC), the other one contains the percentatges of the variance explained by each PC. 

# Expression PCA
root@a45e24868601:/home/eqtl/eQTL# QTLtools pca --bed input/processed/genes.chr22.norm.bed.gz --scale --center --out result/expression 

# Genotypes PCA
root@a45e24868601:/home/eqtl/eQTL# QTLtools pca --vcf input/processed/genotypes.chr22.vcf.gz --scale --center --maf 0.05 --distance 50000 --out result/genotypes

#Now let's plot the first two PCs in each case:

# Note the input should coincide with the output of QTLtools pca on each case
root@a45e24868601:/home/eqtl/eQTL# pcaPlot.R -i result/expression -o result/plots/expression.pca.pdf
root@a45e24868601:/home/eqtl/eQTL# pcaPlot.R -i result/genotypes -o result/plots/genotypes.pca.pdf

#Q3: What we can observe in the plot? In expression plot we see an homogenous distribution, while in the genotype we can see two different groups. But for more information the better option is to color the samples. 

#You can color the samples by the different values of a field in the corresponding metadata. For example, let's plot the first two genotype PCs and color the points by the super_pop field. Try other combinations!

root@a45e24868601:/home/eqtl/eQTL# pcaPlot.R -i result/genotypes --metadata input/unprocessed/1000g/1000g.phase3_metadata.txt --color super_pop --out result/plots/genotypes.pca.super_pop.pdf

#Q4: With this information, which covariates seem more relevant to explain the variability in the data?
#The Characteristics [population]. 

#We can also try to determine which fraction of the total variance in the expression data is explained by each potential covariate:

# Generate a common metadata with info about the population, gender and laboratory.
root@a45e24868601:/home/eqtl/eQTL# join -j 1 -t $'\t' <(sort -k1,1 input/unprocessed/1000g/1000g.phase3_metadata.txt) <(cut -f1,20 input/unprocessed/geuvadis/geuvadis.metadata.txt | sort -k1,1 | uniq) > tmp/metadata.txt

# Set names for the new metadata                                                                          
root@a45e24868601:/home/eqtl/eQTL# sed -i '1s/^/sampleID\tpop\tsuper_pop\tgender\tlab\n/' tmp/metadata.txt

# Build a linear model and plot the contribution of each factor in the metadata to the total variance
root@a45e24868601:/home/eqtl/eQTL# var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m tmp/metadata.txt --formula "~ (1|gender) + (1|pop) + (1|lab)" -o result/plots/vp.pdf


#Q5: Which are the factors that explain more variance? First the lab, then the population and at the end the gender.

#Task 8: Use PEER to infer hidden covariates from the expression matrix.
# Compute 10 PEER factors
root@a45e24868601:/home/eqtl/eQTL# peer.R -i input/processed/genes.chr22.norm.bed.gz -p 10 -o tmp/peer.tsv

# Check how much variance do the first 5 PEER explain in comparison with the known factors
root@a45e24868601:/home/eqtl/eQTL# var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m <(paste tmp/peer.tsv tmp/metadata.txt) -f "~ (1|pop) + (1|lab) + PEER1 + PEER2 + PEER3 + PEER4 + PEER5" -o result/plots/vp.peer.pdf

#Q1: How much variance do they explain? On average is it more or less than the explained by the known factors? In the case of PEER1 less than 10%, in PEER3 and PEER2 they have more or less the same percentatge but it-s also lower than 10%. And in PEER4 and PEER5 they maybe explain just 1% (more or less)It's mote than the explained by the known factors. 

#Finally generate a covariate file in the format required by QTLtools.

# 'Rscript -e' is just a trick to run an R script without opening an interactive R session in the console. ;)
root@a45e24868601:/home/eqtl/eQTL# join -j 1 -t $'\t' tmp/metadata.txt tmp/peer.tsv  | Rscript -e 'write.table(t(read.table(file("stdin", open = "r", blocking = T), h = F)), file = "input/processed/covariates.tsv", quote = F, sep = "\t", col.names = F, row.names = F)'

# Compress it
root@a45e24868601:/home/eqtl/eQTL# gzip input/processed/covariates.tsv

#Task 9: Run a nominal pass reporting all p-values in (0,1] and check the output file nominals.txt. 
root@7bc6d75f3030:/home/eqtl/eQTL# QTLtools cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --nominal 1 --out result/nominals.txt

#Have a look at QTLtools documentation for a detailed description of the options and the output format. Which information contains each field in the ouptut file? 1. The phenotype ID, 2. The chromosome ID of the phenotype, 3. The start position of the phenotype, 4. The end position of the phenotype, 5. The strand orientation of the phenotype, 6. The total number of variants tested in cis, 7. The distance between the phenotype and the tested variant (accounting for strand orientation),  8. The ID of the tested variant, 9. The chromosome ID of the variant, 10. The start position of the variant, 11. The end position of the variant, 12. The nominal P-value of association between the variant and the phenotype, 13. The corresponding regression slope, 14. A binary flag equal to 1 is the variant is the top variant in cis

#Q1: Are there pairs genotype-phenotype with exactly the same p-value and effect size (β)? How is this possible?
#p-value
root@7bc6d75f3030:/home/eqtl/eQTL# wc -l <(awk '{print $12}' result/nominals.txt)
1406669 /dev/fd/63
root@7bc6d75f3030:/home/eqtl/eQTL# wc -l <(awk '{print $12}' result/nominals.txt | uniq)
1250448 /dev/fd/63
#effect size
root@7bc6d75f3030:/home/eqtl/eQTL# wc -l <(awk '{print $13}' result/nominals.txt)
1406669 /dev/fd/63
root@7bc6d75f3030:/home/eqtl/eQTL# wc -l <(awk '{print $13}' result/nominals.txt | uniq)
1251254 /dev/fd/63

#Yes, maybe beacuse they are really related, like they are part from the same haplotype. 

#Have a look at the p-value distribution. Q2: What do you observe?
root@94200c766902:/home/eqtl/eQTL# pvdist.R -i result/nominals.txt --col 12 -o result/plots/pvdist.pdf

#The flat distribution along the bottom are the null p-values, which are uniformly distributed between 0 and 1, because that's part of a definition of a p-value: under the null, it has a 5% chance of being less than 0.05, a 10% chance of being less than 0.1, etc. This describes a uniform distribution. 

#Calculate LD between a pair of variants with exactly the same p-value and effect size using PLINK. Q3: Which SNPs did you select? What do you observe?
#For example we select: 
root@12ff02c99285:/home/eqtl/eQTL# grep 0.122366 result/nominals.txt 
ENSG00000100299.12 22 51066607 51066607 - 2548 94976 rs131797 22 50971631 50971631 0.00137699 0.122366 0
ENSG00000100299.12 22 51066607 51066607 - 2548 94968 rs131796 22 50971639 50971639 0.00137699 0.122366 0

#The LD tell us the haplotypes of this two snps, and its frequencies. --ld rs131797 rs131796:

#   R-sq = 1              D' = 1

   Haplotype     Frequency    Expectation under LE
   ---------     ---------    --------------------
          TG      0.265169                0.070314
     TAAAAAG     -0                       0.194854
         TGA     -0                       0.194854
    TAAAAAGA      0.734831                0.539977

   In phase alleles are TG/TAAAAAGA

#Task 10: Run a permutation pass and check the output file permutation.txt. Have a look at QTLtools documentation for a detailed description of the output format. Hint: This can be done faster using a 4 threads:
root@12ff02c99285:/home/eqtl/eQTL# QTLtools cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --permute 1000 --out result/permutations.txt


for j in $(seq 1 16); do
  echo "cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --permute 1000 --chunk $j 16 --out result/permutations_$j.txt"
done | xargs -P4 -n14 QTLtools
cat result/permutations_*.txt > result/permutations.txt; rm result/permutations_*.txt 

#With the plot we observe that the unsignificant beta approximated p-values are well calibrated given the empirical ones. 

#Task 11: Use the script mtc.R to perform multiple testing correction through different methods: i) Bonferroni on all tests, ii) FDR on all tests, iii) two-level: permutation + FDR. Set α to 0.05 (default). 
#Bonferroni
root@db25c5edc1f2:/home/eqtl/eQTL# mtc.R -n result/nominals.txt -p result/permutations.txt --method 'bonferroni' --alpha 0.05 --out tmp/bonferroni.txt
#FDR
root@db25c5edc1f2:/home/eqtl/eQTL# mtc.R -n result/nominals.txt -p result/permutations.txt --method 'fdr' --alpha 0.05 --out tmp/fdr.txt
#permutation + FDR
root@db25c5edc1f2:/home/eqtl/eQTL# mtc.R -n result/nominals.txt -p result/permutations.txt --method 'perm-fdr' --alpha 0.05 --out result/eqtls.tsv
Global empirical P-value threshold = 1.20e-02

#Q1: How many significant eQTLs do we find in each case in comparison with the nominal pass? Hint: Run Rscript mtc.R --help. Note the output file now has a header, and \t as column separator.
#Nominals (tiene un header)
root@7bc6d75f3030:/home/eqtl/eQTL# wc -l result/nominals.txt
1406669 result/nominals.txt
#Permutation + FDR
root@7bc6d75f3030:/home/eqtl/eQTL# wc -l result/eqtls.tsv
12045 result/eqtls.tsv

#Task 12: Now have a look at your results! Use eQTLviewer.R to make a plot for the top 10 eQTLs
root@7bc6d75f3030:/home/eqtl/eQTL# eQTLviewer.R -i <(head -n 10 result/eqtls.tsv) -g input/processed/genotypes.chr22.vcf.gz -e input/processed/genes.chr22.norm.bed.gz -o result/plots/eQTLs_head.pdf --verbose
Plot 1/9: gene ENSG00000100321.10 and snp rs909685
Plot 2/9: gene ENSG00000100321.10 and snp rs2069235
Plot 3/9: gene ENSG00000260065.1 and snp rs4275
Plot 4/9: gene ENSG00000260065.1 and snp rs147404731
Plot 5/9: gene ENSG00000100376.7 and snp rs738177
Plot 6/9: gene ENSG00000184674.7 and snp rs114323341
Plot 7/9: gene ENSG00000184674.7 and snp rs192096839
Plot 8/9: gene ENSG00000184674.7 and snp rs141325606
Plot 9/9: gene ENSG00000172404.4 and snp rs103197

#Task 13: Use the Ensembl Regulatory Build to assess the enrichment of eQTLs in annotated functional features.

# Download from ftp server
root@db25c5edc1f2:/home/eqtl/eQTL# rsync -av rsync://ftp.ensembl.org/ensembl/pub/grch37/release-86/regulation/homo_sapiens/AnnotatedFeatures.gff.gz input/unprocessed/ensembl

# Get chr, start, end and feature name in BED format
root@db25c5edc1f2:/home/eqtl/eQTL# zcat input/unprocessed/ensembl/AnnotatedFeatures.gff.gz | awk 'BEGIN{FS=OFS="\t"}{print $1, $4-1, $5, $9}' | sed -r 's/Name=([^;]+);.*/\1/' | grep -v '^GL' | sort -V > tmp/ERB.bed

# Merge overlapping features of the same type 
# e.g. chr1 100 200 feat1            chr1 100 300 feat1
#      chr1 150 300 feat1     =>     chr1 100 250 feat2
#      chr1 100 250 feat2

root@db25c5edc1f2:/home/eqtl/eQTL# for feat in $(cut -f4 tmp/ERB.bed | sort | uniq); do 
>   bedtools merge -i <(grep -Fw $feat tmp/ERB.bed) -c 4 -o distinct
> done > input/processed/ERB.collapsed.bed

# Remove 'chr' from chromosome names (-i option to modify the file 'in place')
sed -i "s/^chr//" input/processed/ERB.collapsed.bed

#Perform the enrichment of top eQTLs:
root@db25c5edc1f2:/home/eqtl/eQTL# for feat in $(cut -f4 tmp/ERB.bed | sort | uniq); do 
>   bedtools merge -i <(grep -Fw $feat tmp/ERB.bed) -c 4 -o distinct
> done > input/processed/ERB.collapsed.bed
root@db25c5edc1f2:/home/eqtl/eQTL# sed -i "s/^chr//" input/processed/ERB.collapsed.bed
root@db25c5edc1f2:/home/eqtl/eQTL# for feat in $(cut -f4 input/processed/ERB.collapsed.bed | sort | uniq); do 
>   QTLtools fenrich --qtl <(sed '1d' result/eqtls.tsv | awk '{print $9, $10-1, $10, $8, $1, "."}') --tss tmp/gencode.annotation.bed  --bed <(grep -Fw $feat input/processed/ERB.collapsed.bed) --out tmp/enrich.txt > /dev/null; echo "$(cat tmp/enrich.txt) $feat" 
> done | grep -Fwv inf | grep -Fwv nan > result/enrichments.txt


root@db25c5edc1f2:/home/eqtl/eQTL# plot.enrich.R -i result/enrichments.txt -o result/plots/enrich.pdf

#Q1: Which are the top enriched features? Which kind of factors are they? El que tiene un mayor odd ratio es H3K36me3 y es un marcador de histonas. Seguido de este enconctramos otros marcadores de histonas: H3K79me2, H3K4me1, H4K20me1, H3K9ac, H2AZ, H3K4me2, H3K4me3, H3K27ac, H3K27me3, etc.
#Q2: What does an odds ratio lower than one mean? SE han observado menos de lo que se esperava, o sea hay un enrichment bajo. 

#Task 14: Now use the Variant Effect Predictor (VEP) to assess the impact of the eQTL variants. Q1: Which kind of consequences have they, according to the VEP? In which proportion? intron_variant: 47% / upstream_gene_variant: 15% / downstream_gene_variant:14% / non_coding_transcript_variant: 11% / NMD_transcript_variant: 6% / regulatory_region_variant: 2% / intergenic_variant: 2% / non_coding_transcript_exon_variant: 1% / 3_prime_UTR_variant: 1%

#Q2: How many eQTLs are high impact variants? Which consequences are related to those high impact variants? 
root@d6add0c3151e:/home/eqtl/eQTL# wc -l MOCNDEmRPOJs68Io.IMPACT_is_HIGH\(1\).txt 
24 MOCNDEmRPOJs68Io.IMPACT_is_HIGH(1).txt #23 variants i +1 linea de cabezera

root@d6add0c3151e:/home/eqtl/eQTL# uniq <(cut -f4 MOCNDEmRPOJs68Io.IMPACT_is_HIGH\(1\).txt)
Consequence
splice_acceptor_variant,non_coding_transcript_variant
stop_gained
splice_donor_variant,frameshift_variant
frameshift_variant
stop_gained
splice_acceptor_variant
splice_acceptor_variant,NMD_transcript_variant
splice_acceptor_variant
splice_acceptor_variant,non_coding_transcript_variant
stop_gained

#Q3: Out of all high impact variants, how many of them are falling in acceptor splice sites of protein coding genes? Hint: Use the command below to get the variants that you will upload to the VEP.
#6

#Task 15: And what about eGenes? Perform a GO enrichment to learn more about their function. 
#Q1: In which biological processes are your eGenes enriched? 
#Response to stimulus: response to molecule of bacterial origin and response to lipopolysaccharide.
#Which molecular functions and components correspond to those processes? 
#The function is Ras guanyl-nucleotide exchange factor activity and the component is endoplasmic reticulum. 

#In input/unprocessed/gwas you will find the file gwas.catalog.hg19.bed, which contains a parsed version of the GWAS Catalog. We will use the Regulatory Trait Concordance (RTC) method to co-localize our eQTLs with the GWAS variants in the Catalog.

#Task 16: Perform a co-localization analysis.

#Q1: How many pairs of variants have a RTC value above 0.9? 
root@aa0862d815df:/home/eqtl/eQTL# wc -l <(awk '$20 > "0.9"' result/rtc.txt)
39 /dev/fd/63
#39 pairs

#Q2: For each pair, we have a GWAS hit and an eQTL. 
#Find one example so that the gene to which the eQTL is associated is relevant for the trait/disease to which the GWAS variant is associated. 
#Explore the literature and the biological databases that you know to gather more information. 
#Como más cerca esté de 1 el RTS score, más relevante será. 
root@50c24249d0ba:/home/eqtl/eQTL# awk '$20=="1"' result/rtc.txt
rs909685 rs909685 ENSG00000100321.10 ENSG00000100321.10 22 39747671 0 22 39747671 0 22 39745930 0 1741 65647 65647 39660501 39757500 243 1 1 1

Variant ID 	Location 	Variant type 			Gene 	Molecular consequences 	1000G MAF
				
rs909685 	39,351,666 	single nucleotide variant 	SYNGR1 	intron variant 		T = 0.4938 			

Alleles associated with rs909685
Allele information 
Variant allele 	Transcript change 	RefSeq 	          	Molecular consequence
A 	        c.99+1557T>A 	        NM_004711.4 		intron variant 					
A 		c.99+1557T>A 		NM_145731.3 		intron variant

#

#Q3: Which consequences, according to the variant effect predictor, do these variants have?
#100% intergenic variant
#http://grch37.ensembl.org/Homo_sapiens/Tools/VEP/Results?db=core;tl=9C5BKm6odJBkr81z-4725536

#Prepare input files for CAVIAR:

# Generate the ID/Z-scores input file. Select your favourite gene (e.g. gene=ENS00000000000.0).
# Set k (number of variants) to 50
root@e308405d0d6d:/home/eqtl/eQTL# gene=ENSG00000198951.6
root@e308405d0d6d:/home/eqtl/eQTL# compZscore.R --gene $gene --nominal result/nominals.txt -k 50 --output tmp/$gene.rs_z

# Generate the LD matrix 
root@e308405d0d6d:/home/eqtl/eQTL# plink --r square --snps $(cut -f1 tmp/$gene.rs_z) --vcf input/processed/genotypes.chr22.vcf.gz --out tmp/$gene

#Task 17: Obtain a credible set of causal variants with probability ρ=0.95.
root@e308405d0d6d:/home/eqtl/eQTL# CAVIAR -z tmp/$gene.rs_z -l tmp/$gene.ld -o result/$gene

#Explore the output files
root@50c24249d0ba:/home/eqtl/eQTL# head result/ENSG00000198951.6.log
2.02474e+50
2.02474e+50
2.02474e+50

root@50c24249d0ba:/home/eqtl/eQTL# head result/ENSG00000198951.6_post
SNP_ID	Prob_in_pCausalSet	Causal_Post._Prob.
rs133335	0.0415479	0.0830959
rs133339	0.438681	0.877362
rs112466745	0.5	0.999999
rs1800754	0.0197637	0.0395274
rs134876	7.13652e-06	1.4273e-05
rs5758670	2.87375e-08	5.7475e-08
rs5751250	1.55829e-08	3.11658e-08
rs134902	5.79104e-09	1.15821e-08
rs36093924	5.72266e-09	1.14453e-08

root@50c24249d0ba:/home/eqtl/eQTL# head result/ENSG00000198951.6_set 
rs133335
rs133339
rs112466745

#Q1: How many variants are there in the credible (ρ=0.95) set? For each of these variants, which is the probability to be causal? 
root@e308405d0d6d:/home/eqtl/eQTL# wc -l result/ENSG00000198951.6_set
3 result/ENSG00000198951.6_set
#rs133335 
root@e308405d0d6d:/home/eqtl/eQTL# awk '$1=="rs133335"' result/ENSG00000198951.6_post
rs133335	0.0415479	0.0830959
#rs133339
root@e308405d0d6d:/home/eqtl/eQTL# awk '$1=="rs133339"' result/ENSG00000198951.6_post
rs133339	0.438681	0.877362
#rs112466745 
root@d6add0c3151e:/home/eqtl/eQTL# awk '$1=="rs112466745"' result/ENSG00000198951.6_post
rs112466745	0.5	0.999999

#Q2: Which are the p-values and effect sizes of these variants? 
#rs133335
root@5ebcf2ed1aac:/home/eqtl/eQTL# awk '{print $1, $12, $13}' <(awk '$8=="rs133335"' <(awk '$1=="ENSG00000198951.6"' result/nominals.txt))
ENSG00000198951.6 7.02635e-08 0.223708

#rs133339
root@5ebcf2ed1aac:/home/eqtl/eQTL# awk '{print $1, $12, $13}' <(awk '$8=="rs133339"' <(awk '$1=="ENSG00000198951.6"' result/nominals.txt))
ENSG00000198951.6 7.02635e-08 0.223708

#rs112466745 
root@5ebcf2ed1aac:/home/eqtl/eQTL# awk '{print $1, $12, $13}' <(awk '$8=="rs112466745"' <(awk '$1=="ENSG00000198951.6"' result/nominals.txt))
ENSG00000198951.6 2.08494e-07 -0.209815


#How are they in comparison to the p-values and effect sizes of other variants tested for the same gene? 
#The p-value of our varients are smaller and the effect size is bigger. 
root@5ebcf2ed1aac:/home/eqtl/eQTL# awk '{print $1, $12, $13}' <(awk '$1=="ENSG00000198951.6"' result/nominals.txt)

#Q3: Which consequences, according to the variant effect predictor, do these variants have?
#intron_variant: 38% / non_coding_transcript_variant: 26% / NMD_transcript_variant: 17% / missense_variant: 11% / upstream_gene_variant: 4% / non_coding_transcript_exon_variant: 2% / 3_prime_UTR_variant: 2%
#http://grch37.ensembl.org/Homo_sapiens/Tools/VEP/Results?db=core;tl=IyrEm821GXgh4IAk-4705162

#Task 18: Generate a LocusZoom plot. Use as 'SNP' any of the colocalized or fine-mapped variants. 
# Define the gene corresponding to the co-localized or fine-mapped variants of interest
root@5ebcf2ed1aac:/home/eqtl/eQTL# gene=ENSG00000198951.6
root@5ebcf2ed1aac:/home/eqtl/eQTL# cat <(echo "MarkerName P.value") <(grep $gene result/nominals.txt | cut -d " " -f8,12) > tmp/metal.$gene
#plot: /home/eqtl/eQTL/result/plots/EUR.rs35696419.400kb.pdf
