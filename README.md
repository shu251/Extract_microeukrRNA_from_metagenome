## Extracting 18S/28S rRNA gene signatures from metagenomic data
*Last updated 11/02/2018*
### Overview
Meta'omic approaches to explore microbial eukaryotes *in situ* typically lag behind prokaryotes. This is largely due to the complexity of eukaryotic genomes and the lack of cultured representatives. However, we are closer than ever to being able to more rapidly sequencing and analyze microbial eukaryotic metatranscriptomes and metagenomes.

Here, I've outlined a pipeline to probe metagenome data (shot gun sequencing of all DNA collected), originally intended for bacterial metagenomic work, for signs of eukaryotic life.

### Workflow Summary & Requirements
1) Trim and QC raw reads
    * Raw fastq reads!
    * [Trimmomatic v0.33](http://www.usadellab.org/cms/?page=trimmomatic)
2) Sort reads (select for reads that align to rRNA databases)
    * [SortMeRNA v2.1](http://bioinfo.lifl.fr/RNA/sortmerna/)
3) Merge paired end reads
    * [fastq-join](http://qiime.org/scripts/join_paired_ends.html) in QIIMEv1 or on its [own](https://expressionanalysis.github.io/ea-utils/)
4) Assign taxonomy
    * [QIIMEv1](http://qiime.org/scripts/assign_taxonomy.html) or another taxonomy assignment pipeline of your choice (i.e. mothur, blast)
5) Compile output
    * [R](https://www.r-project.org/)
    * [jupyter notebook](http://jupyter.org/) to run R


#### 1) Trim and QC raw reads
Supply Trimmomatic with your raw reads, both R1 and R2 and then name the outputs, which will be those that passed the QC. 'Unpaired' output reads will be those that will not pair with the other R1/2.

```
java -jar [PATHtoTrimmomatic]/Trimmomatic-0.33/trimmomatic-0.33.jar PE -phred33 [read_1.fastq.gz] [read_2.fastq.gz] [read_1_paired.fastq.gz] [read_1_unpaired.fastq.gz] [read_2_paired.fastq.gz] [read_2_unpaired.fastq.gz] ILLUMINACLIP:TruSeq3-PE.fa:2:30:20 LEADING:5 TRAILING:5 SLIDINGWINDOW:25:10 MINLEN:50
```
Use the 'paired' output from Trimmomatic

#### 2) Sort reads (select for reads that align to rRNA databases)

First, make sure SortMeRNA is downloaded properly and you have unpacked the databases it comes with. I'm working with version 2.1 now - future versions will also probably work.

Paired reads from Trimmomatic, R1_paired and R2_paired need to to be interleaved. SortMeRNA installs a script to do this called 'merge-paired-reads.sh'.

```
# Script located: ~/sortmerna-2.1/scripts/
# usage: bash merge-paired-reads.sh file1.fastq file2.fastq outputfile.fastq

# May need to unzip from trimmomatic output if you used fastq.gz input/output read files

bash merge-paired-reads.sh [read_1_paired.fastq] [read_2_paired.fastq] [outputfile_merged.fastq]
```

Next, you want to take the interleaved paired end reads and sort for rRNA reads that align to whichever databases you want (that SortMeRNA supplies).

Ahead of time, make sure you've got the database you want in **~/sortmerna-2.1/rRNA_databases/** and the associated index that you've built here: **~/sortmerna-2.1/index/**
rRNA databases and index files are listed one after another (if you want to use more than 1) and separated by **":\\"**

To sort reads using one database (16S), list the location of the rRNA database and then a comma, and then the location of the 16S index
```
~/sortmerna-2.1/sortmerna --ref ~/sortmerna-2.1/rRNA_databases/silva-arc-16s-id95.fasta,~/sortmerna-2.1/index/silva-arc-16s-db --reads [outputfile_merged.fastq] --aligned [merged_rrna] --other [merged_nonrrna] --log --fastx -a 6
# -a flag specifies running in parallel, read sortmerna manual for more information
```

To sort reads using all the databases available (my preferred method):
```
~/sortmerna-2.1/sortmerna --ref ~/sortmerna-2.1/rRNA_databases/silva-arc-16s-id95.fasta,~/sortmerna-2.1/index/silva-arc-16s-db:\~/sortmerna-2.1/rRNA_databases/silva-arc-23s-id98.fasta,~/sortmerna-2.1/index/silva-arc-23s-db:\~/sortmerna-2.1/rRNA_databases/silva-bac-16s-id90.fasta,~/sortmerna-2.1/index/silva-bac-16s-db:\~/sortmerna-2.1/rRNA_databases/silva-bac-23s-id98.fasta,~/sortmerna-2.1/index/silva-bac-23s-db:\~/sortmerna-2.1/rRNA_databases/silva-euk-18s-id95.fasta,~/sortmerna-2.1/index/silva-euk-18s-db:\~/sortmerna-2.1/rRNA_databases/silva-euk-28s-id98.fasta,~/sortmerna-2.1/index/silva-euk-28s-db --reads [outputfile_merged.fastq] --aligned [merged_rrna] --other [merged_nonrrna] --log --fastx -a 6
```
*See output log file from SortMeRNA. At the end of this log file, it lists the percent of your input reads that aligned to each rRNA database

After SortMeRNA is run, interleaved reads need to be unmerged. Using a similar script to the one used above to merge the read pairs.

```
bash ~/sortmerna-2.1/scripts/unmerge-paired-reads.sh [merged_rrna] [R1_rrna] [R2_rrna]
```

The **R1_** and **R2_rrna** products are paired end reads which were identified to align to rRNA databases.

#### 3) Merge paired end reads
Taxonomy assignment will be much improved if paired end reads are merged together. There are other programs to do this, but I'm using a script that comes with QIIME v1 (fastq-join).

```
join_paired_ends.py -f R1_rRNA.fastq -r R2_rRNA.fastq -o joined
```
Merged reads are still in fastq format. Depending on how you plan to perform the rest of your workflow, you may want to continue conducting quality control or, here, I'm going to convert to fasta format so I can go directly to taxonomy assignment.
```
# output from join_paired_ends.py is placed in an output directory.
## Be careful to not overwrite your files if you're processing more than one sample at a time
convert_fastaqual_fastq.py -f joined/fastqjoin.join.fastq -c fastq_to_fastaqual -o fasta_dir
# again output is placed in a directory
# Here I move and rename output merged file so it reflects the sample ID
mv fasta_dir/fastqjoin.join.fna [sample_rRNA_joined.fna]

```

***Alternative workflow***
To continue increasing read length, you can also use emirge. [See Ben Tully's approach here](https://www.protocols.io/view/Detecting-16S-rRNA-Gene-Fragments-from-a-Metagenom-d7u9nv). As expected, the use of emirge for exploring microbial eukaryotes was not found to be helpful. There is too much diversity and not enough reads to do successfully get longer reads from emirge. This means we need to be careful about how we interpret diversity.

#### 4) Assign taxonomy

Here I'm using QIIMEv1 to run a uclust search against the SILVA v132 database. An alternative is to use [PR2 for protists](https://github.com/vaulot/pr2database). But here, I found many reads to hit both the 18S and 28S rRNA databases, so I wanted to use both the Silva SSU and LSU databases.

To format the SILVA database ahead of time, I followed the [instructions from here](https://github.com/mikerobeson/Misc_Code/tree/master/SILVA_to_RDP).

```
# LSU database:
assign_taxonomy.py -i sample_rRNA_joined.fna -t ~/db/SILVA_132_LSURef_99qiime1.tax -r ~/db/SILVA_132_LSURef_99qiime1.fasta -m uclust -o taxassign_uclust/sample_taxLSU
# SSU database:
assign_taxonomy.py -i sample_rRNA_joined.fna -t ~/db/SILVA_132_SSURef_99qiime1.tax -r ~/db/SILVA_132_SSURef_99qiime1.fasta -m uclust -o taxassign_uclust/sample_taxSSU\n";

# Move and change name of taxonomy files
mv taxassign_uclust/sample_taxLSU/*.txt taxassign_uclust/sample_LSU_tax_assignment.txt

mv taxassign_uclust/sample_taxSSU/*.txt taxassign_uclust/sample_SSU_tax_assignment.txt
```

#### 5) Compile output

See R script in github repo. **"Compile_tax_assignments.ipynb"**
Script compiles series of output "sample_LSU_tax_assignment.txt" and "sample_SSU_tax_assignment.txt" files and generates a count table. Then visually explores the eukaryotic diversity in each sample.