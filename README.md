# HeteroGenesis
## Introduction
HeteroGenesis is used to generate genomes for multiple related clones in a heterogeneous tumour, along with a matched germline genome. For each clone and germline sample, it provides FASTA files containing the sequences for each copy of a chromosome in the genome, and files detailing the variants incorporated. 

HeteroGenesis can also then be used to combine the variant profiles outputs of each clone to give overall bulk tumour outputs that reflect user defined proportions of each clone in a tumour, and its purity. This is useful, for example, when the user intends to carry out *in silico* sequencing of each clone and combine the reads to form a bulk tumour dataset.

For more information, see "Simulation of Heterogeneous Tumour Genomes with HeteroGenesis and In Silico Whole Exome Sequencing, Tanner G et al., 2019" - https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty1063/5273483 . 
Please cite this when using HeteroGenesis.

## Versions
v1.5 - Improved speed. Allows single chromosome to be processed at a time in varincorp. Bug fixed where deletion indels were allowed to overlap other snvs/indels. (17/02/20)

v1.4 - Known variants taken from vcf file instead of flat files and improved speed when using known variants. Fixed default germline SNV rate and somatic CNV number. (02/01/20)

v1.3 - Bug fixed in freqcalc that caused variant allele frequencies to be calculated incorrectly. (26/03/19)

v1.2 - Allows the user to give lists of SNVs, indels, or CNVs for germline or somatic variants to be taken from. (12/02/19)

v1.1 - Release at paper acceptance. (21/12/18)

## Requirements

Python3 and numpy are required to run HeteroGenesis. Python 3.5.2 and numpy 1.12.0 and 1.12.1 have been tested succesfully with it.

**heterogenesis\_vargen** takes ~6hrs and 5GB RAM on a single thread to run under default parameters, which includes a germline and 2 somatic clones. 
**heterogenesis\_varincorp** takes ~6hr and 8GB RAM on a single thread to process chr1 of ’clone1’ of this output. This can be run in parallel for all chromosomes and clones.


## Installation

```
git clone https://github.com/GeorgetteTanner/HeteroGenesis.git
cd HeteroGenesis/
python setup.py install
```
 
## Overview

HeteroGenesis is implemented in three parts: 

The first, **heterogenesis\_vargen**, takes: i) a FASTA genome sequence, ii) a .fai index file for the genome sequence, iii) an optional file containing known germline SNV and InDel locations and minor allele frequencies from dbSNP, and iv) a JSON file containing a set of parameters. It outputs a JSON file with lists of variants for the germline and each clone in the simulated tumour, as well as files containing the order that mutations occurred in each.

The second part, **heterogenesis\_varincorp** is then run, once for each clone, and incorporates the list of variants for a clone into a reference genome. It outputs: i) the FASTA genome sequence (one file for each copy of a chromosome), ii) a VCF file of SNV and InDel positions and frequencies, and iii) a file containing the copy numbers along the genome.

The last part, **freqcalc**, can then be run to combine outputs from clones to generate bulk tumour outputs. It takes a file containing the proportions of each clone in a tumour, along with the outputs from heterogenesis\_varincorp, and outputs equivalent files for the bulk tumour.


## Implementation

### heterogenesis_vargen 

```
heterogenesis_vargen -j example.json

```
-v/--version : Version 

-j/--json : JSON file containing parameters. 

### heterogenesis_varincorp

```
heterogenesis_varincorp -j example.json -c clone

```
-v/--version : Version 

-j/--json : JSON file containing parameters (the same file used for heterogenesis_vargen).  

-c/--clone : Name of clone to generate genomes for.

-x/--chromosome : Optional - Name of a single chromosome to process. Output files can be combined for multiple chromosomes after running.

### freqcalc


```
freqcalc -c clones.txt -d {directory of heterogenesis_varincorp outputs} -p {prefix} -n {name}

```
-v/-—version : Version 

-c/--clones : File with clone proportions in format: 'clone name' \t 'fraction’.

-d/--directory : Directory containing outputs of heterogenesis\_varincorp. This should be the same as what was provided for the ‘directory’ parameter with heterogenesis_varincorp.

-p/--prefix : Prefix of heterogenesis\_varincorp output file names. This should be the same as what was provided for the ‘prefix’ parameter with heterogenesis_varincorp.

(If the -x option was used in varincorp to process individual chromosoms separately, the vcf an cnv output files must be combined and chromosome names removed from file names before running freqcalc.)

## Inputs

### heterogenesis_vargen

1. **Reference Genome:**
The starting genome sequence, in FASTA format, that variants will be incorporated into. 
2. **Reference Genome Index:**
A .fai index file for the reference genome, created with samtools faidx. This should be saved in the same directory as the reference genome.
3. **dbSNP vcf File:**
A vcf file of known germline SNVs and InDels from dbSNP (uncompressed). Eg. ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz. This may be filtered for lines containing "CAF" and subsampled to around 20,000,000 lines to reduce disk space or memory requirements if necessary. Fewer lines than this may be used but that will likely start to reduce the effect of more common known SNPs being incorporated more frequently than rarer known SNPs.

4. **Parameters File:**
A JSON file containing run parameters and locations of other inputs. Any parameter that is missing from the file will be set at its default value:
 
	(An example parameters file is provided in the repository - 'example.json')

|Parameter|Description|Default Value| 
|---|---|---|
|prefix	|String added to output file names.|""|
|reference|FASTA file containing the sequence of a reference or other input genome. Must have a .fai index file located in the same directory. |Required|
|dbsnp|A vcf file of known germline SNPs and InDels from dbSNP.|none|
|directory|Directory to output all files to.|"./"|
|structure|Structure of clones in the tumour, in the format: “clone1\_name, clone1\_distance\_from\_parent, clone1\_parent\_name, clone2_name, clone2\_distance\_from\_parent, clone2\_parent\_name…”. All parent clone names must also be listed as a separate clone, ie. if clone2’s parent clone is clone1, then clone1 must also be listed as a clone with a parent clone. The exception to this is when the parent clone is ‘germline’, and this must occur at least once as the parent clone for the root clone of the tumour. Loops in the lineage will cause the program to never end, ie. clone1->clone2->clone3->clone1. Distances from parent clones can be any fraction or number as they are used relative to each other.|"clone1,0.2,germline,clone2,0.8,clone1"|
|snvgermline|Rate of germline SNVs per base.|0.0014|
|indgermline|Rate of germline indels per base|0.00014|
|cnvrepgermline|Number of germline replication CNVs.|160|
|cnvdelgermline|Number of germline deletion CNVs.|1000|
|aneuploid|Number of somatic aneuploid events. i.e. replication or deletion of chromosomes. These can either be whole genome duplication or individual chromosome duplication or deletion. Aneuploid events are prevented from deleting all copies of a chromosome. Germline aneuploid events are not available.|2|
|wgdprob|Probability that each aneuploid event is a whole-genome duplication.|0.0|
|snvsomatic|Rate of somatic SNVs per base.|0.00001|
|indsomatic|Rate of somatic indels per base.|0.000002|
|cnvrepsomatic|Number of somatic replication CNVs.|250|
|cnvdelsomatic|Number of somatic deletion CNVs.|250|
|dbsnpsnvproportion|Proportion of germline SNVs taken from dbSNP. |0.99|
|dbsnpindelproportion|Proportion of germline InDels taken from dbSNP. |0.97|
|chromosomes|List of chromosomes to include in the model. Alternatively, "all" can be given, in which case chromosomes 1-22 will be used. This only works for genomes for which chromosomes are labelled 'chr1','chr2'... (Also note that X and Y are not included with "all")|”all”|
|givengermlinesnvs, givengermlineindels, givengermlinecnvs, givensomaticsnvs, givensomaticindels, givensomaticcnvs|Used to provide lists of variants for when the user wishes to sample from given variants instead of randomly generating them. See 'examplegivenXXX.txt' files for formatting. Only variants that fit into the genome (eg. not in deleted regions etc. will be used). Note: when a CNV is sampled from a given list, the distinction between replication and deletion CNVs (eg. cnvrepsomatic vs cnvdelsomatic) is ignored and the copy number is instead just taken from the given list. |''|
|givengermlinesnvsproportion, givengermlineindelsproportion, givengermlinecnvsproportion, givensomaticsnvsproportion, givensomaticindelsproportion, givensomaticcnvsproportion|The proportion of variants taken from given lists. For germline SNVs/InDels the proportion of randomly generated variants is 1-(dbsnpsnvproportion + givengermlinesnvsproportion).|0.0|

CNV lengths and copy numbers, and indel lengths are taken from lognormal distributions, that are defined by the mean and variance of the underlying normal distribution. Values from these distributions are then scaled up by a multiplication factor for cnv lengths. Indel length distributions are the same for germline and somatic.

| | | |
|---|---|---|
|cnvgermlinemean|Germline CNV length lognormal mean.|-10|
|cnvgermlinevariance|Germline CNV length lognormal variance|3|
|cnvgermlinemultiply|Germline CNV length multiplication factor.|1000000|
|cnvsomaticmean	|Somatic CNV length lognormal mean.|-1|
|cnvsomaticvariance|Somatic CNV length lognormal variance.|3|
|cnvsomaticmultiply|Somatic CNV length multiplication factor.|1000000|
|indmean|Indel length lognormal mean.|-2|
|indvariance|Indel length lognormal variance.|2|
|indmultiply|Indel length multiplication factor.|1|
|cnvcopiesmean|CNV copies lognormal mean.|1|
|cnvcopiesvariance|CNV copies lognormal variance.|0.5|


### heterogenesis_varincorp

1. **Variants File:** From heterogenesis\_vargen. 

2. **Parameters File:**
The same JSON file as used for heterogenesis\_vargen can be given but only the following parameters are used. These should contain the same values as given for heterogenesis\_vargen: 

|Parameter|Description|Default Value| 
|---|---|---|
|prefix	|String added to output file names.|""|
|reference|FASTA file containing the sequence of a reference or other input genome. Must have a .fai index file located in the same directory. |Required|
|directory|Directory containing the JSON variants file output from heterogenesis\_vargen and where output files will be written to.|"./"|
|chromosomes|List of chromosomes included in the model. Alternatively, "all" can be given, in which case chromosomes 1-22 will be used. This only works for genomes for which chromosomes are labelled 'chr1','chr2'... (Also note that X and Y are not included with "all")|”all”|

### freqcalc

1. **Clones File:**
File with clone proportions in the format: 'clone name' \t 'fraction’ \n.

2. **Outputs From heterogenesi\_varincorp**

## Outputs

### heterogenesis_vargen
1. **_prefix_varaints.json:** A JSON file containing information from a python dictionary in the format: [clone][chromosome][variants, SNV/InDel positions, CNV breakpoints, deleted regions]. This is for use by heterogenesis_varincorp and not intended to be manulally viewed.
2. **_prefixcloneX_variants.txt:** This file lists every variant that occured in the clone. 

### heterogenesis_varincorp
1. **{prefix}{cloneX}cnv.txt:** This records the copy number status along the genome, allong with phased major/minor alleles.(Positions are 1 based.)

2. **{prefix}{cloneX}.vcf:** This records the position and variant allele frequency (VAF) for each SNV/InDel, allong with the number of occurences on each copy of a chromosome and the overall copy number at that position.

3. **{prefix}{cloneX}{chrXX}.fasta:** The genome sequence in FASTA format. (One file for each copy of each chromosome.)

### freqcalc
1. **{prefix}{sample}cnv.txt:** This records the combined copy number status along the genome, allong with phased major/minor alleles, for the bulk tumour sample.(Positions are 1 based.)
 
2. **{prefix}{sample}.vcf:** This records the combined position and variant allele frequency (VAF) for each SNV/InDel, allong with the phasing of each variant (ie. if the variant occured on an A or B copy of a chromosome) and the overall copy number at that position, for a bulk tumour sample.


## Example


This example demonstrates how to run the entire HeteroGenesis process on test parameters. This limits the simulation to use only chromosomes 21 and 22 in order to reduce run time to a few minutes. For full runs, most parameters can be deleted from the json files and instead ran as defaults.

```bash
#Install:
git clone https://github.com/GeorgetteTanner/HeteroGenesis.git
cd HeteroGenesis/
python setup.py install

#EITHER:

#1. If wanting germline variants from dbSNP, download a dbsnp file and filter it to reduce memory requirement. 
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz
gunzip dbsnp_146.hg38.vcf.gz | grep "CAF" | gshuf -n 20000000 > dbsnp.hg38.vcf

#OR:

#2. If not wanting to include dbSNP variants, remove the 
#'"dbsnp":"./dbsnp.hg38.vcf",' line from example.json
awk '!/dbsnp/' example.json > temp ; mv temp example.json


#Download reference genome (or copy from locally saved reference genome to save time):
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai
gunzip Homo_sapiens_assembly38.fasta.gz

#make test directory:
mkdir ../test1
cd ../test1

#Run heterogenesis_vargen: ~1min
heterogenesis_vargen -j ../HeteroGenesis/example.json

#Run heterogenesis_varincorp on each clone and germline: ~5min
for clone in clone1 clone2 germline ; do heterogenesis_varincorp -j ../HeteroGenesis/example.json -c ${clone} ; done

#Run freqcalc to create bulk sample variant profiles: 
freqcalc -c ../HeteroGenesis/example_clones.txt -d . -p test1 -n sample1

```

If you want to carry out _in silico_ whole-exome sequencing of the created tumour, this can be achieved with the following:

(See https://github.com/GeorgetteTanner/w-Wessim for further details.)

This example demonstrates how to use w-Wessim for _in silico_ whole-exome sequencing with a "real reads" probe sequences set. This method creates the most realistic sequencing dataset, but is a very time and memory consuming process and not feasible without access to a high performance computing system. Therefore the option of using a subsampled probe set (with 1/1000th of the probes) is available if the user wishes to run the programs on a standard computer for testing. The resulting sequencing data set from this will look very patchy and is not intended for any use. Instructions for both options are included below.

Alternatively, the probe sequences from an exon capture kit can be used. This is much quicker and less memory intensive but results in a less realistic distribution of reads. Instructions for this can be found at https://github.com/GeorgetteTanner/w-Wessim.



```bash
#Download programs - (you may get a few warnings during pblat installation that can be ignored):
cd ..
git clone https://github.com/GeorgetteTanner/w-Wessim.git
git clone https://github.com/icebert/pblat.git
cd pblat
make
cd ..
#(need to find the correct binary file for your operating system:)
wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/faToTwoBit
chmod a+x ./fatotwobit

#EITHER:

#1. Download full set of probe sequences and convert from fastq to fasta:
cd w-Wessim
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA157/ERA1574375/fastq/real_wes_reads_probes.fastq.gz
gunzip real_wes_reads_probes.fastq.gz
paste - - - - < real_wes_reads_probes.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > real_wes_reads_probes.fa

#OR

#2.Download subsampled probe sequences:
cd w-Wessim 
wget https://github.com/GeorgetteTanner/data/raw/master/real_wes_reads_probes_subsampled.fa.gz
gunzip real_wes_reads_probes_subsampled.fa.gz
mv real_wes_reads_probes_subsampled.fa real_wes_reads_probes.fa

#Combine fasta files. If using a whole genome run, these need to be
#grouped into files of no more than around 2GB to avoid errors with 
#pblat - we recommend grouping into 4. When using only a couple of 
#chromosomes, all copies can be grouped together but for 
#demostration purposes we group into 2 files here. Alternatively, 
#code for grouping files into 4 is given below. (More than 4 files 
#will probably be needed on whole genome runs with whole-genome 
#duplication events.) Each run of pblat takes a similar length of 
#time regardless of genome length so its best to group into as few 
#files as possible.: ~2sec

#EITHER:

#1.
cd ../test1
for clo in clone1 clone2 germline ; do cat test1${clo}chr*A*fasta > test1${clo}_1.fasta ; done
for clo in clone1 clone2 germline ; do cat test1${clo}chr*B*fasta > test1${clo}_2.fasta ; done

#OR:

#2.
cd ../test1
for clo in clone1 clone2 germline ; do cat test1${clo}chr1*A*fasta > test1${clo}_1.fasta ; done
for clo in clone1 clone2 germline ; do cat test1${clo}chr1*B*fasta > test1${clo}_2.fasta ; done
for clo in clone1 clone2 germline ; do cat test1${clo}chr[^1]*A*fasta > test1${clo}_3.fasta ; done
for clo in clone1 clone2 germline ; do cat test1${clo}chr[^1]*B*fasta > test1${clo}_4.fasta ; done

#Convert each chromosome fasta to 2bit: 
for f in test1*_*.fasta ; do faToTwoBit $f ${f}.2bit ; done

#pblat: ~1h for subsampled probes
for f in test1*.fasta.2bit ; do ../pblat/pblat $f ../w-Wessim/real_wes_reads_probes.fa $(basename $f .fasta.2bit).psl -minScore=95 -minIdentity=95 -threads=8; done

#Combine .psl files:
#Save the header:
head -n 5 $(ls *.psl | head -1 ) > pslheader.txt
#Remove the headers:
for f in *.psl ; do tail -n+6 $f > noheader_$f ; done 
#Combine noheader*.psl files and sort combined file on column 10: 
for clone in clone1 clone2 germline ; do cat noheader_test1${clone}*.psl | sort -k 10 -n > sorted_combined_noheader_test1${clone}.psl ; done
#Add the header:
for clone in clone1 clone2 germline ; do cat pslheader.txt sorted_combined_noheader_test1${clone}.psl > ${clone}.psl ; done

#Combine fasta files into full genomes
for clone in clone1 clone2 germline ; do cat test1${clone}*.fasta > test1${clone}.fasta ; done

#w-Wessim:
cd ../w-Wessim
for clo in clone1 clone2 germline ; do python2 w-wessim.py -R ../test1/test1${clo}.fasta -P real_wes_reads_probes.fa -B ${clo}.psl -n 100000 -l d -M lib/hs2000p.gzip -pz -o ../test1/w-wessimoutput_{clo} -t 1 -v -m 20 -f 170 -d 35 ; done

#Align reads using own pipelines. (Cleaning reads is not recomended 
#if using the error model provided with w-Wessim as that was 
#trained on a pre-cleaned dataset to improve alignment accuracy of 
#the training set.)

#Subsample BAM files in proportions listed in the clones file:
samtools view -b -h -s 0.15 -o germline_0.20.bam germline.bam 
samtools view -b -h -s 0.65 -o clone2_0.65.bam clone2.bam 
samtools view -b -h -s 0.15 -o clone1_0.15.bam clone1.bam 

#Combine subsampled bam files to create bulk data.
samtools merge -c -p example_bulk.bam germline_0.20.bam clone2_0.65.bam clone1_0.15.bam
```