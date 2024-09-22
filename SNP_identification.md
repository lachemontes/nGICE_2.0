# SNP identification molestus vs pipiens 

Path to the the project folder ` /proj/snic2022-23-541/Rohan` in Rackham

## Conda env for Rohan

```bash
# Create an enviroment to install the programs
conda create --name Rohan

# activate the enviroment
conda activate Rohan
```

## Install Programs via conda

### BWA, Samtools,BCF tools, GATK4

```bash
# Install BWA
conda install bioconda::bwa
#Install Picard to remove remove PCR duplicates
conda install bioconda::picard
# Install BCF
conda install bioconda::bcftools
#libgsl.so.25 as requeried to run bcftools
conda install -c conda-forge gsl
# Install GATK4
conda install bioconda::gatk4
# Install samtools
conda install bioconda::samtools
# Install VCF tools
conda install bioconda::vcftools

#
conda install conda-forge::py-bgzip


```

## Data pre-processing

```bash
# remove everthing after the first speace in genome fasta file to make more simple to read

# Save original header in a txt file
grep -e ">" VectorBase-66_CquinquefasciatusJHB2020_Genome.fasta > Headers_genome.txt

# Remove everything in the header of a fasta after the first space:
awk '{print $1;next}1' file.fa > output.fa
```



## BWA to map read to refGenome

#### Index genome

```bash
#!/bin/bash
#SBATCH -A naiss2023-5-461
#SBATCH -p core n-8
#SBATCH -t 12:00:00
#SBATCH -J BWA_Index
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se
#SBATCH --output=hisat_%A_%a.out
#SBATCH --error=hisat_%A_%a.err

# Build and Index before mapping

bwa index ../../Data/Genome/VectorBase-66_CquinquefasciatusJHB2020_Genome_headers.fasta


```

#### Run main BWA

```bash
#!/bin/bash
#SBATCH -A naiss2023-5-461
#SBATCH -p core -n 12
#SBATCH -t 12:00:00
#SBATCH -J BWA_Mapping
#SBATCH --array=1-14
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se
#SBATCH --output=bwa_%A_%a.out
#SBATCH --error=bwa_%A_%a.err


# The programs are installed in Rohan env

# Define input directories

input_dir="/proj/snic2022-23-541/Rohan/Data/RNAseq"
output_dir="/proj/snic2022-23-541/Rohan/Analysis/BWA"
genome_index="/proj/snic2022-23-541/Rohan/Data/Genome/Index/VectorBase-66_CquinquefasciatusJHB2020_Genome_headers.fasta"


# Define sample name
sample=$(sed -n "${SLURM_ARRAY_TASK_ID} p" List.txt)

# Define input files
r1="${input_dir}/${sample}_R1.fastq"
r2="${input_dir}/${sample}_R2.fastq"


# Define output file
output="${output_dir}/${sample}.sam"

# Run BWA for Illumina/454/IonTorrent paired-end reads longer than ~70bp:

bwa mem "${genome_index}" "${r1}" "${r2}" -t 8 > "${output}"

```

```bash
  bwa mem -t 8 "${genome_index}" "${r1}" -2 "${r2}" > "${output}"

```

### From Sam to Bam 


```bash
#!/bin/bash
#SBATCH -A naiss2023-5-461
#SBATCH -p core -n 12
#SBATCH -t 4-00:00:00
#SBATCH -J SAMToBAM
#SBATCH --mail-type=All
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se
#SBATCH --array=1-7

# Load necessary modules or set necessary environment variables here if needed (we have them in the conda env)

filename=$(sed -n "${SLURM_ARRAY_TASK_ID}p" List_stb.txt)
input_file="/proj/snic2022-23-541/Rohan/Analysis/BWA/${filename}.sam"
output_file="/proj/snic2022-23-541/Rohan/Analysis/BAM/${filename}.bam"
sorted_output_file="/proj/snic2022-23-541/Rohan/Analysis/BAM/${filename}.sorted.bam"

# Convert SAM to BAM
samtools view -bS $input_file > $output_file

# Sort the BAM file
samtools sort $output_file -o $sorted_output_file

# Index the sorted BAM file if needed
samtools index $sorted_output_file

# statistics about the sorted bam file

samtools flagstat $sorted_output_file

# List all the sorted BAM files
bam_files=$(ls *.sorted.bam)

# Merge the BAM files
samtools merge merged_output.bam $bam_files

# Index the merged BAM file
samtools index merged_output.bam

```

### MarkDuplicates from Picard to remove PCR duplicates

Sometimes the same DNA fragment is sequenced multiple times, which leads to multiple reads from the same fragment in the fastq file. This can occur due to PCR amplification in the library preparation, or if one read cluster is incorrectly detected as multiple clusters by the sequencing instrument. If a duplicated read contains a genetic variant, the ratio of the two alleles might be obscured, which can lead to incorrect genotyping. It is therefore recommended (in most cases) to mark duplicate reads so that they are counted as one during genotyping.
from `https://nbisweden.github.io/workshop-ngsintro/2403/topics/vc/lab_vc.html#mark-duplicates`



```bash
#!/bin/bash
#SBATCH -A naiss2023-5-461
#SBATCH -p core -n 12
#SBATCH -t 2:00:00
#SBATCH -J MarkDuplicates
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se
#SBATCH --output=mark_duplicates_%j.out
#SBATCH --error=mark_duplicates_%j.err
#SBATCH --array=1-7


# Define input and output files

filename=$(sed -n "${SLURM_ARRAY_TASK_ID}p" List_md.txt)
input_file="/proj/snic2022-23-541/Rohan/Analysis/BAM/${filename}.bam"
output_file="/proj/snic2022-23-541/Rohan/Analysis/MarkDuplicates/${filename}.md.bam"
metrics_file="/proj/snic2022-23-541/Rohan/Analysis/MarkDuplicates/${filename}.md.metrics.txt"


# Run Picard MarkDuplicates
gatk --java-options -Xmx7g MarkDuplicates \
    I="${input_file}" \
    O="${output_file}" \
    M="${metrics_file}" \
    REMOVE_DUPLICATES=true \
    CREATE_INDEX=true

echo "MarkDuplicates completed successfully."




```
### Calculate the read coverage of positions in the genome (Skiped)

Do the first pass on variant calling by counting read coverage with bcftools. We will use the command mpileup. The flag -O b tells bcftools to generate a bcf format output file, -o specifies where to write the output file, and -f flags the path to the reference genome:
from `https://training.galaxyproject.org/training-material/topics/data-science/tutorials/bash-variant-calling/tutorial.html`

Example:

```bash
bcftools mpileup -O b -o raw.bcf -f ../../Data/Genome/VectorBase-66_CquinquefasciatusJHB2020_Genome.fasta --threads 2 -q 20 -Q 30 ../MarkDuplicates/Molestus_1_paired_trimmed_pairs.sorted.md.bam

````

**BCFtools** by uppmax, **via conda doesnt work

```bash
#!/bin/bash
#SBATCH -A naiss2023-5-461
#SBATCH -p core -n 12
#SBATCH -t 2-00:00:00
#SBATCH -J bcftools_mpileup
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se
#SBATCH --output=BCF_tools_%j.out
#SBATCH --error=BCF_tools_%j.err
#SBATCH --array=1-7

module load bioinfo-tools
module load bcftools/1.9

filename=$(sed -n "${SLURM_ARRAY_TASK_ID}p" List_md.txt)
input_file="/proj/snic2022-23-541/Rohan/Analysis/MarkDuplicates/${filename}.md.bam"
output_file="/proj/snic2022-23-541/Rohan/Analysis/BCFtools/${filename}.raw.bcf"

bcftools mpileup -O b -o "${output_file}" -f ../../Data/Genome/VectorBase-66_CquinquefasciatusJHB2020_Genome.fasta "${input_file}"

```

```bash

#!/bin/bash
#SBATCH -A naiss2023-5-461
#SBATCH -p core -n 12
#SBATCH -t 2-00:00:00
#SBATCH -J bcftools_mpileup
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se
#SBATCH --output=VCF_tools_%j.out
#SBATCH --error=VCF_tools_%j.err
#SBATCH --array=1-7

module load bioinfo-tools
module load bcftools/1.9

filename=$(sed -n "${SLURM_ARRAY_TASK_ID}p" List_md.txt)
input_file="/proj/snic2022-23-541/Rohan/Analysis/BCFtools/${filename}.raw.bcf"
output_file="/proj/snic2022-23-541/Rohan/Analysis/VCF_files/${filename}.vcf"
input_file2="/proj/snic2022-23-541/Rohan/Analysis/BCFtools/${filename}.vcf"
output_file_final="/proj/snic2022-23-541/Rohan/Analysis/VCF_files/${filename}_final.vcf"

bcftools call --ploidy 1 -m -v -o "${output_file}" "${input_file}"
vcfutils.pl varFilter "${input_file2} > "${output_file_final}"

```

#### Remove indels

```bash
#!/bin/bash

for vcf in *.vcf; do
	echo "Procesando el archivo $vcf"
	bcftools view -v snps -O v -o ${vcf%.vcf}_snps.vcf $vcf

done
```
# Variant calling following RNAseq short variant discovery (SNPs + Indels)

`https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indelsz`


Since we learned that RNAseq its not the best strategy for variant calling but it is still possible we decided to give it a try.
I will follow this tutorial "Tutorial: RNA-seq short variant calling using GATK4" available at `https://github.com/x-zang/tutorial-gatk4-rnaseq-germline-snps-indels`

** Some reasons why is not recommended RNAseq for variant calling**


**Coverage and Bias:** RNA-seq data is inherently biased toward the expression levels of genes. Highly expressed regions may have high coverage, while lowly expressed or non-expressed regions may have low or no coverage, leading to uneven and non-uniform coverage across the genome.

**Splicing and Transcript Variability:** RNA-seq captures only the exonic regions of the genome, missing out on intronic and intergenic variants. Additionally, alternative splicing can complicate the alignment of reads, making it challenging to determine which transcript a variant belongs to.

**Allele-Specific Expression:** Some variants may appear homozygous due to allele-specific expression when in fact they are heterozygous. This can lead to incorrect variant calling.

**Post-Transcriptional Modifications:** RNA is subject to post-transcriptional modifications such as RNA editing, which can lead to false positives in variant calling as these modifications can be mistaken for genomic variants.

**Quality of Reads:** RNA-seq libraries can have high error rates or biases due to reverse transcription, PCR amplification, and sequencing errors, which can compromise the accuracy of variant detection.




```bash
wget https://github.com/broadinstitute/gatk/releases/download/4.1.8.0/gatk-4.1.8.0.zip
unzip gatk-4.1.8.0.zip
cd gatk-4.1.8.0
```

## Mapping with STAR

### Create genome index
```bash

#!/bin/sh
#SBATCH -A naiss2023-5-461
#SBATCH -p core -n 12
#SBATCH -t 5-00:00:00
#SBATCH -J Genome_Index_node
#SBATCH --mail-type=all
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se

# Module load
module load bioinfo-tools star/2.7.9a


#Genome Index

STAR --runThreadN 60 --runMode genomeGenerate --genomeDir /proj/snic2022-23-541/Ticks_project/Analysis/Trinity/STAR/Genome_Index_Cqui --genomeFastaFiles /proj/snic2022-23-541/Rohan/Data/Genom
e/VectorBase-66_CquinquefasciatusJHB2020_Genome_headers.fasta --limitGenomeGenerateRAM 142967046410 --genomeSAindexNbases 13

```

```bash
#!/bin/sh
#SBATCH -A naiss2023-5-461
#SBATCH -p core -n 16
#SBATCH -t 4-00:00:00
#SBATCH -J STAR_Mapping
#SBATCH --mail-type=all
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se

# Module load
module load bioinfo-tools star/2.7.9a


#Genome Index


# Run the STAR command with the modified options
STAR --runThreadN 20 --genomeDir /proj/snic2022-23-541/Rohan/Analysis/STAR/Genome_Index_Cqui --readFilesIn Molestus_1_paired_trimm
ed_pairs_R1.fastq,Molestus_1_paired_trimmed_pairs_R2.fastq,Molestus_2_paired_trimmed_pairs_R1.fastq,Molestus_2_paired_trimmed_pair
s_R2.fastq,Molestus_4_paired_trimmed_pairs_R1.fastq,Molestus_4_paired_trimmed_pairs_R2.fastq,Pipiens_1_paired_trimmed_pairs_R1.fas
tq,Pipiens_1_paired_trimmed_pairs_R2.fastq,Pipiens_2_paired_trimmed_pairs_R1.fastq,Pipiens_2_paired_trimmed_pairs_R2.fastq,Pipiens
_3_paired_trimmed_pairs_R1.fastq,Pipiens_3_paired_trimmed_pairs_R2.fastq --limitBAMsortRAM 25965203113 -c  --outSAMattrRGline ID:M
olestus_1_paired_trimmed_pairs_R1.fastq, ID:Molestus_1_paired_trimmed_pairs_R2.fastq, ID:Molestus_2_paired_trimmed_pairs_R1.fastq,
 ID:Molestus_2_paired_trimmed_pairs_R2.fastq, ID:Molestus_4_paired_trimmed_pairs_R1.fastq, ID:Molestus_4_paired_trimmed_pairs_R2.f
astq, ID:Pipiens_1_paired_trimmed_pairs_R1.fastq, ID:Pipiens_1_paired_trimmed_pairs_R2.fastq, ID:Pipiens_2_paired_trimmed_pairs_R1
.fastq, ID:Pipiens_2_paired_trimmed_pairs_R2.fastq, ID:Pipiens_3_paired_trimmed_pairs_R1.fastq, ID:Pipiens_3_paired_trimmed_pairs_
R2.fastq -outSAMtype BAM SortedByCoordinate --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0.1 --outFilterMatchN
minOverLread 0.1



```

```
 finishedjobinfo -j 50307420
2024-09-17 15:14:55 jobid=50307420 jobstate=COMPLETED username=zaidemo account=naiss2023-5-461 nodes=r169 procs=16 partition=core qos=normal jobname=STAR_Mapping maxmemory_in_GiB=8.0 maxmemory_node=r169 timelimit=4-00:00:00 submit_time=2024-09-17T14:58:43 start_time=2024-09-17T15:00:23 end_time=2024-09-17T15:14:55 runtime=00:14:32 margin=3-23:45:28 queuetime=00:01:40
```

### Mapping by reads Indv

```bash
#!/bin/sh
#SBATCH -A naiss2023-5-461
#SBATCH -p core -n 16
#SBATCH -t 2-00:00:00
#SBATCH -J STAR_Mapping
#SBATCH --mail-type=all
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se

# Module load
module load bioinfo-tools star/2.7.9a

# Function to run STAR for each library
run_star() {
  local r1=$1
  local r2=$2
  local output_prefix=$3

  STAR --runThreadN 16 \
       --genomeDir /proj/snic2022-23-541/Rohan/Analysis/STAR/Genome_Index_Cqui \
       --readFilesIn $r1 $r2 \
       --outFileNamePrefix $output_prefix \
       --limitBAMsortRAM 25965203113 \
       --outSAMattrRGline ID:$output_prefix \
       --outSAMtype BAM SortedByCoordinate \
       --outFilterScoreMinOverLread 0.1 \
       --outFilterMatchNminOverLread 0.1
}

# Base path to the RNAseq reads
READS_PATH="/proj/snic2022-23-541/Rohan/Data/RNAseq"

# Call the function for each set of paired files with path
run_star $READS_PATH/Molestus_1_paired_trimmed_pairs_R1.fastq $READS_PATH/Molestus_1_paired_trimmed_pairs_R2.fastq Molestus_1
run_star $READS_PATH/Molestus_2_paired_trimmed_pairs_R1.fastq $READS_PATH/Molestus_2_paired_trimmed_pairs_R2.fastq Molestus_2
run_star $READS_PATH/Molestus_4_paired_trimmed_pairs_R1.fastq $READS_PATH/Molestus_4_paired_trimmed_pairs_R2.fastq Molestus_4
run_star $READS_PATH/Pipiens_1_paired_trimmed_pairs_R1.fastq $READS_PATH/Pipiens_1_paired_trimmed_pairs_R2.fastq Pipiens_1
run_star $READS_PATH/Pipiens_2_paired_trimmed_pairs_R1.fastq $READS_PATH/Pipiens_2_paired_trimmed_pairs_R2.fastq Pipiens_2
run_star $READS_PATH/Pipiens_3_paired_trimmed_pairs_R1.fastq $READS_PATH/Pipiens_3_paired_trimmed_pairs_R2.fastq Pipiens_3
```


## Remove duplicated

`/proj/snic2022-23-541/Rohan/Data/RNAseq/Aligned.sortedByCoord.out.bam`

```bash
#!/bin/bash
#SBATCH -A naiss2023-5-461
#SBATCH -p core -n 12
#SBATCH -t 2:00:00
#SBATCH -J MarkDuplicates
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se
#SBATCH --output=mark_duplicates_%j.out
#SBATCH --error=mark_duplicates_%j.err
#SBATCH --array=1-7


# Define input and output files

filename=$(sed -n "${SLURM_ARRAY_TASK_ID}p" List_md_2.txt)
input_file="/proj/snic2022-23-541/Rohan/Data/RNAseq/Indiv/${filename}.bam"
output_file="/proj/snic2022-23-541/Rohan/Analysis/MarkDuplicates/Run_2/${filename}.md.bam"
metrics_file="/proj/snic2022-23-541/Rohan/Analysis/MarkDuplicates/Run_2/${filename}.md.metrics.txt"


# Run Picard MarkDuplicates
picard MarkDuplicates \
    I="${input_file}" \
    O="${output_file}" \
    M="${metrics_file}" \
    REMOVE_DUPLICATES=true \
    CREATE_INDEX=true

echo "MarkDuplicates completed successfully."



```

### Create an GATK genome index

```bash
gatk CreateSequenceDictionary -R VectorBase-66_CquinquefasciatusJHB2020_Genome_headers.fasta
```
### Split'N'Trim and reassign mapping qualities

```bash
# command line in the git hub, but doesnt work 
java -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R ref.fasta -I dedupped.bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

```bash
#Test
# This worked!

gatk SplitNCigarReads -R ../../../Data/Genome/VectorBase-66_CquinquefasciatusJHB2020_Genome_headers.fasta -I Molestus_1Aligned.sortedByCoord.out.md.bam -O split.bam

```

### Script for SplitNcigars

```bash
#!/bin/bash
#SBATCH -A naiss2023-5-461
#SBATCH -p core -n 12
#SBATCH -t 2:00:00
#SBATCH -J MarkDuplicates
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se
#SBATCH --output=mark_duplicates_%j.out
#SBATCH --error=mark_duplicates_%j.err
#SBATCH --array=1-7


# Define input and output files

filename=$(sed -n "${SLURM_ARRAY_TASK_ID}p" List_md_2.txt)
input_file="/proj/snic2022-23-541/Rohan/Data/RNAseq/Indiv/${filename}.bam"
output_file="/proj/snic2022-23-541/Rohan/Analysis/MarkDuplicates/Run_2/${filename}.md.bam"
metrics_file="/proj/snic2022-23-541/Rohan/Analysis/MarkDuplicates/Run_2/${filename}.md.metrics.txt"


# Run Picard MarkDuplicates
picard MarkDuplicates \
    I="${input_file}" \
    O="${output_file}" \
    M="${metrics_file}" \
    REMOVE_DUPLICATES=true \
    CREATE_INDEX=true

echo "MarkDuplicates completed successfully."

```

### Variant calling with gatk and HaplotypeCaller

```bash
#!/bin/bash
#SBATCH -A naiss2023-5-461
#SBATCH -p core -n 12
#SBATCH -t 4-00:00:00
#SBATCH -J HaplotypeCaller
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se
#SBATCH --output=haplotypecaller_%j.out
#SBATCH --error=haplotypecaller_%j.err
#SBATCH --array=1-6

# Load necessary modules
module load bioinfo-tools
module load GATK/4.0.12.0  # Adjust version based on your environment

# Define input and output paths
input_dir="/proj/snic2022-23-541/Rohan/Analysis/MarkDuplicates/Run_2"
reference="/proj/snic2022-23-541/Rohan/Data/Genome/VectorBase-66_CquinquefasciatusJHB2020_Genome_headers.fasta"

# Read the file names from a list
filename=$(sed -n "${SLURM_ARRAY_TASK_ID}p" List_split_bams.txt)
input_file="${input_dir}/${filename}"
# Define output file name
output_file="${input_dir}/${filename%.bam}.vcf"

# Run GATK HaplotypeCaller
gatk HaplotypeCaller \
    -R "$reference" \
    -I "$input_file" \
    -O "$output_file"

echo "HaplotypeCaller completed successfully for ${filename}."

```












```
# Variant calling
https://github.com/broadgsa/gatk/blob/master/doc_archive/methods/Calling_variants_in_RNAseq.md

https://evomics.org/wp-content/uploads/2024/01/Variant-calling-Workshop-on-Genomics-2024-Cesky-Krumlov.pdf

https://learn.gencore.bio.nyu.edu/variant-calling/

https://qcb.ucla.edu/wp-content/uploads/sites/14/2020/03/VariantCallingWithGATK_WINTER2020.pdf

https://nbisweden.github.io/workshop-ngsintro/2209/slide_vc.pdf

https://nbisweden.github.io/workshop-ngsintro/2403/topics/vc/lab_vc.html


**VIDEO**

https://duckduckgo.com/?q=varuancalling&atb=v437-1&iax=videos&ia=videos&iai=https%3A%2F%2Fwww.youtube.com%2Fwatch%3Fv%3DmKqdfdtv0cI

```