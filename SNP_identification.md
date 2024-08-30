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
# Install GATK4
conda install bioconda::gatk4
# Install samtools
conda install bioconda::samtools
# Install VCF tools
conda install bioconda::vcftools

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
### Calculate the read coverage of positions in the genome

Do the first pass on variant calling by counting read coverage with bcftools. We will use the command mpileup. The flag -O b tells bcftools to generate a bcf format output file, -o specifies where to write the output file, and -f flags the path to the reference genome:
from `https://training.galaxyproject.org/training-material/topics/data-science/tutorials/bash-variant-calling/tutorial.html`

Example:

```bash
$ bcftools mpileup -O b -o "${output_file}" -f data/ref_genome/ecoli_rel606.fasta results/bam/SRR2584866.aligned.sorted.bam
````

```bash
#!/bin/bash
#SBATCH -A naiss2023-5-461
#SBATCH -p core -n 12
#SBATCH -t 2:00:00
#SBATCH -J bcftools_mpileup
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se
#SBATCH --output=mark_duplicates_%j.out
#SBATCH --error=mark_duplicates_%j.err
#SBATCH --array=1-7


filename=$(sed -n "${SLURM_ARRAY_TASK_ID}p" List_md.txt)
input_file="/proj/snic2022-23-541/Rohan/Analysis/MarkDuplicates/${filename}.md.bam"
output_file="/proj/snic2022-23-541/Rohan/Analysis/BCFtools/${filename}.raw.bcf"

bcftools mpileup -O b -o "${output_file}" -f /proj/snic2022-23-541/Rohan/Data/Genome/Index/VectorBase-66_CquinquefasciatusJHB2020_Genome_headers.fasta "${input_file}"

````



### Generate g.vcf files

```bash
gatk --java-options -Xmx7g HaplotypeCaller \
-R human_g1k_v37_chr2.fasta \
-ERC GVCF \
-I sample.bam \
-O sample.g.vcf
```







```
# Variant calling
https://github.com/broadgsa/gatk/blob/master/doc_archive/methods/Calling_variants_in_RNAseq.md

https://evomics.org/wp-content/uploads/2024/01/Variant-calling-Workshop-on-Genomics-2024-Cesky-Krumlov.pdf

https://learn.gencore.bio.nyu.edu/variant-calling/

https://qcb.ucla.edu/wp-content/uploads/sites/14/2020/03/VariantCallingWithGATK_WINTER2020.pdf

https://nbisweden.github.io/workshop-ngsintro/2209/slide_vc.pdf

https://nbisweden.github.io/workshop-ngsintro/2403/topics/vc/lab_vc.html
```