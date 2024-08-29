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
grep -e ">" VectorBase-66_CquinquefasciatusJHB2020_Genome.fasta > Headers_genome.txt
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
#SBATCH --output=hisat_%A_%a.out
#SBATCH --error=hisat_%A_%a.err


# The programs are installed in Rohan env

# Define input directories

input_dir="/proj/snic2022-23-541/Rohan/Data/RNAseq"
output_dir="/proj/snic2022-23-541/Rohan/Analysis/BWA"
genome_index="/proj/snic2022-23-541/Rohan/Data/Genome/Index/VectorBase-66_CquinquefasciatusJHB2020_Genome_headers.fasta"


# Define sample name
sample=$(sed -n "${SLURM_ARRAY_TASK_ID} p" List.txt)

# Define input files
r1="${input_dir}/${sample}_R1_001_val_1.fq.gz"
r2="${input_dir}/${sample}_R2_001_val_2.fq.gz"

# Define output file
output="${output_dir}/${sample}.sam"

# Run BWA for Illumina/454/IonTorrent paired-end reads longer than ~70bp:

bwa mem "${genome_index}" "${r1}" "${r2}" -t 8 > "${output}"

```

```bash
  bwa mem -t 8 "${genome_index}" "${r1}" -2 "${r2}" > "${output}"

```

### MarkDuplicates from Picard to remove PCR duplicates

```bash
#!/bin/bash
#SBATCH -A naiss2023-5-461
#SBATCH -p core -n 4
#SBATCH -t 2:00:00
#SBATCH -J MarkDuplicates
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se
#SBATCH --output=mark_duplicates_%j.out
#SBATCH --error=mark_duplicates_%j.err


# Define input and output files
input_bam="path/to/your/input.bam"
output_bam="path/to/your/output.dedup.bam"
metrics_file="path/to/your/metrics.txt"

# Run Picard MarkDuplicates
picard MarkDuplicates \
    I="${input_bam}" \
    O="${output_bam}" \
    M="${metrics_file}" \
    REMOVE_DUPLICATES=true \
    CREATE_INDEX=true

echo "MarkDuplicates completed successfully."



```

```
https://github.com/broadgsa/gatk/blob/master/doc_archive/methods/Calling_variants_in_RNAseq.md
```