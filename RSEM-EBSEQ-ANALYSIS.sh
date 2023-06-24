#!/bin/bash
#Comment -  to be submitted by: sbatch RSEM-EBSEQ-ANALYSIS.sh
#SBATCH --time=30:00:00		# 30 hours max
#SBATCH --nodes=1 		# 1 compute node
#SBATCH --ntasks-per-node=24 	# 24 cpu cores
#SBATCH --partition=batch	# requests a compute node
#SBATCH --job-name=RSEM-EBSEQ-ANALYSIS
#SBATCH --output=RSEM-EBSEQ-ANALYSIS.out

#Comment - batch job setup complete

#Comment – load a program module, for example Python
cd
source loadAnaconda.sh

#Comment – path to execution directory, for example
mkdir -p ~/rsem/
cd $HOME/rsem/

#Comment – Proceed with pipeline

mkdir -p ~/rsem/data
mkdir -p ~/rsem/ref

#Comment - download sequence data in FASTQ format (This pipeline fixes an issue with the other notebook, this will now handle paired data whereas the kallisto pipeline cant handle those properly and treated them as single ended reads!)
cd ~/rsem/data
export PATH="~/software/sratoolkit.3.0.0-centos_linux64/bin:$PATH"
prefetch SRR24062652
prefetch SRR24062653
prefetch SRR24062654
prefetch SRR24062655
prefetch SRR24062656
prefetch SRR24062657

fastq-dump --split-files SRR24062652/
fastq-dump --split-files SRR24062653/
fastq-dump --split-files SRR24062654/
fastq-dump --split-files SRR24062655/
fastq-dump --split-files SRR24062656/
fastq-dump --split-files SRR24062657/

rm -r SRR24062652
rm -r SRR24062653
rm -r SRR24062654
rm -r SRR24062655
rm -r SRR24062656
rm -r SRR24062657

#Comment - create mouse genome reference for RSEM
cd ~/rsem/ref

wget ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-82/gtf/mus_musculus/Mus_musculus.GRCm38.82.chr.gtf.gz
gunzip Mus_musculus.GRCm38.dna.toplevel.fa.gz
gunzip Mus_musculus.GRCm38.82.chr.gtf.gz
rsem-prepare-reference --gtf Mus_musculus.GRCm38.82.chr.gtf --bowtie2 Mus_musculus.GRCm38.dna.toplevel.fa mouse_ref

#Comment - Clean sequences with Trimmomatic
cd ~/rsem/data
java -jar ~/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
  -threads 24 \
  SRR24062652_1.fastq SRR24062652_2.fastq \
  SRR24062652_1.trimmed.fastq SRR24062652_1.trimmedOrphan.fastq \
  SRR24062652_2.trimmed.fastq SRR24062652_2.trimmedOrphan.fastq \
  SLIDINGWINDOW:4:20 MINLEN:75


java -jar ~/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
  -threads 24 \
  SRR24062653_1.fastq SRR24062653_2.fastq \
  SRR24062653_1.trimmed.fastq SRR24062653_1.trimmedOrphan.fastq \
  SRR24062653_2.trimmed.fastq SRR24062653_2.trimmedOrphan.fastq \
  SLIDINGWINDOW:4:20 MINLEN:75

java -jar ~/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
  -threads 24 \
  SRR24062654_1.fastq SRR24062654_2.fastq \
  SRR24062654_1.trimmed.fastq SRR24062654_1.trimmedOrphan.fastq \
  SRR24062654_2.trimmed.fastq SRR24062654_2.trimmedOrphan.fastq \
  SLIDINGWINDOW:4:20 MINLEN:75

java -jar ~/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
  -threads 24 \
  SRR24062655_1.fastq SRR24062655_2.fastq \
  SRR24062655_1.trimmed.fastq SRR24062655_1.trimmedOrphan.fastq \
  SRR24062655_2.trimmed.fastq SRR24062655_2.trimmedOrphan.fastq \
  SLIDINGWINDOW:4:20 MINLEN:75

java -jar ~/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
  -threads 24 \
  SRR24062656_1.fastq SRR24062656_2.fastq \
  SRR24062656_1.trimmed.fastq SRR24062656_1.trimmedOrphan.fastq \
  SRR24062656_2.trimmed.fastq SRR24062656_2.trimmedOrphan.fastq \
  SLIDINGWINDOW:4:20 MINLEN:75

java -jar ~/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
  -threads 24 \
  SRR24062657_1.fastq SRR24062657_2.fastq \
  SRR24062657_1.trimmed.fastq SRR24062657_1.trimmedOrphan.fastq \
  SRR24062657_2.trimmed.fastq SRR24062657_2.trimmedOrphan.fastq \
  SLIDINGWINDOW:4:20 MINLEN:75

rm SRR24062652_1.fastq
rm SRR24062652_2.fastq
rm SRR24062653_1.fastq
rm SRR24062653_2.fastq
rm SRR24062654_1.fastq
rm SRR24062654_2.fastq
rm SRR24062655_1.fastq
rm SRR24062655_2.fastq
rm SRR24062656_1.fastq
rm SRR24062656_2.fastq
rm SRR24062657_1.fastq
rm SRR24062657_2.fastq

#Comment - Calculate expression counts with RSEM
mkdir -p ~/rsem/exp
cd ~/rsem

rsem-calculate-expression -p 24 \
  --paired-end --bowtie2 --estimate-rspd --append-names --output-genome-bam \
  data/SRR24062652_1.trimmed.fastq \
  data/SRR24062652_2.trimmed.fastq \
  ref/mouse_ref exp/KOWT_52
rsem-calculate-expression -p 24 \
  --paired-end --bowtie2 --estimate-rspd --append-names --output-genome-bam \
  data/SRR24062653_1.trimmed.fastq \
  data/SRR24062653_2.trimmed.fastq \
  ref/mouse_ref exp/KOWT_53

rsem-calculate-expression -p 24 \
  --paired-end --bowtie2 --estimate-rspd --append-names --output-genome-bam \
  data/SRR24062654_1.trimmed.fastq \
  data/SRR24062654_2.trimmed.fastq \
  ref/mouse_ref exp/KOWT_54

rsem-calculate-expression -p 24 \
  --paired-end --bowtie2 --estimate-rspd --append-names --output-genome-bam \
  data/SRR24062655_1.trimmed.fastq \
  data/SRR24062655_2.trimmed.fastq \
  ref/mouse_ref exp/KOWT_55

rsem-calculate-expression -p 24 \
  --paired-end --bowtie2 --estimate-rspd --append-names --output-genome-bam \
  data/SRR24062656_1.trimmed.fastq \
  data/SRR24062656_2.trimmed.fastq \
  ref/mouse_ref exp/KOWT_56

rsem-calculate-expression -p 24 \
  --paired-end --bowtie2 --estimate-rspd --append-names --output-genome-bam \
  data/SRR24062657_1.trimmed.fastq \
  data/SRR24062657_2.trimmed.fastq \
  ref/mouse_ref exp/KOWT_57

#Comment - Plot expression results
cd ~/rsem/exp
rsem-plot-model KOWT_52 KOWT_52.pdf

#Comment - Differential expression

rsem-generate-data-matrix  \
  KOWT_52.genes.results \
  KOWT_53.genes.results \
  KOWT_54.genes.results \
  KOWT_55.genes.results \
  KOWT_56.genes.results \
  KOWT_57.genes.results \
  > geneMat.txt
rsem-run-ebseq geneMat.txt 3,3 geneMat.results
rsem-control-fdr geneMat.results 0.05 geneMat.de.txt

rsem-generate-ngvector ../ref/mouse_ref.transcripts.fa mouse_ref
rsem-generate-data-matrix  \
  KOWT_52.isoforms.results \
  KOWT_53.isoforms.results \
  KOWT_54.isoforms.results \
  KOWT_55.isoforms.results \
  KOWT_56.isoforms.results \
  KOWT_57.isoforms.results \
  > isoMat.txt
rsem-run-ebseq --ngvector mouse_ref.ngvec isoMat.txt 3,3 isoMat.results
rsem-control-fdr isoMat.results 0.05 isoMat.de.txt

mv ~/get6MostDE_genes.R ~/rsem
cd ~/rsem
Rscript get6MostDE_genes.R