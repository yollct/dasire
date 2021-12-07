#!bin/bash/
# I run this script like this: nohup bash preprocessing_single_end.sh [options] > preprocessing.out
# The preprocessing.out file is used to get the alignment statistics of kallisto.
#
# Requirements:
#   - STAR version 020201
#   - samtools version 1.7 (using htslib 1.7-2)
#   - dexseq_count.py: this script is included in DEXseq package and can be obtained as described here:
#       https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#23_HTSeq
#   - kallisto version 0.46.1
#   - Trimmomatic-0.39
#   - MultiQC version 1.10.dev0
#   - Genome annotation: Gencode Homo sapiens v.36
#
# This script works asuming it is run form the directory that contains: the gzipped raw fastq files used as input (i.e. input_SE.fastq.gz)
#
# Options:
#   --exon_custom_annotation    It is a custom annotation file with exons only, it is created from the original 
#                               annotation for that organism as input, in my case I used the gencode annotation for human (*gencode.v36.annotation.gtf* 
#                               to generate the custom exon annotation named *gencode.v36.annotation.dexseq.gff*). This file has to be created only once 
#                               for each different annotation. The DEXseq package has the python script to create this custom file. More info here:
#                               https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#24_Preparing_the_annotationAnnotation
#
#   --STAR_genome_index         The genome index folder created with STAR, requiered for alignment to genome with STAR
# 
#   --transcripts_index         The transcripts index file created with kallisto
#
#   --fragment_length           Fragment length of the sequencing library   

if [ $# -eq 0 ]
  then
    echo "No arguments supplied. See --help"
    exit 1
fi

usage() { echo "Usage: $0 [--fragment_length <fragment_length> --exon_custom_annotation <exon_custom_annotation> --STAR_genome_index <STAR_genome_index> --transcripts_index <transcripts_index>]"; exit 1; }

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --fragment_length) fragment_length="$2"; shift ;;
        --exon_custom_annotation) exon_custom_annotation="$2"; shift ;;
        --STAR_genome_index) STAR_genome_index="$2"; shift ;;
        --transcripts_index) transcripts_index="$2"; shift ;;
        --help) usage; exit 1 ;;
        *) echo "Unknown parameter passed: $1"; usage; exit 1 ;;
    esac
    shift
done

fragment_length=${fragment_length:-.}
exon_custom_annotation=${exon_custom_annotation:-.}
STAR_genome_index=${STAR_genome_index:-.}
transcripts_index=${transcripts_index:-.}

echo "step 1 of 5: trimming of reads with trimmomatic"
for i in *.fastq.gz
do
echo "trimming $i"
java -jar trimmomatic-0.39.jar SE input_SE.fastq.gz input_SE_trimmed.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

echo "uncompressing"
for i in *gz
do
echo "uncompressing $i"
gunzip $i
done

echo "step 2 of 5: alignment to the genome with star"
for i in *_trimmed.fastq
do 
echo "initiating alignment for $i"
STAR --genomeDir $STAR_genome_index --readFilesIn $i --outSAMattributes NH HI AS nM NM MD --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastxm --outFileNamePrefix $i.
done

echo "step 3 of 5: generate bai indexes"
for i in *.bam
do
	echo "initializing $i"
	samtools index $i $i.bai
	echo "finished $i"
done

echo "step 4 of 5: counts for ech exon"

mkdir exon_counts
for i in *.bam
do
echo "counting $i"
	python dexseq_count.py -s no -r pos -f bam $exon_custom_annotation $i exon_counts/$i.exon_counts.txt
done


echo "step 5 of 5: creates the pseudocounts with kallisto"
mkdir pseudocounts
for i in *_trimmed.fastq
do 
echo "counting $i"
kallisto quant -i $transcripts_index -o pseudocounts/$i --bias --single -l $fragment_length -s 20 $i
done

multiqc .

echo "all finished"
