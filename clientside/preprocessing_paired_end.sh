#!bin/bash/
#
# Requirements:
#   - STAR version 020201
#   - samtools version 1.7 (using htslib 1.7-2)
#   - dexseq_count.py: this script is included in DEXseq package and can be obtained as described here:
#       https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#23_HTSeq
#   - kallisto version 0.46.1
#   - Trimmomatic-0.39
#   - Genome annotation: Gencode Homo sapiens v.36
#
# This script works asuming it is run form the a directory that contains: 
#   - the gzipped raw fastq files used as input: input_SE.fastq.gz
#   - It also works assuming 6 samples (12 files), each sample having the forward (r1) and reverse (r2) files.
#
# Options:
#   --prefix1,prefix2...        It is the name of the fastq file that includes both forward and reverse reads, example: for a sequencing file with forward
#                               and reverse files named like this: sample_r1.fastq.gz and sample_r2.fastq.gz, the prefix would be "sample".
#   --exon_custom_annotation    It is a custom annotation file with exons only, it is created from the original 
#                               annotation for that organism as input, in my case I used the gencode annotation for human (*gencode.v36.annotation.gtf* 
#                               to generate the custom exon annotation named *gencode.v36.annotation.dexseq.gff*). This file has to be created only one 
#                               for each different annotation. The DEXseq package has the python script to create this custom file. More info here:
#                               https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#24_Preparing_the_annotationAnnotation
#
#   --STAR_genome_index         The genome index folder created with STAR, requiered for alignment to genome with STAR
# 
#   --transcripts_index         The transcripts index file created with kallisto

if [ $# -eq 0 ]
  then
    echo "No arguments supplied. See --help"
    exit 1
fi

usage() { echo "Usage for paired end: $0 [--read_length <read_length>] [--prefix1 <prefix1>] [--prefix2 <prefix2>] [--prefix3 <prefix3>] [--prefix4 <prefix4>] [--prefix5 <prefix5>] [--prefix6 <prefix6>]"; exit 1; }

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --read_length) read_length="$2"; shift ;;
        --prefix1) prefix1="$2"; shift ;;
        --prefix2) prefix2="$2"; shift ;;
        --prefix3) prefix3="$2"; shift ;;
        --prefix4) prefix4="$2"; shift ;;
        --prefix5) prefix5="$2"; shift ;;
        --prefix6) prefix6="$2"; shift ;;
        --help) usage; exit 1 ;;
        *) echo "Unknown parameter passed: $1"; usage; exit 1 ;;
    esac
    shift
done

read_length=${read_length:-.}
prefix1=${prefix1:-.}
prefix2=${prefix2:-.}
prefix3=${prefix3:-.}
prefix4=${prefix4:-.}
prefix5=${prefix5:-.}
prefix6=${prefix6:-.}

echo "uncompressing"
for i in *.fastq.gz
do
echo "uncompressing $i"
gunzip $i
done

echo "step 1 of 5: trimming of reads"
read1=$prefix1\_r1.fastq.gz
read2=$prefix1\_r2.fastq.gz
echo "$read1"
echo "$read2" 
echo "initiating trimming of $prefix1"
java -jar trimmomatic-0.39.jar PE $read1 $read2 $prefix1\_r1_trimmed.fastq.gz $prefix1\_r1_unpaired.fastq.gz $prefix1\_r2_trimmed.fastq.gz $prefix1\_r2_unpaired.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

read1=$prefix2\_r1.fastq.gz
read2=$prefix2\_r2.fastq.gz
echo "$read1"
echo "$read2" 
echo "initiating trimming of $prefix2"
java -jar trimmomatic-0.39.jar PE $read1 $read2 $prefix2\_r1_trimmed.fastq.gz $prefix2\_r1_unpaired.fastq.gz $prefix2\_r2_trimmed.fastq.gz $prefix2\_r2_unpaired.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

read1=$prefix3\_r1.fastq.gz
read2=$prefix3\_r2.fastq.gz
echo "$read1"
echo "$read2" 
echo "initiating trimming of $prefix3"
java -jar trimmomatic-0.39.jar PE $read1 $read2 $prefix3\_r1_trimmed.fastq.gz $prefix3\_r1_unpaired.fastq.gz $prefix3\_r2_trimmed.fastq.gz $prefix3\_r2_unpaired.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

read1=$prefix4\_r1.fastq.gz
read2=$prefix4\_r2.fastq.gz
echo "$read1"
echo "$read2" 
echo "initiating trimming of $prefix4"
java -jar trimmomatic-0.39.jar PE $read1 $read2 $prefix4\_r1_trimmed.fastq.gz $prefix4\_r1_unpaired.fastq.gz $prefix4\_r2_trimmed.fastq.gz $prefix4\_r2_unpaired.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

read1=$prefix5\_r1.fastq.gz
read2=$prefix5\_r2.fastq.gz
echo "$read1"
echo "$read2" 
echo "initiating trimming of $prefix5"
java -jar trimmomatic-0.39.jar PE $read1 $read2 $prefix5\_r1_trimmed.fastq.gz $prefix5\_r1_unpaired.fastq.gz $prefix5\_r2_trimmed.fastq.gz $prefix5\_r2_unpaired.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

read1=$prefix6\_r1.fastq.gz
read2=$prefix6\_r2.fastq.gz
echo "$read1"
echo "$read2" 
echo "initiating trimming of $prefix6"
java -jar trimmomatic-0.39.jar PE $read1 $read2 $prefix6\_r1_trimmed.fastq.gz $prefix6\_r1_unpaired.fastq.gz $prefix6\_r2_trimmed.fastq.gz $prefix6\_r2_unpaired.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "step 2 of 5: alignment to the genome with star"
read1=$prefix1\_r1_trimmed.fastq
read2=$prefix1\_r2_trimmed.fastq
echo "$read1"
echo "$read2" 
echo "initiating alignment fo $prefix1"
STAR --genomeDir $STAR_genome_index --readFilesIn $read1 $read2 --outSAMattributes NH HI AS nM NM MD --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastxm --outFileNamePrefix $prefix1.

read1=$prefix2\_r1_trimmed.fastq
read2=$prefix2\_r2_trimmed.fastq
echo "$read1"
echo "$read2" 
echo "initiating alignment fo $prefix2"
STAR --genomeDir $STAR_genome_index --readFilesIn $read1 $read2 --outSAMattributes NH HI AS nM NM MD --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastxm --outFileNamePrefix $prefix2.

read1=$prefix3\_r1_trimmed.fastq
read2=$prefix3\_r2_trimmed.fastq
echo "$read1"
echo "$read2"  
echo "initiating alignment fo $prefix3"
STAR --genomeDir $STAR_genome_index --readFilesIn $read1 $read2 --outSAMattributes NH HI AS nM NM MD --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastxm --outFileNamePrefix $prefix3.

read1=$prefix4\_r1_trimmed.fastq
read2=$prefix4\_r2_trimmed.fastq
echo "$read1"
echo "$read2"  
echo "initiating alignment fo $prefix4"
STAR --genomeDir $STAR_genome_index --readFilesIn $read1 $read2 --outSAMattributes NH HI AS nM NM MD --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastxm --outFileNamePrefix $prefix4.

read1=$prefix5\_r1_trimmed.fastq
read2=$prefix5\_r2_trimmed.fastq
echo "$read1"
echo "$read2"  
echo "initiating alignment fo $prefix5"
STAR --genomeDir $STAR_genome_index --readFilesIn $read1 $read2 --outSAMattributes NH HI AS nM NM MD --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastxm --outFileNamePrefix $prefix5.

read1=$prefix6\_r1_trimmed.fastq
read2=$prefix6\_r2_trimmed.fastq
echo "$read1"
echo "$read2"  
echo "initiating alignment fo $prefix6"
STAR --genomeDir $STAR_genome_index --readFilesIn $read1 $read2 --outSAMattributes NH HI AS nM NM MD --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastxm --outFileNamePrefix $prefix6.


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

echo "counting $prefix1"
read1=$prefix1\_r1_trimmed.fastq
read2=$prefix1\_r2_trimmed.fastq
echo "$read1"
echo "$read2" 
kallisto quant -i $transcripts_index -o pseudocounts/$prefix1 --bias $read1 $read2

echo "counting $prefix2"
read1=$prefix2\_r1_trimmed.fastq
read2=$prefix2\_r2_trimmed.fastq
echo "$read1"
echo "$read2" 
kallisto quant -i $transcripts_index -o pseudocounts/$prefix2 --bias $read1 $read2

echo "counting $prefix3"
read1=$prefix3\_r1_trimmed.fastq
read2=$prefix3\_r2_trimmed.fastq
echo "$read1"
echo "$read2"  
kallisto quant -i $transcripts_index -o pseudocounts/$prefix3 --bias $read1 $read2

echo "counting $prefix4"
read1=$prefix4\_r1_trimmed.fastq
read2=$prefix4\_r2_trimmed.fastq
echo "$read1"
echo "$read2"  
kallisto quant -i $transcripts_index -o pseudocounts/$prefix4 --bias $read1 $read2

echo "counting $prefix5"
read1=$prefix5\_r1_trimmed.fastq
read2=$prefix5\_r2_trimmed.fastq
echo "$read1"
echo "$read2"  
kallisto quant -i $transcripts_index -o pseudocounts/$prefix5 --bias $read1 $read2

echo "counting $prefix6"
read1=$prefix6\_r1_trimmed.fastq
read2=$prefix6\_r2_trimmed.fastq
echo "$read1"
echo "$read2"  
kallisto quant -i $transcripts_index -o pseudocounts/$prefix6 --bias $read1 $read2

echo "all finished"
