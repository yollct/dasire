#!/bin/bash

# Enabling bash strict mode: (Read more about it at: http://redsymbol.net/articles/unofficial-bash-strict-mode/#expect-nonzero-exit-status)
set -euo pipefail

# Importing configuration
source /MOUNT/config.sh

# uncompressing all fastq.gz files
for i in $(find input/ -name "*fastq.gz" -nowarn )
do echo "uncompressing $i"
pigz -d $i
done
wait

# Building Indeces [ Duration: 1.5 hours (on Human genome)]
# ---------------------------------------------------------------------------------------------------------
echo "DASiRe Step 1 of 4: Building indeces for STAR and Kallisto"

# Building Kallisto's Index
if ! test -f /MOUNT/index/kallisto-index
 then echo "Indexing reference genome with Kallisto"
 	kallisto index -i /MOUNT/index/kallisto-index $transcripts_fasta
 else echo "Index found at index/kallisto-index, reusing it:"
 	ls /MOUNT/index/kallisto-index -lah
fi

# Building STAR's Index
if ! test -f /MOUNT/index/starindex/genomeParameters.txt
	then mkdir -p /MOUNT/index/starindex || true # to allow mkdir to fail gracefully, we add "|| true"
		echo "Indexing reference genome with STAR"
		STAR \
		--runMode genomeGenerate \
		--genomeDir /MOUNT/index/starindex \
		--genomeFastaFiles $fasta \
		--runThreadN 4 \
		--sjdbGTFfile $gtf \
		-sjdbOverhang 100 \
		--outFileNamePrefix output/STAR/
	else echo "Indexes found at index/starindex, reusing them:"
	ls /MOUNT/index/starindex/ -lah
 fi 

# Alignments
# ---------------------------------------------------------------------------------------------------------

# find input/ -name "*fastq" -nowarn  | sed 's/..fastq$//g' | sed 's/^input\///g'| sort | uniq | parallel -j 4 "
# fastqcount=$(find input/ -name "{}?.fastq" | wc -l)
# 	echo "Sample: input/{}*  \| Number of Fastqfiles: $fastqcount"
# 	sample_basename=$(echo $i| sed 's/_$//g')
# 	# -[ ]  check if fastqcount is 2| 0: check input name | case: greater: single  # currently assuming paired

# 	echo "DASiRe Step 2 of 4: trimming reads for Sample: $sample_basename  & Aligning to the genome with STAR & Kallisto"
# 	if ! test -f output/{}r1_trimmed.fastq
# 		then trimmomatic PE input/{}1.fastq input/{}2.fastq output/{}r1_trimmed.fastq output/{}r1_unpaired.fastq output/{}r2_trimmed.fastq output/{}r2_unpaired.fastq ILLUMINACLIP:$adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
# 	fi

# 	echo "Transcript quantification with Kallisto for Sample: $sample_basename"
# 	mkdir -p output/pseudocounts/$sample_basename || true # to allow mkdir to fail gracefully, we add "|| true"
# 	if ! test -f output/pseudocounts/$sample_basename/abundance.tsv 
# 		then kallisto quant -i /MOUNT/index/kallisto-index -o output/pseudocounts/$sample_basename --bias output/{}r1_trimmed.fastq output/{}r2_trimmed.fastq
# 	fi

# 	echo "Genome Alignment with STAR for Sample: $sample_basename" ; mkdir -p output/STAR/$sample_basename || true # to allow mkdir to fail gracefully, we add "|| true"
# 	if ! test -f output/STAR/${sample_basename}/Aligned.sortedByCoord.out.bam
# 		then STAR --genomeDir /MOUNT/index/starindex --readFilesIn output/{}r1_trimmed.fastq output/{}r2_trimmed.fastq --outSAMattributes NH HI AS nM NM MD --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastxm --outFileNamePrefix output/STAR/$sample_basename/
# 	fi
# "

for i in $(find input/ -name "*fastq" -nowarn  | sed 's/..fastq$//g' | sed 's/^input\///g'| sort | uniq)
	do fastqcount=$(find input/ -name "${i}?.fastq" | wc -l)
	echo "Sample: input/${i}*  \| Number of Fastqfiles: $fastqcount"
	sample_basename=$(echo $i| sed 's/_$//g')
	# -[ ]  check if fastqcount is 2| 0: check input name | case: greater: single  # currently assuming paired

	echo "DASiRe Step 2 of 4: trimming reads for Sample: $sample_basename  & Aligning to the genome with STAR & Kallisto"
	if ! test -f output/${i}r1_trimmed.fastq
		then trimmomatic PE input/${i}1.fastq input/${i}2.fastq output/${i}r1_trimmed.fastq output/${i}r1_unpaired.fastq output/${i}r2_trimmed.fastq output/${i}r2_unpaired.fastq ILLUMINACLIP:$adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	fi

	echo "Transcript quantification with Kallisto for Sample: $sample_basename"
	mkdir -p output/pseudocounts/$sample_basename || true # to allow mkdir to fail gracefully, we add "|| true"
	if ! test -f output/pseudocounts/$sample_basename/abundance.tsv 
		then kallisto quant -i /MOUNT/index/kallisto-index -o output/pseudocounts/$sample_basename --bias output/${i}r1_trimmed.fastq output/${i}r2_trimmed.fastq
	fi

	echo "Genome Alignment with STAR for Sample: $sample_basename" ; mkdir -p output/STAR/$sample_basename || true # to allow mkdir to fail gracefully, we add "|| true"
	if ! test -f output/STAR/${sample_basename}/Aligned.sortedByCoord.out.bam
		then STAR --genomeDir /MOUNT/index/starindex --readFilesIn output/${i}r1_trimmed.fastq output/${i}r2_trimmed.fastq --outSAMattributes NH HI AS nM NM MD --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastxm --outFileNamePrefix output/STAR/$sample_basename/
	fi
done

# Exp independent analysis
# ---------------------------------------------------------------------------------------------------------
mkdir -p output/exon_counts || true # to allow mkdir to fail gracefully, we add "|| true"
echo "DASiRe Step 3 of 4: Quantifying exons with DEXSeq"

find output/STAR/ -name "*.bam" -nowarn | xargs dirname {} | grep output/|  sed 's/output\/STAR\///g' | parallel -j 4 \
	"if ! test -f output/STAR/{}/Aligned.sortedByCoord.out.bam.bai
	then echo "indexing {}"
		samtools index output/STAR/{}/Aligned.sortedByCoord.out.bam output/STAR/{}/Aligned.sortedByCoord.out.bam.bai
	fi
	if ! test -f output/exon_counts/{}_exon_counts.txt;
	then echo "counting exons {}"
    	python /opt/conda/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -s no -r pos -f bam $exon_custom_annotation output/STAR/{}/Aligned.sortedByCoord.out.bam output/exon_counts/{}_exon_counts.txt
	fi"
wait


# Analysis with Experiment Design.
# ---------------------------------------------------------------------------------------------------------
echo "DASiRe Step 4 of 4: DESeq, DEXSeq and IsoformSwitchAnalyzeR in R"

Rscript deseq2_dexseq_isoformswitchanalyzer.Rscript -b $bamdir -l paired -gtf $gtf -o output -m $metadata -gff $exon_custom_annotation -f $transcripts_fasta

#Majiq
# # 		# build sperate config file for each BAM file
# # 		majiq_basename=$(basename -s .bam $i)
		outdir_name=$(basename -s .bam $i)
		mkdir -p /MOUNT/output/MAJIQ/$outdir_name|| true
# #
# # 		MajiqConfig=/MOUNT/output/MAJIQ/$outdir_name/config.txt
# # 		echo "[info]" > $MajiqConfig
# # 		echo "readlen=75" >> $MajiqConfig
# # 		echo "bamdirs=/MOUNT/output" >> $MajiqConfig
# # 		echo "genome=hg38" >> $MajiqConfig
# # 		echo "strandness=None" >> $MajiqConfig
# # 		echo "[experiments]" >> $MajiqConfig
# # 		echo "BAM=$majiq_basename" >> $MajiqConfig

		echo "building MAJIQ reference ..."
		majiq build $gff -c $MajiqConfig -j 4 -o /MOUNT/output/MAJIQ/$outdir_name/build

		#get all .majiq files which were created with build
	        majiqlist=$(ls -1p /MOUNT/output/MAJIQ/$outdir_name/build/*.majiq | xargs echo)
		majiq psi $majiqlist -j 4 -o /MOUNT/output/MAJIQ/$outdir_name/psi -n "BAM"
		# create voila.tsv outputfiles
		voila tsv /MOUNT/output/MAJIQ/$outdir_name/build/splicegraph.sql /MOUNT/output/MAJIQ/$outdir_name/psi/*.voila -f /MOUNT/output/MAJIQ/$outdir_name/voila.tsv

echo "All finished visit https://exbio.wzw.tum.de/dasire/ for your next steps. Upload files from directory MOUNT/output to the webserver."