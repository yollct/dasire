#!/bin/bash

# Enabling bash strict mode: (http://redsymbol.net/articles/unofficial-bash-strict-mode/#expect-nonzero-exit-status)
set -euo pipefail

source /MOUNT/config.sh
for i in $(find input/ -name "*fastq.gz" -nowarn )
do echo "uncompressing $i"
gunzip $i
done
wait

echo "step 1 of 4: Building indeces for STAR and Kallisto"
# -[ ] check if indeces are already built
echo "Indexing reference genome with STAR"
if ! test -f /MOUNT/index/starindex/genomeParameters.txt
	then mkdir -p /MOUNT/index/starindex
		STAR \
		--runMode genomeGenerate \
		--genomeDir /MOUNT/index/starindex \
		--genomeFastaFiles $fasta \
		--runThreadN 4 \
		--sjdbGTFfile $gtf \
		-sjdbOverhang 100 \
		--outFileNamePrefix output/STAR/
 fi 
# -[ ] check if indeces are already built
echo "Indexing reference genome with Kallisto"
if ! test -f /MOUNT/index/kallisto-index
 then kallisto index -i /MOUNT/index/kallisto-index $transcripts_fasta
fi


for i in $(find input/ -name "*fastq" -nowarn  | sed 's/..fastq$//g' | sed 's/^input\///g')
do fastqcount=$(find input/ -name "${i}?.fastq" | wc -l)
echo Sample: input/${i}*  \| Number of Fastqfiles: $fastqcount

# -[ ]  check if fastqcount is 2| 0: check input name | case: greater: single  # currently assuming paired


echo "step 2 of 4: trimming reads for Sample: $i  & Aligning to the genome with STAR & Kallisto"
trimmomatic PE input/${i}1.fastq input/${i}2.fastq output/${i}r1_trimmed.fastq.gz output/${i}r1_unpaired.fastq.gz output/${i}r2_trimmed.fastq.gz output/${i}r2_unpaired.fastq.gz ILLUMINACLIP:$adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Genome Alignment with STAR for Sample: $i"
STAR --genomeDir /MOUNT/index/starindex --readFilesIn output/${i}r1_trimmed.fastq.gz output/${i}r2_trimmed.fastq.gz --outSAMattributes NH HI AS nM NM MD --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastxm --outFileNamePrefix output/$i

echo "Transcript quantification with Kallisto for Sample: $i"
kallisto quant -i /MOUNT/index/kallisto-index -o pseudocounts/$i --bias ${i}r1_trimmed.fastq ${i}r2_trimmed.fastq
done



mkdir output/exon_counts
echo "step 3 of 4: Quantifying exons with DEXSeq and AS Events with Majiq"
# for i in $(find output/ -name "*.bam" -nowarn  )
# do
# 	echo "indexing $i"
# 	samtools index $i $i.bai
# 	echo "counting exons $i"
#     python dexseq_count.py -s no -r pos -f bam $exon_custom_annotation $i output/exon_counts/$i.exon_counts.txt
#     #Majiq
# # 		# build sperate config file for each BAM file
# # 		majiq_basename=$(basename -s .bam $i)
# # 		outdir_name=$(basename -s .bam $i)_output
# # 		mkdir -p /MOUNT/output/MAJIQ/$outdir_name
# #
# # 		config=/MOUNT/output/MAJIQ/$outdir_name/config.txt
# # 		echo "[info]" > $config
# # 		echo "readlen=75" >> $config
# # 		echo "bamdirs=/MOUNT/output" >> $config
# # 		echo "genome=hg38" >> $config
# # 		echo "strandness=None" >> $config
# # 		echo "[experiments]" >> $config
# # 		echo "BAM=$majiq_basename" >> $config

# # 		echo "building MAJIQ reference ..."
# # 		majiq build $gff -c $config -j 4 -o /MOUNT/output/MAJIQ/$outdir_name/build
# # 		wait

# # 		#get all .majiq files which were created with build
# # 	        majiqlist=$(ls -1p /MOUNT/output/MAJIQ/$outdir_name/build/*.majiq | xargs echo)

# # 		majiq psi $majiqlist -j 4 -o /MOUNT/output/MAJIQ/$outdir_name/psi -n "BAM"

# # 		# create voila.tsv outputfiles
# # 		voila tsv /MOUNT/output/MAJIQ/$outdir_name/build/splicegraph.sql /MOUNT/output/MAJIQ/$outdir_name/psi/*.voila -f /MOUNT/output/MAJIQ/$outdir_name/voila.tsv
# done
# wait



echo "step 4 of 4: DESeq, DEXSeq and IsoformSwitchAnalyzeR in R"
#Rscript deseq2_dexseq_isoformswicthanalyzer.Rscript -b $bamdir -l paired -gtf $gtf -o output -m $metadata -gff $exon_custom_annotation -f $transcripts_fasta

echo "All finished visit https://exbio.wzw.tum.de/dasire/ for your next steps. Upload files from directory MOUNT/output to the webserver."