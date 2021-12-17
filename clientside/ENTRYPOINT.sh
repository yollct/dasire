#!/bin/bash

# Enabling bash strict mode: (http://redsymbol.net/articles/unofficial-bash-strict-mode/#expect-nonzero-exit-status)
set -euo pipefail

source /MOUNT/config.sh
for i in $(find input/ -name "*fastq.gz" -nowarn )
do echo "uncompressing $i"
gunzip $i
done
wait

echo "step 1 of 6: Building indeces for STAR and Kallisto"
# -[ ] check if indeces are already built
echo "Indexing reference genome with STAR"
mkdir -p /MOUNT/index/starindex
	    STAR \
		--runMode genomeGenerate \
		--genomeDir /MOUNT/index/starindex \
		--genomeFastaFiles $fasta \
		--runThreadN 4 \
		--sjdbGTFfile $gtf \
		-sjdbOverhang 100 \
		--outFileNamePrefix output/STAR/

# -[ ] check if indeces are already built
echo "Indexing reference genome with Kallisto"
kallisto index -i /MOUNT/index/kallisto-index $transcripts_fasta

for i in $(find input/ -name "*fastq" -nowarn  | sed 's/..fastq$//g')
do fastqcount=$(find input/ -name "${i}?.fastq" | wc -l)
echo Sample: input/${i}*  \| Number of Fastqfiles: $fastqcount

# -[ ]  check if fastqcount is 2| 0: check input name | case: greater: single  # currently assuming paired


echo "step 2 of 6: trimming reads for Sample: $i  & Aligning to the genome with STAR"
trimmomatic PE input/${i}1.fastq input/${i}2.fastq output/${i}r1_trimmed.fastq.gz output/${i}r1_unpaired.fastq.gz output/${i}r2_trimmed.fastq.gz output/${i}r2_unpaired.fastq.gz ILLUMINACLIP:$adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

STAR --genomeDir /MOUNT/index/starindex --readFilesIn output/${i}r1_trimmed.fastq.gz output/${i}r2_trimmed.fastq.gz --outSAMattributes NH HI AS nM NM MD --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastxm --outFileNamePrefix output/$i
done
wait

mkdir output/exon_counts
echo "step 3 of 6: Counting exons with DEXSeq and Majiq"
for i in $(find output/ -name "*.bam" -nowarn  )
do
	echo "indexing $i"
	samtools index $i $i.bai
	echo "counting exons $i"
    python dexseq_count.py -s no -r pos -f bam $exon_custom_annotation $i output/exon_counts/$i.exon_counts.txt
    #Majiq
# for filename in $(cat /tmp/controlbamlist)
# 	do
# 		# build sperate config file for each BAM file
# 		majiq_basename=$(basename -s .bam $filename)
# 		outdir_name=$(basename -s .bam $filename)_output
# 		mkdir -p $outdir/$outdir_name
# 		config=$outdir/$outdir_name/config.txt
# 		echo "[info]" > $config
# 		echo "readlen=$read_length" >> $config
# 		#bam-directory depending on type of run
# 		echo "bamdirs=$controlbam" >> $config
# 		echo "genome=hg38" >> $config
# 		echo "strandness=None" >> $config
# 		echo "[experiments]" >> $config
# 		echo "BAM=$majiq_basename" >> $config

# 		echo "building MAJIQ reference ..."
# 		majiq build $gff -c $config -j $ncores -o $outdir/$outdir_name/build
# 		wait

# 		#get all .majiq files which were created with build
# 	        majiqlist=$(ls -1p $outdir/$outdir_name/build/*.majiq | xargs echo)

# 		majiq psi $majiqlist -j $ncores -o $outdir/$outdir_name/psi -n "BAM"

# 		# create voila.tsv outputfiles
# 		voila tsv $outdir/$outdir_name/build/splicegraph.sql $outdir/$outdir_name/psi/*.voila -f $outdir/$outdir_name/voila.tsv
# 		wait
# done

done
wait

echo "step 4 of 6: Pseudoaligning with kallisto"
mkdir output/pseudocounts

for i in $(find input/ -name "output/${i}r?_trimmed.fastq" | sed 's/r._trimmed.fastq//g' )
do echo "Transcript quantification with Kallisto for Sample: $i"
kallisto quant -i /MOUNT/index/kallisto-index -o pseudocounts/$i --bias ${i}r1_trimmed.fastq ${i}r2_trimmed.fastq
done



echo "step 6 of 6: DESeq, DEXSeq and IsoformSwitchAnalyzeR in R"
#Rscript deseq2_dexseq_isoformswicthanalyzer.Rscript -b $bamdir -l paired -gtf $gtf -o output -m $metadata -gff $exon_custom_annotation -f $transcripts_fasta

echo "All finished visit https://exbio.wzw.tum.de/dasire/ for your next steps. Upload files from directory MOUNT/output to the webserver."