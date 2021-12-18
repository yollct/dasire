#!/bin/bash

#Filepaths all begin with '/'' suggesting that these are absolute filepaths for processes within the docker.
# /MOUNT == DASiRe/MOUNT/ directory in your local machine.

#This file hosts the variables for DASiRe's preprocessing.
fasta=/MOUNT/input/GRCh38.primary_assembly.genome.fa
gtf=/MOUNT/input/gencode.v36.annotation.gtf
gff=/MOUNT/input/gencode.v36.annotation.gff
transcripts_fasta=/MOUNT/input/gencode.v36.transcripts.fa
exon_custom_annotation=/MOUNT/input/gencode.v36.annotation.dexseq.gff
metadata=/MOUNT/input/metadata.txt
adapters=/MOUNT/input/adapters_pe.fa
majiq-config=/MOUNT/input/majiq-config.txt
