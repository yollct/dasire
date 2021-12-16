library(shinyjs)
library(shiny)
library(plyr)
#setwd("/nfs/home/students/chit/DASiRe/R")

readfile <- function(x, dir){
  f <- read.csv(paste0(dir,"/",x), sep="\t", header=T)
  return(f)
}

make_gene_matrix <- function(dir){
  genefiles <- list.files(dir)
  genefiles <- genefiles[grepl("raw_count_matrix", genefiles)]
  
  gene_counts.list<-list()
  for (sample in genefiles){
    gene_counts.list[[sample]]<-read.csv(paste0(dir,"/",sample), sep="\t", header=T)
  }
  
  count_matrix <- join_all(gene_counts.list, by='ensembl_gene_id_version', type='left')
  colnames(count_matrix) <- lapply(colnames(count_matrix), function(x){strsplit(x,"[.]")[[1]][1]})
  rownames(count_matrix) <- count_matrix$ensembl_gene_id_version
  return(count_matrix)
}


# we want to use CHIP experiments of known splicing factors from encode to show as references in the plots.

# Import splicing factors list 
sf_list <- read.delim(file = "examples/splicing_factors_list.txt",header = T,sep = "\t",fill = TRUE)

# import reference ENCODE datasets 
ENCODE_metadata <- read.delim("examples/encode_bedNarrowPeak_files/metadata.tsv",sep="\t",header = TRUE)




