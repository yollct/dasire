library(shinyjs)
library(shiny)
setwd("/nfs/home/students/chit/DASiRe/R")

readfile <- function(x, dir){
  f <- read.csv(paste0(dir,"/",x), sep="\t", header=T)
  return(f)
}
make_gene_matrix <- function(dir){
  genefiles <- list.files(dir)
  genefiles <- genefiles[grepl("raw_count_matrix", genefiles)]
  
  rr <- lapply(genefiles, readfile, "../examples/gene_counts") 
  r <- rr %>% reduce(inner_join, by=colnames(rr[1][[1]])[1])
  
  colnames(r) <- lapply(colnames(r), function(x){strsplit(x,"[.]")[[1]][1]})
  return(r)
}

