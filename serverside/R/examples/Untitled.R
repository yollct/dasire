library(tidyverse)

gsc <- read.csv("./sorted_pos_tpm.csv", sep="\t", )
gene <- read.csv("/Users/chit/Desktop/spycone/spycone/data/covid_hs/genelist.csv",sep="\t")
gsc$gene <- gene$symb

genecount <- aggregate(.~gene,data=gsc,FUN=sum) %>% filter(gene!="")
write.table(genecount, "./sorted_pos_tpm.csv", sep="\t")
