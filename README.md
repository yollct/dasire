# DASiRe

## What is DASiRe?

**Direct Alternative Splicing Regulator predictor (DASiRe)** is a web application that allows non-expert users to perform different types of splicing analysis from RNA-seq experiments and also incorporates ChIP-seq data of a DNA-binding protein of interest to evaluate whether its presence is associated with the splicing changes detected in the RNA-seq dataset. 

DASiRe is an accessible web-based platform that performs the analysis of raw RNA-seq and ChIP-seq data to study the relationship between DNA-binding proteins and alternative splicing regulation. It provides a fully integrated pipeline that takes raw reads from RNA-seq and performs extensive splicing analysis by incorporating the three current methodological approaches to study alternative splicing: isoform switching, exon and event-level. Once the initial splicing analysis is finished, DASiRe performs ChIP-seq peak enrichment in the spliced genes detected by each one of the three approaches. 

## Serverside
Start the serverside app by running the following command from the top directory of this git.

`docker run -p 3838:3838 -v $(pwd)/serverside/R:/srv/shiny-server/ -d --rm dasire-server:0.1`
