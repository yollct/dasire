FROM rocker/shiny-verse

RUN R -e "install.packages(c(\"BiocManager\",\"optparse\",\"reshape\"))"

RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y curl libcurl4-openssl-dev libssl-dev libxml2-dev vim samtools unar libbz2-dev liblzma-dev

RUN wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 && unar htslib-1.10.2.tar.bz2 && cd htslib-1.10.2/ &&  ./configure --prefix=/usr/ && make && make install && cd ..  && rm htslib*.bz2

RUN R -e "BiocManager::install(\"Rsamtools\");  BiocManager::install(\"Rsubread\")" 
RUN R -e "BiocManager::install(\"DEXSeq\"); BiocManager::install(\"DESeq2\"); BiocManager::install(\"IsoformSwitchAnalyzeR\")"

EXPOSE 3838
CMD ["/usr/bin/shiny-server.sh"]
