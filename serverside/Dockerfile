FROM rocker/shiny-verse

RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y curl libcurl4-openssl-dev libssl-dev libxml2-dev vim samtools unar libbz2-dev liblzma-dev
RUN R -e "tryCatch(install.packages(\"BiocManager\"),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"optparse\", dependencies = TRUE ),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"reshape\", dependencies = TRUE ),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"shinydashboard\", dependencies = TRUE ),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"dashboardthemes\", dependencies = TRUE ),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"shinyjs\", dependencies = TRUE ),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"shinycssloaders\", dependencies = TRUE ),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"shinyWidgets\", dependencies = TRUE ),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"RMariaDB\", dependencies = TRUE ),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"shinyFiles\", dependencies = TRUE ),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"reshape2\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"ggpubr\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"plotly\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"VennDiagram\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"Rsamtools\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"Rsubread\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"DT\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"biomaRt\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"pheatmap\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"ComplexHeatmap\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"viridis\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"DESeq2\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"DEXSeq\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"GenomicRanges\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"rtracklayer\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"IRanges\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"IsoformSwitchAnalyzeR\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"genomation\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"ChIPseeker\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"ChIPpeakAnno\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"EnsDb.Hsapiens.v86\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"org.Hs.eg.db\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"GenomicFeatures\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"Gviz\",ask=FALSE, force=TRUE),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 && unar htslib-1.10.2.tar.bz2 && cd htslib-1.10.2/ &&  ./configure --prefix=/usr/ && make && make install && cd ..  && rm htslib*.bz2

RUN wget https://raw.githubusercontent.com/rocker-org/shiny/master/shiny-server.sh -O /usr/bin/shiny-server.sh; chmod +x /usr/bin/shiny-server.sh

COPY /localscratch/marisol/DASiRe/databases/encode_bedNarrowPeak_files/ /examples/databases/
COPY /localscratch/marisol/DASiRe/databases/mart_export.txt /databases/
COPY /localscratch/marisol/DASiRe/databases/splicing_factors_list.txt /examples/databases/
RUN chmod 777 /databases -R

EXPOSE 3838
CMD ["/usr/bin/shiny-server.sh"]


