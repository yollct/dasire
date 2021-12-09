library(shiny)
library(shinydashboard)
library(dashboardthemes)
library(shinyjs)
library(tidyverse)
library(shinycssloaders)
library(ggplot2)
library(plotly)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(DESeq2)
library(ggpubr)
library(DT)
library(reshape2)
library(dplyr)
library(viridis)
library(IsoformSwitchAnalyzeR)
library(DEXSeq)
library(ComplexHeatmap)
`%notin%` <- Negate(`%in%`)

options(shiny.maxRequestSize=30*1024^2)

ui <- dashboardPage(
    dashboardHeader(title="DASiRe"),
    dashboardSidebar(
        sidebarMenu(
            id="tabs",
            menuItem("Start", 
                     tabName="start", 
                     icon=icon("cat")),
            menuItem("Import",
                     tabName = "import",
                     icon=icon('spinner')
            ),
            uiOutput("rnaseq_page"),
            uiOutput("chipseq_page")
        )
    ),
    dashboardBody(
        shinyDashboardThemes(
            theme = "blue_gradient"
        ),
        useShinyjs(),
        tabItems(
            tabItem(
                tabName = "start",
                fluidRow(
                    div(
                        id="start_panel",
                        column(12,
                               h2("Welcome!")
                        )
                    )
                )
            )
            ,
            tabItem(
                tabName="import",
                fluidRow(
                    div(
                        id="import_panel",
                        fluidRow(
                            column(12,
                                   box(
                                       title="Select your data",
                                       div(
                                           checkboxInput(inputId="useexamples", label="Use example data", value=TRUE),
                                           fileInput(
                                               inputId = "rnaseq_csv",
                                               label="Choose RNAseq data you want to run.",
                                               accept = "text/csv"
                                           ),
                                           checkboxInput("header_rna", "Header", value=TRUE),
                                           textInput("rnaseq_sep", "Delimiter", value="\t", placeholder="Type in delimiter."),
                                           fileInput(
                                               inputId = "rnaseq_meta",
                                               label="Choose meta data for your data.",
                                               accept="text/csv"
                                           ),
                                           checkboxInput("header_rnameta", "Header", value=TRUE),
                                           textInput("rnaseq_meta_sep", "Delimiter", value="\t", placeholder="Type in delimiter."),
                                       ),
                                       p("Some bugs here, press upload two times so you can choose sample columns"),
                                       div(
                                           actionButton("renderimport_rna", label="upload", icon=icon("file-import"))
                                       )
                                   ),
                                   box(
                                       title="Select CHIP-seq data you want to run.",
                                       div(
                                           fileInput(inputId = "chipseq_csv", label="Choose your CHIPseq data", accept=".csv"),
                                           actionButton("renderimport_chip", label="upload", icon=icon("file-import") )
                                       )
                                   )
                            ),
                            column(12,
                                   box(
                                       title="Select annotation",
                                       div(
                                           checkboxInput("gtf", "Use GENCODE gtf.", value=TRUE)
                                       )
                                   )
                            )
                        )
                    )
                )
            ),
            tabItem(
                tabName="deseq2",
                uiOutput("deseq2_panels")
            ),
            tabItem(
                tabName="rna_qc",
                uiOutput("rna_qc_panels")
            ),
            tabItem(
                tabName = "rna_iso",
                uiOutput("rna_iso_panels")
            ),
            tabItem(
                tabName = "rna_exon",
                uiOutput("rna_exon_panels")
            ),
            tabItem(
                tabName="rna_splice",
                uiOutput("rna_splice_panels")
            )
        )
    )
)




server <- function(input, output, session) {
    ##render UI    
    output$rnaseq_page <- renderMenu({
        #show menu only if rnaseq data uploaded
        if(input$renderimport_rna == 0) return()
        
        #ui
        sidebarMenu(
            menuItem("RNAseq",
                     icon=icon("dna"),
                     tabName = "rnaseq",
                     menuItem("Quality control",
                              tabName = "rna_qc"),
                     menuItem("Differential expression",
                              tabName="deseq2"),
                     menuItem("Differential exon",
                              tabName="rna_exon"),
                     menuItem("Isoform switch",
                              tabName="rna_iso"),
                     menuItem("Splicing event",
                              tabName="rna_splice")
            )
        )
    })
    
    output$chipseq_page <- renderMenu({
        #show menu only if chipseq data uploaded
        if(input$renderimport_chip == 0) return()
        
        #ui
        sidebarMenu(
            menuItem("CHIP-seq",
                     tabName="chipseq",
                     icon=icon("dna"),
                     menuItem("Quality control"),
                     menuItem("Peak enrichment",
                              menuItem("Promoter level"),
                              menuItem("Gene level")))
        )
    })
    
    output$deseq2_panels <- renderUI({
            fluidRow(
                column(
                    box(
                        h4("Select parameters:"),
                        selectInput("rna_meta_var1", "Select condition (choose 'time')", choices=c()),
                        selectizeInput("rna_meta_var2", "Select a variable", choices=c()),
                        checkboxGroupInput("rna_samples", "Select samples for analysis.", choices = c()),
                        actionButton("load_deseq2", "Load Analysis", icon = icon("play-circle")),
                    ),
                    width=4
                ),
                tabBox(
                    title = "Visualization",
                    width=8,
                    tabPanel("PCA",
                             box(
                                 withSpinner(
                                     plotlyOutput("rna_pca"), type=4
                                 )
                             ),
                             box(
                                 withSpinner(
                                     plotOutput("rna_heatmap"), type=4
                                 )
                             )
                             
                    ),
                    tabPanel("gene",
                             box(
                                 selectizeInput("gene_name", label = "Gene", choices = c()),
                                 withSpinner(
                                     plotOutput("rna_genecount"), type=4
                                 )
                             ),
                    ),
                    tabPanel("DESeq2 result",
                             box(
                                 div(
                                     selectInput("deseq_result", label = "Result", choices = c()),
                                     withSpinner(
                                         plotOutput("rna_volcano"), type=4
                                     )
                                 )
                             ),
                             div(
                                 withSpinner(
                                     DT::dataTableOutput("rna_deseq_table"), type=4
                                 )
                             )
                    )
                )
            )
    })
    
    output$rna_qc_panels <- renderUI({
        fluidRow(
            column(
                box(
                    actionButton("load_qc", "Load Analysis", icon = icon("play-circle"))
                ),
                width=3
            ),
            tabBox(
                id="Quality control",
                width=10,
                tabPanel("STAR output",
                         box(
                            plotOutput("star_qc_percentplot")
                        ),
                        box(
                            plotOutput("star_qc_readlenplot")
                        )
                ),
                tabPanel("Kallisto",
                         box(
                             plotOutput("kallisto_percentplot")
                        )
                )
            )
        )
    })
    
    output$rna_iso_panels <- renderUI({
        fluidRow(
            column(
                box(
                    actionButton("load_iso", "Load Analysis", icon = icon("play-circle"))
                ),
                width=3
            ),
            tabBox(
                id="Isoform switch analysis",
                width=9,
                tabPanel("Switch summary",
                         box(
                             plotOutput("switch_sum")
                         )),
                tabPanel("Switch plots",
                         selectInput("iso_gene_name", "Select a gene", choices=c()),
                         withSpinner(
                            plotOutput("switch_plot"), type=4
                         )
                )
            )
        )
    })
    
    output$rna_exon_panels <- renderUI({
        fluidRow(
            column(
                box(
                    h4("Select parameters:"),
                    textInput("exons_pval_thres", "Select a p-value threshold", value = 0.1),
                    textInput("exons_fc_thres", "Select a fold change threshold", value = 1.5),
                    actionButton("load_exon", "Load Analysis", icon = icon("play-circle"))
                ),
                width=4
            ),
            tabBox(
                title="DEXSeq",
                width=8,
                tabPanel("Summary of DEXSeq",
                         fluidRow(
                             box(
                                 withSpinner(
                                    textOutput("exons_evaluated"),type=4,size=0.2, proxy.height="5px"
                                    ),
                                 p("Number of exons evaluated")
                             ),
                             box(
                                 withSpinner(
                                    textOutput("exons_de"),type=4,size=0.2, proxy.height="5px"),
                                 p("Number of exons differentially expressed")),
                             box(
                                 withSpinner(
                                     textOutput("exons_gene"),type=4,size=0.2, proxy.height="5px"),
                                 p("Number of genes evaluated")
                             ),
                             box(
                                 withSpinner(
                                    textOutput("exons_gene_de"),type=4,size=0.2, proxy.height="5px"),
                                 p("Number of genes with at least one exon differentially expressed")
                             )
                         ),
                         fluidRow(
                             box(
                                 selectizeInput("exon_gene_name", label="Gene", choices=c()),
                                 withSpinner(
                                    plotOutput("exons_per_gene"), type=4
                                 )
                             )
                         )
                )
            )
        )
    })
    
    output$rna_splice_panels <- renderUI({
        fluidRow(
            column(
                box(
                    actionButton("load_majiq", "Load Analysis", icon = icon("play-circle"))
                ),
                width=3
            ),
            tabBox(
                title="Splicing events",
                tabPanel("Number of events",
                         box(
                            plotOutput("event_sum")
                         ),
                         box(
                             plotOutput("majiq_heatmap")
                         )
                )
            )
            
        )
    })
        # fluidRow(
        #     id="deseq2_panel",
        #     box(
        #         verbatimTextOutput("check_meta")
        #     )
        # )
    
    #################### data handler / modal ##########################
    import_data_modal <- function(failed=F){
        if (failed==F){
            if (input$useexamples == TRUE){
                modalDialog(
                    h3("Confirm using examples dataset."),
                    p("If you are importing your own data, please click 'cancel' and uncheck 'Use example data'."),
                    selectInput("rnaseq_metacol", "Automatic sample column:", selected = "sampleID", choices=c("sampleID")),
                    selectInput("rnaseq_genecol", "Automatic sample column:", selected = "ensembl_gene_id_version", choices=c("ensembl_gene_id_version")),
                    footer = tagList(
                        modalButton("Cancel"),
                        actionButton("ok_import_rna", "OK")
                    )
                )
            } else {
                modalDialog(
                    h4("Select gene column for RNA-seq table: (Now you can ignore this)"),
                    selectInput("rnaseq_genecol", "", choices=c()),
                    h4("Select sample column for meta table: (choose 'Run')"),
                    selectInput("rnaseq_metacol", "", choices=c()),
                    footer = tagList(
                        modalButton("Cancel"),
                        actionButton("ok_import_rna", "OK")
                    )
                )
            }
        } else {
            modalDialog(
                h2("Failed to import.")
            )
        }
    }
    
    observeEvent(input$renderimport_rna, {
        showModal(import_data_modal())
    })
    
    deseqlevel <- function(failed=F){
        if (is.null(input$rna_meta_var1)){
            modalDialog(
                h4("Failed to set base level"),
                p("Please input the condition column name.")
            ) 
        } else {
            modalDialog(
                selectInput("rna_deseq_level", "Select a condition as base level for DESeq2.", choices=unique(rnameta_df()[,input$rna_meta_var1])),
                footer = tagList(
                    actionButton("ok_base_level", "OK")
                )
            )
        }
    }
    
    # observeEvent(input$run_pca_rna, {
    #     showModal(deseqlevel())
    # })
    
    #################### reactive object #####################3
    ##show import dialog
    rnaseq_df <- reactive({
        rnafile <- input$rnaseq_csv
        if (input$useexamples == TRUE){
            df <- make_gene_matrix("../examples/gene_counts")
        } else {
            if (is.null(rnafile)) {
                showModal(import_data_modal(failed = TRUE))
                return(NULL)
            }
            df <- read.csv(rnafile$datapath, sep=input$rnaseq_sep, header=input$header_rna)

        }
        return(df)
    })
    
    
    # import meta data
    rnameta_df <- reactive({
        metafile <- input$rnaseq_meta
 
        if (input$useexamples  == TRUE){
            df <- read.csv("../examples/metadata.txt", sep=" ", header=T)
        } else {
            if (is.null(metafile)){
                showModal(import_data_modal(failed = TRUE))
                return(NULL)
            } 
            df <- read.csv(metafile$datapath, sep=input$rnaseq_meta_sep, header=input$header_rnameta)
        }
        return(df)
    })
    
    rnaseq_df_clean <- reactive({
        seqdf <- rnaseq_df()
        row.names(seqdf) <- seqdf[,input$rnaseq_genecol]
        return(seqdf %>% dplyr::select(-input$rnaseq_genecol))
    })
    
    rnaseq_meta_clean <- reactive({
        metadf <- rnameta_df()
        row.names(metadf) <- metadf[,input$rnaseq_metacol]
        return(metadf)
    })
    
    ## filter meta data 
    filtered_meta <- reactive({
        if (input$load_deseq2 == 0){return()}
        metadf <- rnaseq_meta_clean()
        samplename <- input$rnaseq_metacol
        sub <- metadf %>% dplyr::select(samplename) 
        return(metadf[sub[,1] %in% input$rna_samples,])
        
    })

    ## generate dds object 
    #### TODO add covariates
    dds_obj <- reactive({
        if (input$useexamples == TRUE){
            load(file = "../examples/gene_counts/dds.RData")
        } else {
            seqdf <- rnaseq_df_clean()
            metadf <- filtered_meta()
            
            dds <- DESeqDataSetFromMatrix(countData=as.matrix(seqdf[,input$rna_samples]),
                                          colData=metadf[input$rna_samples,],
                                          design=formula(c("~", input$rna_meta_var1)))
            
            dds <- DESeq(dds)
        }
        dds 
    })
    # 
    # dds_result <- reactive({
    #     if (input$ok_deseq_level){
    # 
    #         dds <- dds_obj()
    #         dds[,input$rna_meta_var1] <- relevel(dds[,input$rna_meta_var1], input$rna_deseq_level)
    #         dds <- DESeq(dds)
    #         return(dds)
    #     }
    # })

    ### make reactive pca plot
    rnapcaplot <- eventReactive(input$load_deseq2, {
        if (input$load_deseq2 == 0){return()}
        
        # seqdf <- rnaseq_df_clean()
        # metadf <- filtered_meta()
        samplenames <- input$rnaseq_metacol
        
        # pr <- prcomp(t(seqdf[,input$rna_samples]))
        # pca <- data.frame(pr$x)
        # pca$samples <- row.names(pca)
        # allpca <- inner_join(pca, metadf, by=c("samples"=samplenames))
        
        dds <- dds_obj()
        vst <- vst(dds, blind=FALSE)
        PCA_plot <- plotPCA(vst, intgroup=c(input$rna_meta_var1))
        # p<-ggplot(allpca, aes(PC1, PC2, shape=input$rna_meta_var1, color=input$rna_meta_var2))+
        #     geom_point(size=2.5)
        PCA_plot +
            theme_classic(base_size = 12)+
            theme(#legend.position=c(0.1, 0.9),
                aspect.ratio = 1,
                plot.title = element_text(size = 12),#color = "red",, face = "bold"),
                axis.text.x = element_text(size = 8),
                axis.title.x = element_text(size = 10),
                axis.text.y = element_text(size = 8),
                axis.title.y = element_text(size = 10),
                legend.position = "left",
                axis.ticks = element_blank(),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "transparent",colour = NA),
                plot.background = element_rect(fill = "transparent",colour = NA),
                legend.background = element_rect(fill = "transparent", colour = NA),
                legend.box.background = element_rect(fill = "transparent", colour = NA))
        
    })
    
    output$rna_pca <- renderPlotly({
        rnapcaplot()
    })
    
    ## heatmap rnaseq
    rnaheatmap <- eventReactive(input$load_deseq2, {
        if (input$load_deseq2==0){return()}
        dds <- dds_obj()
        vst <- vst(dds, blind=FALSE)
        sampleDists <- dist(t(assay(vst)))
        sampleDistsMatrix <- as.matrix(sampleDists)
        
        colors <- colorRampPalette( rev(brewer.pal(length(dds$sample), "Blues")) )(255)
        
        
        pheatmap(sampleDistsMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors,legend = FALSE,main = "Sample distances",treeheight_row = 10,treeheight_col = 20)
    })
    
    output$rna_heatmap <- renderPlot({
        rnaheatmap()
        })
    
    check_gene <- eventReactive(input$load_deseq2, {
        req(input$rna_meta_var1)
        dds<-dds_obj()
        gene_counts <- plotCounts(dds, 
                                  gene=input$gene_name, 
                                  intgroup=input$rna_meta_var1, 
                                  returnData=TRUE)
        gene_counts
    })
    
    
    ## gene plot
    rna_genecount_plot <- eventReactive(input$gene_name, {
        req(input$rna_meta_var1)
        dds<-dds_obj()
        gene_counts <- plotCounts(dds, 
                                  gene=input$gene_name, 
                                  intgroup=input$rna_meta_var1, 
                                  returnData=TRUE)
        
        ggplot(gene_counts, aes(x = input$rna_meta_var1, y = count, col=input$rna_meta_var1)) +
            geom_point() +
            scale_color_brewer(palette = "Set1",)+
            theme_classic(base_size = 12)+
            # ggtitle(label = i,subtitle = element_blank())+
            theme(#legend.position=c(0.9, 0.9),
                # aspect.ratio = 1,
                # axis.text.x = element_text(size = 8),
                # axis.title.x = element_text(size = 10),
                # axis.text.y = element_text(size = 8),
                # axis.title.y = element_text(size = 10),
                legend.position = "none",
                axis.ticks = element_blank(),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "transparent",colour = NA),
                plot.background = element_rect(fill = "transparent",colour = NA),
                legend.background = element_rect(fill = "transparent", colour = NA),
                legend.box.background = element_rect(fill = "transparent", colour = NA))+
            scale_y_log10()+
            labs(y="log10(normalized_counts)",col="Condition")
        
    })
    
    output$rna_genecount <- renderPlot({
        rna_genecount_plot()
    })
    
    rna_volcano_plot <- eventReactive(input$load_deseq2, {
        dds <- dds_obj()
        req(input$deseq_result)
        
        
        deg <- data.frame(results(dds, name=input$deseq_result))
        
        ggplot(data=deg,aes(x=log2FoldChange,y=-log10(padj))) +
            geom_vline(xintercept=c(-log2(1.5),log2(1.5)), color="red")+ 
            geom_hline(yintercept=-log10(0.05), color="blue")+ 
            geom_point(color="black",alpha=0.5,stat="identity")+
            xlim(-10,10)
    })
    
    output$rna_volcano <- renderPlot({
        rna_volcano_plot()
    })
        
    rna_deseq_res_table <- eventReactive(input$load_deseq2, {
        dds <- dds_obj()
        req(input$deseq_result)
        
        
        deg <- data.frame(DESeq2::results(dds, name=input$deseq_result))
        deg
    })
    
    output$rna_deseq_table <- DT::renderDataTable({rna_deseq_res_table()})
    
    ######################### isoform switch ###########################
    
    isoform_data <- reactive({
        if (input$useexamples == TRUE){
            load(file = "../examples/pseudocounts/exampleSwitchListAnalyzed.RData")
            datatable(extractSwitchSummary(exampleSwitchListAnalyzed))
        }
        return(exampleSwitchListAnalyzed)
    })
    
    isoform_sig_genes <- reactive({
        exampleSwitchListAnalyzed <- isoform_data()
        
        signifcant_genes<- extractTopSwitches(
            exampleSwitchListAnalyzed,
            extractGenes = TRUE,
            filterForConsequences = FALSE,
            n = NA,
            sortByQvals = TRUE
        )
        return(signifcant_genes)
    })
    
    switch_plotplot <- eventReactive(input$load_iso, {
        req(input$iso_gene_name)
        exampleSwitchListAnalyzed <- isoform_data()
        isa_genes <- isoform_sig_genes()
        
        switchPlot(exampleSwitchListAnalyzed,
                   gene=input$iso_gene_name)
        
    })
    
    iso_switch_sumplot <- eventReactive(input$load_iso, {
        exampleSwitchListAnalyzed <- isoform_data()
        
        extractSplicingSummary(
            exampleSwitchListAnalyzed,
            asFractionTotal = FALSE,
            plotGenes=FALSE
        )
        
    })

    
    output$switch_sum <- renderPlot({
        iso_switch_sumplot()
    })
    
    output$switch_plot <- renderPlot({
        switch_plotplot()
    })
    
    ######################### DEXSeq ################################
    
    exon_data <- reactive({
        if (input$useexamples==TRUE){
            load(file = "../examples/exon_counts/dxr1.RData")
        }
        return(dxr1)
    })


    sig_exon_data <- reactive({
        if (is.null(exon_data())) {return()}
        dxr1 <- exon_data()

        significant_exons <- as.data.frame(dxr1) %>%
            dplyr::filter(padj < input$exons_pval_thres)
        return(significant_exons)

    })

   observeEvent(input$load_exon, {
        if (is.null(exon_data())){return()}
        
        output$exons_evaluated <- renderText({
            dxr1 <- exon_data()
            return(sum(table ( dxr1$padj < input$exons_pval_thres )))
        })
        
        output$exons_de <- renderText({
            dxr1 <- exon_data()
            return(table ( dxr1$padj < input$exons_pval_thres )[2])
        })
   
        output$exons_gene <- renderText({
            dxr1 <- exon_data()
            return(sum(table(tapply( dxr1$padj < input$exons_pval_thres, dxr1$groupID, any))))
        })
        output$exons_gene_de <- renderText({
            dxr1 <- exon_data()
            return(table(tapply(dxr1$padj < input$exons_pval_thres, dxr1$groupID, any))[2])
        })
    })
   
   observeEvent(input$load_exon, {
       sig_exons <- sig_exon_data()
       updateSelectizeInput(session, "exon_gene_name", label="Select a gene", choices=unique(sig_exons$groupID), server=TRUE)
   })
   
    exons_per_gene_plot <- eventReactive(input$exon_gene_name, {
        if (is.null(exon_data())){return()}
        dxr1 <- exon_data()
        plotDEXSeq(dxr1, gene=input$gene_name, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
    })
    
    output$exons_per_gene <- renderPlot({exons_per_gene_plot()})
    
    ########################## MAJIQ ####################################
    
    ################
    ## Exon Skipping
    ################
    
    voila_res_all <- reactive({
        if (input$useeamples==TRUE){
            df<-read.delim(file = "../examples/majiq_output/voila_results_all.tsv",header = TRUE,sep = "\t")
        }
        return(df)
    })
    
    voila_res <- reactive({
        if(input$useexamples==TRUE){
            df <- read.delim(file = "../examples/majiq_output/voila_results.tsv",header = TRUE,sep = "\t")
        } 
        return(df)
    })
    
    majiq_es_all <- reactive({
        df <- voila_res_all()
        df %>%
            dplyr::filter(ES =="True" & A5SS =="False" & A3SS =="False" & IR.coords=="")  %>%
            dplyr::filter(!grepl('na', LSV.ID))
    })
    
    majiq_es_pos <- reactive({
        df <- voila_res()
        df %>%
            dplyr::filter(ES =="True" & A5SS =="False" & A3SS =="False" & IR.coords=="")  %>%
            dplyr::filter(!grepl('na', LSV.ID))
    })
    
    ## Get POSITIVE ES events
    majiq.es.pos <- reactive({
        df <- voila_res()
        df %>%
        dplyr::filter(ES =="True" & A5SS =="False" & A3SS =="False" & IR.coords=="")  %>%
        dplyr::filter(!grepl('na', LSV.ID))
    })
    
    ## Get NEGATIVE ES events
    majiq.es.neg <- reactive({
        df <- majiq_es_all()
        temp <- data.frame(do.call("rbind", strsplit(x = df$E.dPSI..per.LSV.junction,split = ";",fixed = TRUE)))
        temp  <- as.data.frame(apply(temp, 2, as.numeric))
        keep <- apply(temp < 0.001 & temp > -0.001, 1, all, na.rm=TRUE)
        temp  <- as.data.frame(cbind(temp, df$LSV.ID))
        df[keep,]
    })
    
    ################
    ## A3SS
    ################
    
    majiq.a3.all <- reactive({ 
        df <- voila_res_all() 
        df%>%
        dplyr::filter(ES =="False" & A5SS =="False" & A3SS =="True" & IR.coords=="")  %>%
        dplyr::filter(!grepl('na', LSV.ID))
    })
    
    ## Get POSITIVE A3 events
    majiq.a3.pos <- reactive({
        df <- voila_res() 
        df %>%
        dplyr::filter(ES =="False" & A5SS =="False" & A3SS =="True" & IR.coords=="")  %>%
        dplyr::filter(!grepl('na', LSV.ID))
    })
    
    ## Get NEGATIVE A3 events
    majiq.a3.neg <- reactive({
        df <- majiq.a3.all()
        temp <- data.frame(do.call("rbind", strsplit(x = df$E.dPSI..per.LSV.junction,split = ";",fixed = TRUE)))
        temp  <- as.data.frame(apply(temp, 2, as.numeric))
        keep <- apply(temp < 0.001 & temp > -0.001, 1, all, na.rm=TRUE)
        temp  <- as.data.frame(cbind(temp, df$LSV.ID))
        df[keep,]
    })
    
    ################
    ## A5SS
    ################
    
    ## Get all A5 events
    majiq.a5.all <- reactive({
        df <- voila_res_all() 
        df %>%
        dplyr::filter(ES =="False" & A5SS =="True" & A3SS =="False" & IR.coords=="")  %>%
        dplyr::filter(!grepl('na', LSV.ID))
    })
    
    ## Get POSITIVE A5 events
    majiq.a5.pos <- reactive({
        df <- voila_res() 
        df %>%
        dplyr::filter(ES =="False" & A5SS =="True" & A3SS =="False" & IR.coords=="")  %>%
        dplyr::filter(!grepl('na', LSV.ID))
    })
    
    ## Get NEGATIVE A5 events
    majiq.a5.neg <- reactive({
        df <- majiq.a5.all()
        temp <- data.frame(do.call("rbind", strsplit(x = df$E.dPSI..per.LSV.junction,split = ";",fixed = TRUE)))
        temp  <- as.data.frame(apply(temp, 2, as.numeric))
        keep <- apply(temp < 0.001 & temp > -0.001, 1, all, na.rm=TRUE)
        temp  <- as.data.frame(cbind(temp, df$LSV.ID))
        df[keep,]
    })
    
    ################
    ## IR
    ################
    
    ## Get all IR events
    majiq.ir.all <- reactive({
        df <- voila_res_all() 
        df %>%
        dplyr::filter(ES =="False" & A5SS =="False" & A3SS =="False" & IR.coords!="")  %>%
        dplyr::filter(!grepl('na', LSV.ID))
    })
    
    ## Get POSITIVE IR events
    majiq.ir.pos <- reactive({
        df <- voila_res() 
        df %>%
        dplyr::filter(ES =="False" & A5SS =="False" & A3SS =="True" & IR.coords!="")  %>%
        dplyr::filter(!grepl('na', LSV.ID))
    })
    
    ## Get NEGATIVE IR events
    majiq.ir.neg <- reactive({
        df <- majiq.ir.all()
        temp <- data.frame(do.call("rbind", strsplit(x = df$E.dPSI..per.LSV.junction,split = ";",fixed = TRUE)))
        temp  <- as.data.frame(apply(temp, 2, as.numeric))
        keep <- apply(temp < 0.001 & temp > -0.001, 1, all, na.rm=TRUE)
        temp  <- as.data.frame(cbind(temp, df$LSV.ID))
        df[keep,]
    })
    
    
    ##more than one event 
    majiq.moreThanOneEvent.pos <- reactive({
        
        
        df <- voila_res()
        majiq.es.pos.df <- majiq.es.pos()
        majiq.a3.pos.df <- majiq.a3.pos()
        majiq.a5.pos.df <- majiq.a5.pos()
        majiq.ir.pos.df <- majiq.ir.pos()
        df %>%
        dplyr::filter(LSV.ID %notin% c(majiq.es.pos.df$LSV.ID, majiq.a3.pos.df$LSV.ID, majiq.a5.pos.df$LSV.ID, majiq.ir.pos.df$LSV.ID))
    })
    
    majiq.NrGenes.df <- reactive({
        majiq.NrGenes <- data.frame(rbind(
            c("ES",nrow(majiq.es.pos())),
            c("IR",nrow(majiq.ir.pos())),
            c("A3SS",nrow(majiq.a3.pos())),
            c("A5SS",nrow(majiq.a5.pos())),
            c("Complex",nrow(majiq.moreThanOneEvent.pos()))
        ))
        
        colnames(majiq.NrGenes) <- c("Event type","Number of genes")
        majiq.NrGenes$`Number of genes` <- as.numeric(majiq.NrGenes$`Number of genes`)
        majiq.NrGenes$fraction <- majiq.NrGenes$`Number of genes` / sum(majiq.NrGenes$`Number of genes`)
        majiq.NrGenes$ymax <- cumsum(majiq.NrGenes$fraction)
        majiq.NrGenes$ymin <- c(0, head(majiq.NrGenes$ymax, n=-1))
        majiq.NrGenes$labelPosition <- (majiq.NrGenes$ymax + majiq.NrGenes$ymin) / 2
        majiq.NrGenes$label <- paste0(majiq.NrGenes$`Event type`, "\n", majiq.NrGenes$`Number of genes`)
        
        return(majiq.NrGenes)
    })
    
    event_sum_plot <- reactive({
        majiq.NrGenes <- majiq.NrGenes.df()
        ggplot(majiq.NrGenes, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=`Event type`)) +
            geom_rect() +
            geom_label( x=c(3.5,3.0,3.5,4,3), aes(y=labelPosition, label=label), size=5) +
            scale_fill_brewer(palette=4) +
            coord_polar(theta="y") +
            xlim(c(2, 4)) +
            theme_void() +
            theme(legend.position = "none")
        
    })
    
    output$event_sum <- renderPlot({event_sum_plot()})
    
    majiq.heatmap.df <- reactive({
        majiq.moreThanOneEvent.pos.df <- majiq.moreThanOneEvent.pos()
        majiq.heatmap <- majiq.moreThanOneEvent.pos.df %>% 
            dplyr::select(A5SS,A3SS,ES,IR.coords) %>%
            mutate(IR.coords=if_else(IR.coords!="", true="True",false="False"))
        
        majiq.heatmap <- as.data.frame(table(majiq.heatmap))
        
        majiq.heatmap <- majiq.heatmap %>%
            dplyr::filter(!(A3SS=="True" & A5SS=="False" & IR.coords=="False" & ES=="False"))%>%
            dplyr::filter(!(A3SS=="False" & A5SS=="True" & IR.coords=="False" & ES=="False"))%>%
            dplyr::filter(!(A3SS=="False" & A5SS=="False" & IR.coords=="True" & ES=="False"))%>%
            dplyr::filter(!(A3SS=="False" & A5SS=="False" & IR.coords=="False" & ES=="True"))%>%
            dplyr::filter(!(A3SS=="False" & A5SS=="False" & IR.coords=="False" & ES=="False"))
        
        return(majiq.heatmap)
    })
    
    majiq_heatmap_plot <- reactive({
        majiq.heatmap <- majiq.heatmap.df()
        majiq.heatmap <- t(majiq.heatmap)
        rownames(majiq.heatmap) <- c("A5SS","A3SS","ES","IR","Freq")
        
        col_ha = HeatmapAnnotation(NrGenes = anno_barplot(as.numeric(majiq.heatmap[5,]),
                                                          add_numbers=TRUE,
                                                          height = unit(7,"cm"),gp=gpar(fill=RColorBrewer::brewer.pal(n = 3,name = "Dark2")[3])),
                                   Genes = anno_simple(as.numeric(majiq.heatmap[5,]),pch = majiq.heatmap[5,]))

        Heatmap(matrix = majiq.heatmap[1:4,],
                col = c("False"="White","True"="grey"),#RColorBrewer::brewer.pal(n = 3,name = "Dark2")[3]),
                border = TRUE,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                height = unit(3,"cm"),
                name = "mat", 
                show_heatmap_legend = FALSE,
                width = unit(11,"cm"),
                top_annotation = col_ha,
                cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.rect(x = x, y = y, width = width, height = height, 
                              gp = gpar(col = "black", fill = NA))})
            
    })
    
    output$majiq_heatmap <- renderPlot({
        majiq_heatmap_plot()
    })
    
    ######################### observer #################################
    
    observeEvent(input$renderimport_rna, {
        if (is.null(rnaseq_df()) | is.null(rnameta_df())){return()}
        if (input$useexamples!=TRUE){
            updateSelectInput(session, "rnaseq_genecol", "Select the gene column", choices=colnames(rnaseq_df()))
            updateSelectInput(session, "rnaseq_metacol", "Select the sample names column", choices=colnames(rnameta_df()))
        } else {
            updateTabItems(session, "RNAseq", selected="deseq2")
            
        }
    })
    
    observeEvent(input$load_iso, {
        x <- isoform_sig_genes()
        updateSelectizeInput(session, "iso_gene_name", "Select a gene with switches", choices=unique(x$gene_name), server=TRUE)
    })
    
    observeEvent(input$ok_import_rna,{
        removeModal()
    })
    
    observe({
        if (input$tabs=="deseq2") {
            updateSelectInput(session, "rna_meta_var1", "Select group column", choices=colnames(rnameta_df()))
            updateSelectizeInput(session, "rna_meta_var2", "Select group column", choices=colnames(rnameta_df()))
            updateCheckboxGroupInput(session, "rna_samples", "Select samples for analysis", choices=colnames(rnaseq_df()), selected=colnames(rnaseq_df()))
        } else {
            updateSelectInput(session, "rna_meta_var1", "Select group column", choices="check again lalalal")
        }
    })
    
    ## gene selection
    # observeEvent(input$rnaseq_genecol,{
    #     req(input$rnaseq_genecol)
    #     if (is.null(dds_obj())) {return()}
    #     dds <- dds_obj()
    #     updateSelectizeInput(session, "gene_name", "Select gene", choices=row.names(counts(dds)), server = TRUE)
    # })
    
    observeEvent(input$load_deseq2, {
        if (input$load_deseq2==0){return()}
        dds <- dds_obj()
        updateSelectInput(session, "deseq_result", "Choose a comparison", choices=resultsNames(dds))
        updateSelectizeInput(session, "gene_name", "Select gene", choices=row.names(counts(dds)), server = TRUE)
    })


    
    ##################### quality control #########################
    
    # star qc
    star_qc <- reactive({
        if (input$useexamples == TRUE){
            alignment_stats_star <- read.table(file = "../examples/multiqc_data/multiqc_star.txt",header = TRUE,sep = "\t")
            alignment_stats_star <- melt(data = alignment_stats_star,id.vars="Sample")#,"total_reads","avg_input_read_length"))
        }
        return(alignment_stats_star)
    })
    
    star_qc_plot <- eventReactive(input$load_qc, {
        alignment_stats_star <- star_qc()
        alignment_stats_star$value <- as.numeric(alignment_stats_star$value)
        
            alignment_stats_star %>%
                dplyr::filter(grepl(pattern = "percent",x = alignment_stats_star$variable)) %>%
                ggplot(aes(x = value, y = Sample)) +
                geom_bar(aes(fill = variable), stat="identity") + 
                scale_fill_viridis(discrete = TRUE,alpha = 0.75,option = "D",direction = 1) +
                theme_classic()+
                ggtitle(label = "Percentage of reads aligned with STAR")+
                theme(#legend.position=c(0.1, 0.7),
                    axis.text.x = element_text(angle = 45,hjust = 1),
                    axis.ticks = element_blank(),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "transparent",colour = NA),
                    plot.background = element_rect(fill = "transparent",colour = NA),
                    legend.background = element_rect(fill = "transparent", colour = NA),
                    legend.box.background = element_rect(fill = "transparent", colour = NA))

    })
    
    star_qc_readlength <- eventReactive(input$load_qc, {
        alignment_stats_star <- star_qc()
        alignment_stats_star$value <- as.numeric(alignment_stats_star$value)
        
        alignment_stats_star %>%
            dplyr::filter(grepl(pattern = "read_length",x = alignment_stats_star$variable)) %>%
            ggplot(aes(x = value, y = variable)) +
            geom_bar(aes(fill=Sample),position = "dodge",stat="identity") + 
            scale_fill_viridis(discrete = TRUE,alpha = 0.75,option = "D",direction = 1) +
            theme_classic()+
            ggtitle(label = "Percentage of reads aligned with STAR")+
            theme(#legend.position=c(0.1, 0.7),
                axis.text.x = element_text(angle = 45,hjust = 1),
                axis.ticks = element_blank(),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "transparent",colour = NA),
                plot.background = element_rect(fill = "transparent",colour = NA),
                legend.background = element_rect(fill = "transparent", colour = NA),
                legend.box.background = element_rect(fill = "transparent", colour = NA))
        
    })
        
    output$star_qc_percentplot <- renderPlot({
        star_qc_plot()
    })
    
    output$star_qc_readlenplot <- renderPlot({
        star_qc_readlength()
    })
    
    pseudo_qc <- eventReactive(input$load_qc, {
        if (input$useexamples == TRUE){
            alignment_stats_kallisto <- read.csv("../examples/multiqc_data/multiqc_kallisto.txt",header = TRUE,sep = "\t")
            alignment_stats_kallisto <- melt(data = alignment_stats_kallisto,id.vars="Sample")#,"total_reads","avg_input_read_length"))
        }
        
        return(alignment_stats_kallisto)
    })
    
    pseudo_qc_plot <- eventReactive(input$load_qc, {
        alignment_stats_kallisto <- pseudo_qc()
        alignment_stats_kallisto %>%
            dplyr::filter(grepl(pattern = "percent",x = alignment_stats_kallisto$variable)) %>%
            ggplot(aes(x = value, y = Sample)) +
            geom_bar(aes(fill = variable), stat="identity",) +
            scale_fill_viridis(discrete = TRUE,alpha = 0.75,option = "D",direction = 1)+
            xlim(0,100)+
            theme_classic()+
            ggtitle(label = "Percentage of reads pseudoaligned with Kallisto")+
            theme(#legend.position=c(0.1, 0.7),
                axis.text.x = element_text(angle = 45,hjust = 1),
                axis.ticks = element_blank(),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "transparent",colour = NA),
                plot.background = element_rect(fill = "transparent",colour = NA),
                legend.background = element_rect(fill = "transparent", colour = NA),
                legend.box.background = element_rect(fill = "transparent", colour = NA))
        
    })
    
    output$kallisto_percentplot <- renderPlot({pseudo_qc_plot()})
    
    #################### dynamic rendering ##########################
    # observeEvent("",{
    #     shinyjs::show("start_panel")
    #     shinyjs::hide("import_panel")
    #     shinyjs::hide("rnaseq")
    #     shinyjs::hide("deseq2_panel")
    # }, once=T)
    # 
    # observeEvent(input$tabs,{
    #     if (input$tabs=="import"){
    #         shinyjs::show("import_panel")
    #         shinyjs::hide("start_panel")
    #         shinyjs::hide("deseq2_panel")
    #     } else if (input$tabs=="start"){
    #         shinyjs::show("start_panel")
    #         shinyjs::hide("import_panel")
    #         shinyjs::hide("deseq2_panel")
    #     } 
    # })
    # 

}

shinyApp(ui, server)