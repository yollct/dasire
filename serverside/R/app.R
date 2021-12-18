library(shiny)
library(shinydashboard)
library(dashboardthemes)
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
library(shinyWidgets)
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(genomation)
library(ChIPseeker)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(GenomicFeatures)
library(RMariaDB)
library(biomaRt)
library(Gviz)
library(ChIPpeakAnno)
library(shinyFiles)
library(tippy)
library(VennDiagram)

source("global.R")
`%notin%` <- Negate(`%in%`)

options(shiny.maxRequestSize=30*1024^2, shiny.trace = TRUE, shiny.reactlog=TRUE)

ui <- dashboardPage(
    dashboardHeader(title="DASiRe", titleWidth = 300),
    dashboardSidebar(
        width = 300,
        sidebarMenu(
            id="tabs",
            menuItem("Start", 
                     tabName="start", 
                     icon=icon("cat")),
            menuItem("Documentation",
                     tabName="docs",
                     icon=icon("book")),
            menuItem("Upload data",
                     tabName = "import",
                     icon=icon('spinner')
            ),
            uiOutput("rnaseq_page"),
            uiOutput("chipseq_page"),
            menuItem("About",
                     tabName="about_tab",
                     icon=icon("fish"))
        )
    ),
    dashboardBody(
        shinyDashboardThemes(
            theme = "blue_gradient"
        ),
        tabItems(
            tabItem(tabName="about_tab",
                    h4("Download example data here: "),
                    downloadButton("download_raw_data", "Download data"),
                    uiOutput("about_page")),
                    
            tabItem(tabName = "docs",
                    withSpinner(htmlOutput("doc_page"), type=4),
                    shiny::tags$head(shiny::tags$style(HTML("
                                
                               body {
                                  width: 100% !important;
                                  max-width: 100% !important;
                               }

                               ")))),
            tabItem(
                tabName = "start",
                fluidRow(
                    div(
                        id="start_panel",
                        column(12,
                               withSpinner(htmlOutput("landing_page"), type=4)
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
                                       title="Upload sample data",
                                       div(h4(HTML("<b>Example data</b>"),span(shiny::icon("info-circle"), id = "rna_ex")),
                                           tippy::tippy_this(elementId = "rna_ex",tooltip = "Check to use our example data and press the upload button.",placement = "right"),
                                           checkboxInput(inputId="useexamples", label="Use example data", value=TRUE),
                                           h4(HTML("<b>Choose your directory</b>"),span(shiny::icon("info-circle"), id = "rna_dir")),
                                           tippy::tippy_this(elementId = "rna_dir",tooltip = "Please input the preprocessing output directory",placement = "right"),
                                           shinyDirButton("dir", "Input directory", "Upload"),
                                       div(
                                           actionButton("renderimport_rna", label="upload", icon=icon("file-import"))
                                       )
                                   )),
                                   box(
                                       title="Upload ChIP data",
                                       div(
                                           h4(HTML("<b>Example data</b>"),span(shiny::icon("info-circle"), id = "chip_ex")),
                                           tippy::tippy_this(elementId = "chip_ex",tooltip = "Check to use our example data and press the upload button.",placement = "right"),
                                           checkboxInput(inputId="useexamples_chip", label="Use example data", value=TRUE),
                                           h4(HTML("<b>Choose your directory</b>"),span(shiny::icon("info-circle"), id = "xhip_dir")),
                                           tippy::tippy_this(elementId = "chip_dir",tooltip = "Please input the preprocessing output directory",placement = "right"),
                                           fileInput(inputId = "chipseq_bed", label="Choose your ChIPseq bed file"),
                                           actionButton("renderimport_chip", label="upload", icon=icon("file-import") )
                                       )
                                   )
                            ),
                            column(12,
                                   box(
                                       title="Select annotation",
                                       div(checkboxInput("gtf", "Use GENCODE gtf.", value=TRUE))
                                   )
                                   
                            )
                        )
                    )
                
            )),
            tabItem(
                tabName="deseq2",
                fluidRow(
                    column(12,
                        box(width=12,
                            column(3, h4("Select parameters:"),
                            fluidRow(selectInput("rna_meta_var1", "Select condition (choose 'time')", choices=c())),
                            fluidRow(selectizeInput("rna_meta_var2", "Select a variable", choices=c()))),
                            column(3, textInput("log2cutoff", "Set a log2 Fold Change cut off", value=1.5)),
                            column(3, fluidRow(checkboxGroupInput("rna_samples", "Select samples for analysis.", choices = c())),
                            fluidRow(actionButton("load_deseq2", "Load Analysis", icon = icon("play-circle")))),
                        )
                    )),
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
                fluidRow(
                    column(12,
                        box(width=12,
                            column(3,h4("Select parameters:")),
                            column(3,textInput("exons_pval_thres", "Select a p-value threshold", value = 0.1)),
                            column(3, textInput("exons_fc_thres", "Select a fold change threshold", value = 1.5)),
                            column(3, actionButton("load_exon", "Load Analysis", icon = icon("play-circle")))
                        )
                    )),
                uiOutput("rna_exon_panels")
            ),
            tabItem(
                tabName="rna_splice",
                uiOutput("rna_splice_panels")
            ),
            tabItem(
                tabName="com_splice",
                uiOutput("com_splice_panels")
            ),
            tabItem(
                tabName="chip_qc",
                fluidRow(column(12,box(
                    column(4, selectInput("fileassembly", label="Choose an assembly version", choices=c())),
                    column(4, selectInput("biosample", label="Choose a sample", choices=c())),
                    column(4, actionButton("load_encode", "Load ENCODE DATA")),
                    )
                )),
                uiOutput("chip_qc_panels")
            ),
            tabItem(
                tabName = "peak_enr",
                fluidRow(column(12,box(width=12,
                    column(3, box(textOutput("show_fileassembly"), p("Current assembly version"))),
                    column(3, box(textOutput("show_biosample"), p("Current biosample"))),
                    column(3, actionButton("run_enr", "Run enrichment"))
                )
                )),
                uiOutput("chip_peak_panels")
            )
        )
    )
)





server <- function(input, output, session) {
    ##render UI    
    output$about_page <- renderUI({
        includeHTML("about_page.html")
    })

    output$landing_page <- renderUI({
        includeHTML("landing_page.html")
    })
    
    output$doc_page <- renderUI({
        #HTML(markdown::markdownToHTML(knit("documentation/tutorial.Rmd", quiet = TRUE), fragment.only=TRUE))
        includeHTML("documentation/tutorial.html")
    })

    output$download_raw_data <- downloadHandler(filename=function(){"example_data.zip"},
                                            content = function(file){
                                                file.copy("examples/example_data.zip", file)
                                            })
    

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
                              tabName="rna_splice"),
                    menuItem("Comparative analysis",
                             tabName="com_splice")
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
                     menuItem("Quality control",
                              tabName="chip_qc"),
                     menuItem("Peak enrichment",
                              tabName="peak_enr")
            )
        )
    })
    
    output$deseq2_panels <- renderUI({
        fluidRow(column(12,
            tabBox(
                title = "Visualization",
                width=12,
                height=8,
                tabPanel("PCA",
                         column(6,div(
                             withSpinner(plotlyOutput("rna_pca"), type=4),
                             div(
                                 style = "position: absolute; left: 1em; bottom: 0.5em;",
                                 dropdown(
                                     selectInput(inputId="deseq2_pca_down_ext", label="Download as...", choices=c("png","svg","pdf","eps","jpeg")),
                                     downloadButton(outputId = "deseq2_pca_down", label = "DOWNLOAD"),
                                     size = "xs",
                                     icon = icon("download", class = "opt"), 
                                     up = TRUE
                                 )
                             )
                         )),
                         column(6, div(
                             withSpinner(
                                 plotOutput("rna_heatmap"), type=4
                             ),
                             div(
                                 style = "position: absolute; left: 1em; bottom: 0.5em;",
                                 dropdown(
                                     selectInput(inputId="deseq2_hm_down_ext", label="Download as...", choices=c("png","svg","pdf","eps","jpeg")),
                                     downloadButton(outputId = "deseq2_hm_down", label = "DOWNLOAD"),
                                     size = "xs",
                                     icon = icon("download", class = "opt"), 
                                     up = TRUE
                                 )
                             )
                         ))
                         
                ),
                tabPanel("Gene Normalized Counts",
                         column(6, div(
                             selectizeInput("gene_name", label = "Gene", choices = c()),
                             withSpinner(
                                 plotOutput("rna_genecount"), type=4
                             ),
                             div(
                                 style = "position: absolute; left: 1em; bottom: 0.5em;",
                                 dropdown(
                                     selectInput(inputId="deseq2_gc_down_ext", label="Download as...", choices=c("png","svg","pdf","eps","jpeg")),
                                     downloadButton(outputId = "deseq2_gc_down", label = "DOWNLOAD"),
                                     size = "xs",
                                     icon = icon("download", class = "opt"), 
                                     up = TRUE
                                 )
                             )
                         )),
                ),
                tabPanel("DESeq2 result",
                         fluidRow(
                             div(
                                 selectInput("deseq_result", label = "Result", choices = c()),
                                 withSpinner(
                                     plotOutput("rna_volcano"), type=4
                                 )
                             ),
                             div(
                                 style = "position: absolute; left: 1em; bottom: 0.5em;",
                                 dropdown(
                                     selectInput(inputId="deseq2_vc_down_ext", label="Download as...", choices=c("png","svg","pdf","eps","jpeg")),
                                     downloadButton(outputId = "deseq2_vc_down", label = "DOWNLOAD"),
                                     size = "xs",
                                     icon = icon("download", class = "opt"), 
                                     up = TRUE
                                 )
                             )
                         ),
                         fluidRow(
                             div(withSpinner(
                                 DT::dataTableOutput("rna_deseq_table"), type=4),
                             div(style = "position: absolute; left: 1em; bottom: 0em;",
                                 dropdown(
                                     downloadButton(outputId = "deseq2_table_csv", label = "CSV"),
                                     downloadButton(outputId = "deseq2_table_excel", label = "EXCEL"),
                                     size = "xs",
                                     icon = icon("download", class = "opt"), 
                                     up = FALSE
                                 ))
                             )
                         )
                ),
                tabPanel("Splicing factors differential expression",
                         fluidRow(div(withSpinner(plotlyOutput("sf_deseq2"), type=4),
                                  div(
                                      style = "position: absolute; left: 1em; bottom: 0.5em;",
                                      dropdown(
                                          selectInput(inputId="sf_deseq2_down_ext", label="Download as...", choices=c("png","svg","pdf","eps","jpeg")),
                                          downloadButton(outputId = "sf_deseq2_down", label = "DOWNLOAD"),
                                          size = "xs",
                                          icon = icon("download", class = "opt"), 
                                          up = TRUE
                                      )
                                  ))))
            )
        ))
    })
    
    output$rna_qc_panels <- renderUI({
        fluidRow(column(12,
            tabBox(
                id="Quality control",
                width=12,
                tabPanel("STAR alignment",
                         column(6,div(
                             withSpinner(
                                 plotOutput("star_qc_percentplot"), type=4
                             ),
                             div(
                                 style = "position: absolute; left: 1em; bottom: 0.5em;",
                                 dropdown(
                                     selectInput(inputId="star_pp_down_ext", label="Download as...", choices=c("png","svg","pdf","eps","jpeg")),
                                     downloadButton(outputId = "star_pp_down", label = "DOWNLOAD"),
                                     size = "xs",
                                     icon = icon("download", class = "opt"), 
                                     up = TRUE
                                 )
                             )
                         )),
                         column(6, div(
                             withSpinner(
                                 plotOutput("star_qc_readlenplot"), type=4
                             ), 
                             div(
                                 style = "position: absolute; left: 1em; bottom: 0.5em;",
                                 dropdown(
                                     selectInput(inputId="star_rp_down_ext", label="Download as...", choices=c("png","svg","pdf","eps","jpeg")),
                                     downloadButton(outputId = "star_rp_down", label = "DOWNLOAD"),
                                     size = "xs",
                                     icon = icon("download", class = "opt"), 
                                     up = TRUE
                                 )
                             )
                         ))
                ),
                tabPanel("Kallisto pseudoalignment",
                         width=12,
                         div(
                             width=6,
                             withSpinner(
                                 plotOutput("kallisto_percentplot"), type=4
                             ), 
                             div(
                                 style = "position: absolute; left: 1em; bottom: 0.5em;",
                                 dropdown(
                                     selectInput(inputId="kall_pp_down_ext", label="Download as...", choices=c("png","svg","pdf","eps","jpeg")),
                                     downloadButton(outputId = "kall_pp_down", label = "DOWNLOAD"),
                                     size = "xs",
                                     icon = icon("download", class = "opt"), 
                                     up = TRUE
                                 )
                             )
                         )
                )
            )
        ))
    })
    
    output$rna_iso_panels <- renderUI({
        fluidRow(column(12,
            tabBox(
                id="Isoform switch analysis",
                width=12,
                height=8,
                tabPanel("Genome-wide isoform splicing analysis",
                         fluidRow(div(width=12,
                             withSpinner(plotOutput("switch_sum"),type=4),
                             div(style = "position: absolute; left: 1em; bottom: 0em;",
                                 dropdown(
                                     selectInput(inputId="iso_switch_sum_down_ext", label="Download as...", choices=c("png","svg","pdf","eps","jpeg")),
                                     downloadButton(outputId = "iso_switch_sum_down", label = "DOWNLOAD"),
                                     size = "xs",
                                     icon = icon("download", class = "opt"), 
                                     up = TRUE))
                         )),
                        fluidRow(div(width=8, withSpinner(plotOutput("splice_genom"),type=4)))),
                # tabPanel("Splicing summary",
                #          div(width=12,
                #              withSpinner(plotOutput("switch_enrich"),type=4))),
                tabPanel("Gene Switch plots",
                         fluidRow(
                            selectInput("iso_gene_name", "Select a gene", choices=c()),
                            actionButton("load_iso", "Load"))
                            ,
                         fluidRow(
                            withSpinner(
                                 plotOutput("switch_plot"), type=4
                             ),
                            div(
                                style = "position: absolute; left: 1em; bottom: 0em;",
                                dropdown(
                                    selectInput(inputId="iso_switch_plot_down_ext", label="Download as...", choices=c("png","svg","pdf","eps","jpeg")),
                                    downloadButton(outputId = "iso_switch_plot_down", label = "DOWNLOAD"),
                                    size = "xs",
                                    icon = icon("download", class = "opt"), 
                                    up = TRUE
                                ))
                         )
                ),
                tabPanel("Genes Switch table",
                         withSpinner(
                             DT::dataTableOutput("switch_table"), type=4
                         ),
                         div(
                             style = "position: absolute; left: 1em; bottom: 0em;",
                             dropdown(
                                 downloadButton(outputId = "iso_table_csv", label = "CSV"),
                                 downloadButton(outputId = "iso_table_excel", label = "EXCEL"),
                                 size = "xs",
                                 icon = icon("download", class = "opt"), 
                                 up = FALSE))
                ))))
    })
    
    output$rna_exon_panels <- renderUI({
        fluidRow(column(12,
            tabBox(
                title="DEXSeq",
                width=12,
                height=8,
                tabPanel("Summary of DEXSeq",
                         fluidRow(
                             box(width=3,
                                 withSpinner(
                                     textOutput("exons_evaluated"),type=4,size=0.2, proxy.height="5px"
                                 ),
                                 p("Number of exons evaluated")
                             ),
                             box(width=3,
                                 withSpinner(
                                     textOutput("exons_de"),type=4,size=0.2, proxy.height="5px"),
                                 p("Number of exons differentially expressed")),
                             box(width=3,
                                 withSpinner(
                                     textOutput("exons_gene"),type=4,size=0.2, proxy.height="5px"),
                                 p("Number of genes evaluated")
                             ),
                             box(width=3,
                                 withSpinner(
                                     textOutput("exons_gene_de"),type=4,size=0.2, proxy.height="5px"),
                                 p("Number of genes with at least one exon differentially expressed")
                             )
                         ),
                         fluidRow(
                             box(width=6,
                                 selectizeInput("exon_gene_name", label="Gene", choices=c()),
                                 withSpinner(
                                     plotOutput("exons_per_gene"), type=4
                                 ),
                                 div(
                                     style = "position: absolute; left: 1em; bottom: 0em;",
                                     dropdown(
                                         selectInput(inputId="exons_gene_plot_down_ext", label="Download as...", choices=c("png","svg","pdf","eps","jpeg")),
                                         downloadButton(outputId = "exons_gene_plot_down_ext", label = "DOWNLOAD"),
                                         size = "xs",
                                         icon = icon("download", class = "opt"), 
                                         up = TRUE))
                                 
                             )
                         )
                ),
                tabPanel("DEXSeq table",
                         box(
                             withSpinner(
                                 DT::dataTableOutput("exons_table"), type=4
                             ),
                             div(
                                 style = "position: absolute; left: 1em; bottom: 0em;",
                                 dropdown(
                                     downloadButton(outputId = "exons_table_csv", label = "CSV"),
                                     downloadButton(outputId = "exons_table_excel", label = "EXCEL"),
                                     size = "xs",
                                     icon = icon("download", class = "opt"), 
                                     up = FALSE))
                             
                         )
                )
            )
        ))
    })
    
    output$rna_splice_panels <- renderUI({
        fluidRow(column(12,
            tabBox(
                title="Splicing events",
                width=12,
                tabPanel("Number of events",
                         column(width=6,
                             withSpinner(
                                 plotOutput("event_sum"),type=4
                             ),
                             div(
                                 style = "position: absolute; left: 1em; bottom: 0em;",
                                 dropdown(
                                     selectInput(inputId="majiq_event_sum_down_ext", label="Download as...", choices=c("png","svg","pdf","eps","jpeg")),
                                     downloadButton(outputId = "majiq_event_sum_down", label = "DOWNLOAD"),
                                     size = "xs",
                                     icon = icon("download", class = "opt"), 
                                     up = TRUE))
                             
                         ),
                         column(width=6,
                             withSpinner(
                                 plotOutput("majiq_heatmap"), type=4
                             ),
                             div(
                                 style = "position: absolute; left: 1em; bottom: 0em;",
                                 dropdown(
                                     selectInput(inputId="exons_hm_down_ext", label="Download as...", choices=c("png","svg","pdf","eps","jpeg")),
                                     downloadButton(outputId = "exons_hm_down", label = "DOWNLOAD"),
                                     size = "xs",
                                     icon = icon("download", class = "opt"), 
                                     up = TRUE
                                 ))
                         )),
                tabPanel("Event table",
                    withSpinner(DT::dataTableOutput("majiq_table"), type=4)
                ))))
    })
    
    output$com_splice_panels <- renderUI({
        fluidRow(column(12,
                        tabBox(title="Comparative analysis of splicing analysis tools",
                               width=12,
                               tabPanel("Overlap of genes",
                                        withSpinner(plotOutput("com_spli_venn"), type=4)))))
    })
    
    output$chip_qc_panels <- renderUI({
        fluidRow(column(12,
            tabBox(
                title="ChIP-Seq ENCODE",
                width=12,
                tabPanel("ENCODE table",
                         withSpinner(
                             DT::dataTableOutput("encode_table"), type=4
                         ),
                         div(
                             style = "position: absolute; left: 1em; bottom: 0em;",
                             dropdown(
                                 downloadButton(outputId = "encode_table_csv", label = "CSV"),
                                 downloadButton(outputId = "encode_table_excel", label = "EXCEL"),
                                 size = "xs",
                                 icon = icon("download", class = "opt"), 
                                 up = FALSE))
                         
                ),
                tabPanel("Peaks visualization",
                         withSpinner(
                             plotOutput("bedpeaks"), type=4
                         ),
                         div(
                             style = "position: absolute; left: 1em; bottom: 0em;",
                             dropdown(
                                 selectInput(inputId="chip_bedpk_down_ext", label="Download as...", choices=c("png","svg","pdf","eps","jpeg")),
                                 downloadButton(outputId = "chip_bedpk_down", label = "DOWNLOAD"),
                                 size = "xs",
                                 icon = icon("download", class = "opt"), 
                                 up = TRUE
                             ))
                         
                ),
                tabPanel("Gene track",
                         selectizeInput("gene_to_display", label="", choices=c()),
                         withSpinner(
                             plotOutput("chipgenetrack"), type=4
                         ),
                         div(
                             style = "position: absolute; left: 1em; bottom: 0em;",
                             dropdown(
                                 selectInput(inputId="chip_genetrack_down_ext", label="Download as...", choices=c("png","svg","pdf","eps","jpeg")),
                                 downloadButton(outputId = "chip_genetrack_down", label = "DOWNLOAD"),
                                 size = "xs",
                                 icon = icon("download", class = "opt"), 
                                 up = TRUE))
                )
            )
            
        ))
    })
    
    output$chip_peak_panels <- renderUI({
        fluidRow(column(12,
            tabBox(
                title="ChIP-seq peak enrichment",
                width=12,
                tabPanel("Gene level peak enrichment",
                         withSpinner(
                             plotOutput("enrichment_plot"),type=4
                         ),div(
                             style = "position: absolute; left: 1em; bottom: 0em;",
                             dropdown(
                                 selectInput(inputId="chip_geneenr_down_ext", label="Download as...", choices=c("png","svg","pdf","eps","jpeg")),
                                 downloadButton(outputId = "chip_geneenr_down", label = "DOWNLOAD"),
                                 size = "xs",
                                 icon = icon("download", class = "opt"), 
                                 up = TRUE
                             ))
                         ,
                         verbatimTextOutput("check")),
                tabPanel("Promoter level peak enrichment",
                         withSpinner(
                             plotOutput("pro_enrichment_plot"), type=4
                         ),
                         div(
                             style = "position: absolute; left: 1em; bottom: 0em;",
                             dropdown(
                                 selectInput(inputId="chip_pro_down_ext", label="Download as...", choices=c("png","svg","pdf","eps","jpeg")),
                                 downloadButton(outputId = "chip_pro_down", label = "DOWNLOAD"),
                                 size = "xs",
                                 icon = icon("download", class = "opt"), 
                                 up = TRUE))
                         
                )
            )
        ))
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
                # modalDialog(
                #     h3("Confirm using examples dataset."),
                #     p("If you are importing your own data, please click 'cancel' and uncheck 'Use example data'."),
                #     selectInput("rnaseq_metacol", "Example sample column:", selected = "sampleID", choices=c("sampleID")),
                #     selectInput("rnaseq_genecol", "Example gene column:", selected = "ensembl_gene_id_version", choices=c("ensembl_gene_id_version")),
                #     footer = tagList(
                #         modalButton("Cancel"),
                #         actionButton("ok_import_rna", "OK")
                #     )
                # )
                
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
    
    observeEvent(input$ok_import_rna,{
        removeModal()
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
            df <- make_gene_matrix("examples/gene_counts")
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
            df <- read.csv("examples/metadata.txt", sep=" ", header=T)
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
        if (input$useexamples==TRUE){
            genecol <- "ensembl_gene_id_version"
        } else {
            genecol <- input$rnaseq_genecol
        }
        seqdf <- rnaseq_df()
        row.names(seqdf) <- seqdf[,genecol]
        return(seqdf %>% dplyr::select(-genecol))
    })
    
    rnaseq_meta_clean <- reactive({
        metadf <- rnameta_df()
        if (input$useexamples==TRUE){
            metacol <- "sampleID"
        } else {
            metacol <- input$rnaseq_metacol
        }
        row.names(metadf) <- metadf[, metacol]
        return(metadf)
    })
    
    ## filter meta data 
    filtered_meta <- reactive({
        if (input$load_deseq2 == 0){return()}
        if (input$useexamples==TRUE){
            samplename <- "sampleID"
        } else {
            samplename <- input$rnaseq_metacol
        }
        metadf <- rnaseq_meta_clean()
        sub <- metadf %>% dplyr::select(samplename) 
        return(metadf[sub[,1] %in% input$rna_samples,])
        
    })

    ## generate dds object 
    #### TODO add covariates
    dds_obj <- reactive({
        if (input$useexamples == TRUE){
            load(file = "examples/gene_counts/dds.RData")
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
        if (input$useexamples==TRUE){
            samplenames <- "sampleID"
        } else {
            samplenames <- input$rnaseq_metacol
        }
        
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
                                  intgroup=c(input$rna_meta_var1), 
                                  returnData=TRUE)
        gene_counts
    })
    
    
    ## gene plot
    rna_genecount_plot <- eventReactive(input$gene_name, {
        req(input$rna_meta_var1)
        dds<-dds_obj()
        gene_counts <- plotCounts(dds, 
                                  gene=input$gene_name, 
                                  intgroup=c(input$rna_meta_var1), 
                                  returnData=TRUE)
        
        ggplot(data=gene_counts, aes(x = .data[[input$rna_meta_var1]], y = count, col=.data[[input$rna_meta_var1]])) +
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
        dds<-dds_obj()
        gene_counts <- plotCounts(dds, 
                                  gene=input$gene_name, 
                                  intgroup=c(input$rna_meta_var1), 
                                  returnData=TRUE)
        
        ggplot(data=gene_counts, aes(x = .data[[input$rna_meta_var1]], y = count, col=.data[[input$rna_meta_var1]])) +
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
    
    rna_volcano_plot <- eventReactive(input$load_deseq2, {
        
    })
    
    output$rna_volcano <- renderPlot({
        dds <- dds_obj()
        req(input$deseq_result)
        
        
        deg <- data.frame(results(dds, name=input$deseq_result))
        
        ggplot(data=deg,aes(x=log2FoldChange,y=-log10(padj))) +
            geom_vline(xintercept=c(-log2(as.numeric(input$log2cutoff)),log2(as.numeric(input$log2cutoff))), color="red")+ 
            geom_hline(yintercept=-log10(0.05), color="blue")+ 
            geom_point(color="black",alpha=0.5,stat="identity")+
            xlim(-10,10)
    })
        
    rna_deseq_res_table <- eventReactive(input$load_deseq2, {
        dds <- dds_obj()
        req(input$deseq_result)
        
        
        deg <- data.frame(DESeq2::results(dds, name=input$deseq_result))
        deg
    })
    
    output$rna_deseq_table <- DT::renderDataTable({rna_deseq_res_table()})
    
    sf_deseq2_plot <- reactive({
        sf_genes<- read.delim(file = "examples/splicing_factors_list.txt")
        
        if (input$useexamples==TRUE){
            degs<-read.table(file = "examples/gene_counts/degs_deseq2.txt",header = TRUE,sep = "\t")
            degs$ensembl_gene_id <- gsub(pattern = "\\..*$",replacement = "",x=degs$ensembl_gene_id_version)
            mart_export <-  mart_export_obj()
            
            degs <- dplyr::left_join(degs,mart_export,by="ensembl_gene_id")

        }
        
        heatmap_sf <- degs[degs$external_gene_name %in% sf_genes$Gene,]
        heatmap_sf$padj[is.na(heatmap_sf$padj)] <- 1
        ggplotly(ggplot(data=heatmap_sf,aes(x=log2FoldChange, y=external_gene_name, size = 1-padj,color=log2FoldChange, text=paste0("SF: ", external_gene_name))) +
            geom_point() +
            # xlim(-1,1)+
            # scale_size(range = c(1, 10), name="padj")+
            theme_classic()+
            geom_vline(xintercept=c(-log2(as.numeric(input$log2cutoff)),log2(as.numeric(input$log2cutoff))), color="black")+
            scale_color_gradient2(midpoint=0, low="#313695", mid="#FFFFBF",
                                  high="#A50026", space ="Lab" ))
            
    })
    
    output$sf_deseq2 <- renderPlotly({ sf_deseq2_plot() })
    
    ######################### isoform switch ###########################
    
    isoform_data <- reactive({
        if (input$useexamples == TRUE){
            load(file = "examples/pseudocounts/exampleSwitchListAnalyzed.RData")
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
    
    switch_plotplot <- reactive({
        if (input$load_iso==0){return()}
        exampleSwitchListAnalyzed <- isoform_data()
        isa_genes <- isoform_sig_genes()
        
        switchPlot(exampleSwitchListAnalyzed,
                   gene=input$iso_gene_name)
        
    })
    
    iso_switch_sumplot <- reactive({
        #if (input$load_iso==0){return()}
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
    
    switch_table_tb <- reactive({
        df <- isoform_sig_genes()
        return(df)
    })
    
    output$switch_table <- renderDataTable({switch_table_tb()})
    
    splice_enrich_plot <- reactive({
        exampleSwitchListAnalyzed <- isoform_data()
        extractSplicingEnrichment(
            exampleSwitchListAnalyzed,
            splicingToAnalyze='all',
            returnResult=TRUE,
            returnSummary=FALSE
        )
        
        
    })
    
    splice_genomewide <- reactive({
        exampleSwitchListAnalyzed <- isoform_data()
        extractSplicingGenomeWide(
            exampleSwitchListAnalyzed,
            featureToExtract = 'all',                 # all isoforms stored in the switchAnalyzeRlist
            splicingToAnalyze = c('A3','MES','ATSS'), # Splice types significantly enriched in COAD
            plot=TRUE,
            returnResult=FALSE  # Preventing the summary statistics to be returned as a data.frame
        )
    })
    
    output$splice_genom <- renderPlot({ splice_genomewide() })
    output$splice_enrich <- renderPlot({ splice_enrich_plot() })
    ######################### DEXSeq ################################
    
    exon_data <- reactive({
        if (input$useexamples==TRUE){
            load(file = "examples/exon_counts/dxr1.RData")
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
    
    exons_table_tb <- reactive({
        dxr1 <- exon_data()
        return(data.frame(dxr1) %>% na.omit())
    })
    
    output$exons_table <- renderDataTable({voila_res_all()})

   observe({
        if (is.null(exon_data())){return()}
        if (input$tabs=="rna_exon"){
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
        }
    })
   
   observe({
       sig_exons <- sig_exon_data()
       if (input$tabs=="rna_exon"){
        updateSelectizeInput(session, "exon_gene_name", label="Select a gene", choices=unique(sig_exons$groupID), server=TRUE)
       }
   })
   
    exons_per_gene_plot <- eventReactive(input$exon_gene_name, {
        if (is.null(exon_data())){return()}
        dxr1 <- exon_data()
        plotDEXSeq(dxr1, gene=input$exon_gene_name, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
    })
    
    output$exons_per_gene <- renderPlot({exons_per_gene_plot()})
    
    ########################## MAJIQ ####################################
    
    ################
    ## Exon Skipping
    ################
    
    voila_res_all <- reactive({
        if (input$useexamples==TRUE){
            df<-read.delim(file = "examples/majiq_output/voila_results_all.tsv",header = TRUE,sep = "\t")
        }
        return(df)
    })
    
    voila_res <- reactive({
        if(input$useexamples==TRUE){
            df <- read.delim(file = "examples/majiq_output/voila_results.tsv",header = TRUE,sep = "\t")
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
    
    output$majiq_table <- renderDataTable({ majiq.NrGenes.df() })
        
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
    
    com_spli_venn_plot <- reactive({
        ### from isa
        
        genes.iso <- read.delim(file = "examples/pseudocounts/isoform_signifcant_genes.txt")
        genes.iso <- unique(genes.iso$gene_id)
        
        ### from dexseq
        dxr1 <- exon_data()
        genes.exo <- as.data.frame(dxr1) %>%
            dplyr::filter(padj <  input$exons_pval_thres )
        genes.exo <- unique(genes.exo$groupID)
        
        ### from majiq
        genes.majiq <- voila_res()
        genes.majiq <- unique(genes.majiq$Gene.ID)
        
        venn.splicing <- venn.diagram(x = list(DEXseq=genes.exo,
                                               IsoformSwitchAnalyzer=genes.iso,
                                               Majiq=genes.majiq),
                                      filename = NULL,
                                      col="black",
                                      fill=RColorBrewer::brewer.pal(n = 3,name = "Dark2"))
        
        grid.newpage(); grid::grid.draw(venn.splicing)
    })
    
    output$com_spli_venn <- renderPlot({ com_spli_venn_plot() })
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
    

    observe({
        if (input$useexamples==TRUE) {
            updateSelectInput(session, "rna_meta_var1", "Select condition column", choices=colnames(rnameta_df()), selected=colnames(rnameta_df())[2])
            updateSelectizeInput(session, "rna_meta_var2", "Select a variable column (if any)", choices=colnames(rnameta_df()))
            updateCheckboxGroupInput(session, "rna_samples", "Select samples for analysis", choices=colnames(rnaseq_df())[!grepl("gene", colnames(rnaseq_df()))], selected=colnames(rnaseq_df()))
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
        updateSelectInput(session, "deseq_result", "Choose a comparison", choices=resultsNames(dds)[2:length(resultsNames(dds))])
        updateSelectizeInput(session, "gene_name", "Select gene", choices=row.names(counts(dds)), server = TRUE)
    })


    
    ##################### quality control #########################
    
    # star qc
    star_qc <- reactive({
        if (input$useexamples == TRUE){
            alignment_stats_star <- read.table(file = "examples/multiqc_data/multiqc_star.txt",header = TRUE,sep = "\t")
            alignment_stats_star <- melt(data = alignment_stats_star,id.vars="Sample")#,"total_reads","avg_input_read_length"))
        }
        return(alignment_stats_star)
    })
    
    star_qc_plot <- reactive({
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
    
    star_qc_readlength <- reactive({
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
    
    pseudo_qc <- reactive({
        if (input$useexamples == TRUE){
            alignment_stats_kallisto <- read.csv("examples/multiqc_data/multiqc_kallisto.txt",header = TRUE,sep = "\t")
            alignment_stats_kallisto <- melt(data = alignment_stats_kallisto,id.vars="Sample")#,"total_reads","avg_input_read_length"))
        }
        
        return(alignment_stats_kallisto)
    })
    
    pseudo_qc_plot <- reactive({
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
    
    #################### chip quality control ####################

    
    
    observeEvent(input$load_encode, {
        
            updateSelectInput(session, "fileassembly", "Select an assembly version", choices=unique(ENCODE_metadata$File.assembly))
    })
    
    observeEvent(input$fileassembly, {
        fas_meta <- ENCODE_metadata %>% dplyr::filter(File.assembly==input$fileassembly) 
        updateSelectInput(session, "biosample", "Select a sample", choices=unique(fas_meta$Biosample.term.name))
    })
    
    filter_encode_meta <- reactive({
        req(input$fileassembly)
        req(input$biosample)
        encode_meta <- ENCODE_metadata %>%
            dplyr::filter(File.assembly==input$fileassembly) %>%
            dplyr::filter(Output.type %in% c("optimal IDR thresholded peaks")) %>%
            dplyr::filter(Biosample.term.name==input$biosample) %>% # Maybe this can be also selected by the user
            dplyr::select(File.accession,Assay,Biosample.term.name,Experiment.target)
        encode_meta$Experiment.target <- gsub(pattern = "-human",replacement = "",x = encode_meta$Experiment.target)
        encode_meta <- encode_meta %>%
            mutate(Assay = if_else(Experiment.target %in% c(sf_list$Gene,"TARDBP","RBM25","RBFOX2","PTBP1","KHSRP","FUS","SF1","CELF1") | grepl(pattern = "HNRNP",x = Experiment.target),
                                   true = "Splicing Factor",
                                   false = Assay))
        input_chips <- encode_meta %>%
            dplyr::filter(Assay == "Splicing Factor")
        
        
        return(input_chips)
        
    })
    
    output$encode_table <- renderDataTable({
        filter_encode_meta()
    })
    
    chip_name <- reactive({
        if (input$useexamples_chip==TRUE){
            user_chip_name <- "YBX1"
        } else {
            user_chip_name <- input$chipseq_name
        }
        return(user_chip_name)
    })

    
    allchips_obj <- reactive({
        if (input$useexamples_chip==TRUE){
            user_chip_file  <- "examples/encode_bedNarrowPeak_files/ENCFF520DIY.bed.gz"
            
        } else {
            user_chip_file <- input$chipseq_bed
        }
        
        user_chip_name <- chip_name()
        input_chips <- filter_encode_meta()
        all.chips <- c()
        for (accession in 1:nrow(input_chips)) {
            all.chips[[input_chips$Experiment.target[accession]]] <- IRanges::reduce(readPeakFile(peakfile = paste0("examples/encode_bedNarrowPeak_files/",input_chips$File.accession[accession],".bed.gz")))
        }
        
        all.chips[[user_chip_name]] <- IRanges::reduce(readPeakFile(peakfile = user_chip_file),)
        
        for (chip in 1:length(all.chips)){
            all.chips[[chip]]@seqnames <- gsub(pattern = "chr",replacement = "",x = all.chips[[chip]]@seqnames)
        }
        return(all.chips)
        
    })
    
    txdb_gff <- reactive({
        if (input$gtf ==TRUE) {
            load("examples/peaksregion.RData")
            return(txdb_object)
        }
    })
    
    df_peakanno_edb <- reactive({
        if (input$useexamples_chip==TRUE){
            load("examples/peaksregion.RData")
            return(df_peakAnno.edb)
        } else {
            allchips <- allchips_obj()
            txdb_object <- txdb_gff()
            peakAnno.edb <- c()
            df_peakAnno.edb <- c()
            
            for (i in 1:length(all.chips)){
                peakAnno.edb[[names(all.chips)[i]]] <- assignChromosomeRegion(all.chips[[i]], nucleotideLevel=FALSE,
                                                                              TxDb=txdb_object)
                data_tmp <- as.data.frame(peakAnno.edb[[i]]$percentage) %>%  mutate(dataset=names(all.chips)[i])
                df_peakAnno.edb <- as.data.frame(rbind(df_peakAnno.edb,data_tmp))
            }
            
            df_peakAnno.edb$subjectHits <- gsub(pattern = "fiveUTRs",replacement = "5'UTR",x = df_peakAnno.edb$subjectHits)
            df_peakAnno.edb$subjectHits <- gsub(pattern = "threeUTRs",replacement = "3'UTR",x = df_peakAnno.edb$subjectHits)
            df_peakAnno.edb$subjectHits <- gsub(pattern = "immediateDownstream",replacement = "Immediate downstream",x = df_peakAnno.edb$subjectHits)
            df_peakAnno.edb$subjectHits <- gsub(pattern = "Intergenic.Region",replacement = "Intergenic",x = df_peakAnno.edb$subjectHits)
            return(df_peakAnno.edb)
        }
    })
    
    df_peakanno_plot <- reactive({
        df_peakAnno.edb <- df_peakanno_edb()
        df_peakAnno.edb %>%
            ggplot( aes(fill=subjectHits, y=Freq, x=dataset)) +
            geom_bar(position="fill", stat="identity",color="black") +
            coord_flip()
    })
    
    output$bedpeaks <- renderPlot({df_peakanno_plot()})
    
    get_mart <- reactive({
        mart = useMart("ensembl",dataset="hsapiens_gene_ensembl")
        fm = Gviz:::.getBMFeatureMap()
        fm["symbol"] = "external_gene_name"
        
        bm = BiomartGeneRegionTrack(biomart=mart,
                                    size=2, name="GENCODE annotation", utr5="red", utr3="red",
                                    protein_coding="black", col.line=NULL, #cex=7,
                                    collapseTranscripts="longest",
                                    featureMap=fm)
    })
    
    mart_export_obj <- reactive({
        mart_export <- read.table(file = "examples/mart_export.txt",sep = "\t",header = T)
        mart_export$chromosome_name <- paste0("chr",mart_export$chromosome_name)
        return(mart_export)
    })
    
    observe({
        if (input$tabs=="chip_qc"){
            mart_export <- mart_export_obj()
            updateSelectizeInput(session, "gene_to_display", label="Select a gene", choices=unique(mart_export$external_gene_name), selected="BTNL10", server=T)
        }
    })
    
    plot_gene_track <- reactive({
        all.chips <- allchips_obj()
        AT = GenomeAxisTrack()
        bm <- get_mart()
        user_chip_name <- chip_name()
        mart_export <- mart_export_obj() 
        
        TrackUserChip = AnnotationTrack(all.chips[[user_chip_name]], name=user_chip_name, shape='box',fill='blue',size=2,col="blue")
        
        plotTracks(c(TrackUserChip, bm, AT),
                   chromosome = mart_export$chromosome_name[mart_export$external_gene_name==input$gene_to_display],
                   from=as.numeric(mart_export$start_position[mart_export$external_gene_name==input$gene_to_display][1]),
                   to=as.numeric(mart_export$end_position[mart_export$external_gene_name==input$gene_to_display][1]),
                   extend.right = 3000,extend.left = 3000,
                   transcriptAnnotation="symbol",
                   window="auto",
                   cex.title=1, fontsize=8,
                   type="histogram" )
        
    })
    
    output$chipgenetrack <- renderPlot({plot_gene_track()})
    
    #################### peak enrichment ##########################
    make_granges_obj <- function(x, tool, type="gene"){
        ####
        astool_obj <- x
        
        mart_export <- mart_export_obj()
        if (tool=="majiq"){
            pos_genes_to_select <- unique(astool_obj$X.Gene.Name)
        } else if (tool=="dexseq_pos"){
            #### dexseq obj
            DEXseq_res <- as.data.frame(x)
            pos_genes_to_select <- gsub(pattern = "\\..*",replacement = "",x = unique(DEXseq_res$groupID[DEXseq_res$padj<0.1]))
        } else if (tool=="dexseq_neg"){
            DEXseq_res <- as.data.frame(x)
            pos_genes_to_select <- gsub(pattern = "\\..*",replacement = "",x = unique(DEXseq_res$groupID[DEXseq_res$padj>0.9]))
        } else if (tool=="isa_pos") {
            exampleSwitchListAnalyzed <- x
            isa_res_pos <- IsoformSwitchAnalyzeR::extractTopSwitches(exampleSwitchListAnalyzed,alpha = 0.05,dIFcutoff = 0.1,n=Inf)
            pos_genes_to_select <- isa_res_pos$gene_name
        } else if (tool=="isa_neg") {
            exampleSwitchListAnalyzed <- x
            isa_res_neg <- IsoformSwitchAnalyzeR::extractTopSwitches(exampleSwitchListAnalyzed,alpha = 1,dIFcutoff = 0,n=Inf)
            isa_res.temp <- IsoformSwitchAnalyzeR::extractTopSwitches(exampleSwitchListAnalyzed,alpha = 0.9,dIFcutoff = 0,n=Inf)
            isa_res_neg <- isa_res_neg %>%
                dplyr::filter(gene_name %notin% isa_res.temp$gene_name)
            rm(isa_res.temp)
            pos_genes_to_select <- isa_res_neg$gene_name
        }
        
        if (tool=="dexseq_pos" | tool=="dexseq_neg"){
            mart_filt <- mart_export %>%
                dplyr::filter(ensembl_gene_id %in% pos_genes_to_select)
        } else {
            mart_filt <- mart_export %>%
                dplyr::filter(external_gene_name %in% pos_genes_to_select)
        }
        mart_filt <- unique(mart_filt[c("ensembl_gene_id","external_gene_name","start_position","end_position","chromosome_name")])
        ####
        if (type=="gene"){
            this.gene <-   GRanges(seqnames = Rle( mart_filt$chromosome_name),
                                          ranges = IRanges( start = mart_filt$start_position,
                                                            end = mart_filt$end_position),
                                          strand = Rle( rep("*", nrow(mart_filt)) ),
                                          gene = mart_filt$external_gene_name)
            return(this.gene)
        } else if (type=="promoter"){
            this.promoter <- GRanges(seqnames = Rle( mart_filt$chromosome_name),
                                            ranges = IRanges( start = mart_filt$start_position -200,
                                                              end = mart_filt$start_position +200),
                                            strand = Rle( rep("*", nrow(mart_filt)) ),
                                            gene = mart_filt$external_gene_name)
            return(this.promoter)
        }
    }
    
    
    majiq.gr.es.gene.ob <- reactive({
        make_granges_obj(majiq.es.pos(), tool="majiq", type="gene")
    })
    
    majiq.gr.es.promoter.ob <- reactive({
        make_granges_obj(majiq.es.pos(), tool="majiq",type="promoter")
    })
    
    majiq.gr.es.gene.neg.ob <- reactive({
        make_granges_obj(majiq.es.neg(), tool="majiq",type="gene")
    })
    
    majiq.gr.es.promoter.neg.ob <- reactive({
        make_granges_obj(majiq.es.neg(), tool="majiq",type="promoter")
    })
    
    majiq.gr.ir.gene.ob <- reactive({
        make_granges_obj(majiq.ir.pos(), tool="majiq",type="gene")
    })
    
    majiq.gr.ir.promoter.ob <- reactive({
        make_granges_obj(majiq.ir.pos(), tool="majiq",type="promoter")
    })
    
    majiq.gr.ir.gene.neg.ob <- reactive({
        make_granges_obj(majiq.ir.neg(), tool="majiq",type="gene")
    })
    
    majiq.gr.ir.promoter.neg.ob <- reactive({
        make_granges_obj(majiq.ir.neg(), tool="majiq",type="promoter")
    })
    
    majiq.gr.a5.gene.ob <- reactive({
        make_granges_obj(majiq.a5.pos(), tool="majiq",type="gene")
    })
    
    majiq.gr.a5.promoter.ob <- reactive({
        make_granges_obj(majiq.a5.pos(), tool="majiq",type="promoter")
    })
    
    majiq.gr.a5.gene.neg.ob <- reactive({
        make_granges_obj(majiq.a5.neg(), tool="majiq",type="gene")
    })
    
    majiq.gr.a5.promoter.neg.ob <- reactive({
        make_granges_obj(majiq.a5.neg(), tool="majiq",type="promoter")
    })
    
    majiq.gr.a3.gene.ob <- reactive({
        make_granges_obj(majiq.a3.pos(), tool="majiq",type="gene")
    })
    
    majiq.gr.a3.promoter.ob <- reactive({
        make_granges_obj(majiq.a3.pos(), tool="majiq",type="promoter")
    })
    
    majiq.gr.a3.gene.neg.ob <- reactive({
        make_granges_obj(majiq.a3.neg(), tool="majiq",type="gene")
    })
    
    majiq.gr.a3.promoter.neg.ob <- reactive({
        make_granges_obj(majiq.a3.neg(), tool="majiq",type="promoter")
    })
    
    dexseq.gr.ex.gene.ob <- reactive({
        make_granges_obj(exon_data(), tool="dexseq_pos", type="gene")
    })
    
    dexseq.gr.ex.promoter.ob <- reactive({
        make_granges_obj(exon_data(), tool="dexseq_pos", type="promoter")
    })
    
    dexseq.gr.ex.gene.neg.ob <- reactive({
        make_granges_obj(exon_data(), tool="dexseq_neg", type="gene")
    })
    
    dexseq.gr.ex.promoter.neg.ob <- reactive({
        make_granges_obj(exon_data(), tool="dexseq_neg", type="promoter")
    })
    
    isa.gr.iso.gene.ob <- reactive({
        make_granges_obj(isoform_data(), tool="isa_pos", type="gene")
    })
    
    isa.gr.iso.promoter.ob <- reactive({
        make_granges_obj(isoform_data(), tool="isa_pos", type="promoter")
    })
    
    isa.gr.iso.gene.neg.ob <- reactive({
        make_granges_obj(isoform_data(), tool="isa_neg", type="gene")
    })
    
    isa.gr.iso.promoter.neg.ob <- reactive({
        make_granges_obj(isoform_data(), tool="isa_neg", type="promoter")
    })
    
    #################### get enrichment result ####################
    
    get_enrich_results <- function(pos, neg, name, enrichment_results) {
        all.chips <- allchips_obj()
        obsFreq_interest <- c()
        for (chip in 1:length(all.chips)){
            # freq_bg=length(subsetByOverlaps(gr.gencode.genes,all_chips[[chip]]))
            freq_as=length(subsetByOverlaps(pos, all.chips[[chip]]))
            freq_nas=length(subsetByOverlaps(neg, all.chips[[chip]]))
            obsFreq_interest <- as.data.frame(rbind(obsFreq_interest,
                                                    c(as.numeric(freq_as),as.numeric(length(pos))-as.numeric(freq_as),
                                                      as.numeric(freq_nas),as.numeric(length(neg))-as.numeric(freq_nas),
                                                      # as.numeric(freq_bg),as.numeric(length(gr.gencode.genes))-as.numeric(freq_bg),
                                                      names(all.chips)[chip])))
        }
        colnames(obsFreq_interest) <- c("freqASwithPeak","freqASnoPeak",
                                        "freqNASwithPeak","freqNASnoPeak",
                                        # "freqBGwithPeak","freqBGnoPeak",
                                        "chip_protein")
        #3.3. Generate the contingency table
        library(reshape)
        final_contingency_table <- c()
        contingency_table <- melt(data = obsFreq_interest,id.vars = "chip_protein")
        contingency_table <- contingency_table %>%
            dplyr::mutate(event=if_else(condition = grepl(pattern = "NAS",x = contingency_table$variable),true = "NAS","AS"))  %>%
            dplyr::mutate(peak=if_else(condition = grepl(pattern = "withPeak",x = contingency_table$variable),true = "Peak","noPeak"))
        contingency_table$variable <- NULL
        for (chip in names(all.chips)){
            df <- reshape::cast(contingency_table,formula = peak ~event,value.var = "value",fun.aggregate = NULL,subset = chip_protein== chip)
            df <- df %>% 
                dplyr::mutate(chip_protein = chip)
            df <- df[c(2,1),]
            final_contingency_table <- as.data.frame(rbind(final_contingency_table,
                                                           df))
        }
        
        #3.4. Perform statistical test (Chi square)
        for (chip in names(all.chips)){
            chip_df <- final_contingency_table %>%
                dplyr::filter(chip_protein == chip)%>%
                dplyr::select(AS,NAS)
            chip_df <- apply(X = chip_df,MARGIN = 2,FUN = as.numeric)
            rownames(chip_df) <- c("noPeak","Peak")
            chisq <- fisher.test(chip_df,alternative = "two.sided")
            pval <- chisq$p.value
            estimate <- chisq$estimate
            enrichment_results <- as.data.frame(rbind(enrichment_results,
                                                      c(pval,estimate,chip, name)))
        }
        
        rm(df)
        return(enrichment_results)
    }
    
    get_gene_enrichment <- reactive({
        enrichment_results <- c()
        enrichment_results <- get_enrich_results(majiq.gr.ir.gene.ob(), majiq.gr.ir.gene.neg.ob(), name="Majiq - Intron Retention", enrichment_results = enrichment_results)
        cat(file=stderr(),"ok1")
        enrichment_results <- get_enrich_results(majiq.gr.es.gene.ob(), majiq.gr.es.gene.neg.ob(), name="Majiq - Exon Skipping", enrichment_results = enrichment_results)
        cat(file=stderr(),"ok2")
        enrichment_results <- get_enrich_results(majiq.gr.a5.gene.ob(), majiq.gr.a5.gene.neg.ob(), name="Majiq - A5SS" , enrichment_results = enrichment_results)  
        cat(file=stderr(),"ok3")
        enrichment_results <- get_enrich_results(majiq.gr.a3.gene.ob(), majiq.gr.a3.gene.neg.ob(), name="Majiq - A3SS" , enrichment_results = enrichment_results) 
        cat(file=stderr(),"ok4")
        enrichment_results <- get_enrich_results(dexseq.gr.ex.gene.ob(), dexseq.gr.ex.gene.neg.ob(), name="DEXseq - Differential exon" , enrichment_results = enrichment_results) 
        cat(file=stderr(),"ok5")
        enrichment_results <- get_enrich_results(isa.gr.iso.gene.ob(), isa.gr.iso.gene.neg.ob(), name="IsoformSwitchAnalyzer - Isoform switch" , enrichment_results = enrichment_results) 
        colnames(enrichment_results) <- c("pvalue","estimate","chip_dataset","Event type - Tool")
        return(enrichment_results)
    })
    
    get_promoter_enrichment <- reactive({
        enrichment_results <- c()
        enrichment_results <- get_enrich_results(majiq.gr.ir.promoter.ob(), majiq.gr.ir.promoter.neg.ob(), name="Majiq - Intron Retention", enrichment_results = enrichment_results)
        cat(file=stderr(),"ok1")
        enrichment_results <- get_enrich_results(majiq.gr.es.promoter.ob(), majiq.gr.es.promoter.neg.ob(), name="Majiq - Exon Skipping", enrichment_results = enrichment_results)
        cat(file=stderr(),"ok2")
        enrichment_results <- get_enrich_results(majiq.gr.a5.promoter.ob(), majiq.gr.a5.promoter.neg.ob(), name="Majiq - A5SS" , enrichment_results = enrichment_results)  
        cat(file=stderr(),"ok3")
        enrichment_results <- get_enrich_results(majiq.gr.a3.promoter.ob(), majiq.gr.a3.promoter.neg.ob(), name="Majiq - A3SS" , enrichment_results = enrichment_results) 
        cat(file=stderr(),"ok4")
        enrichment_results <- get_enrich_results(dexseq.gr.ex.promoter.ob(), dexseq.gr.ex.promoter.neg.ob(), name="DEXseq - Differential exon" , enrichment_results = enrichment_results) 
        cat(file=stderr(),"ok5")
        enrichment_results <- get_enrich_results(isa.gr.iso.promoter.ob(), isa.gr.iso.promoter.neg.ob(), name="IsoformSwitchAnalyzer - Isoform switch" , enrichment_results = enrichment_results) 
        colnames(enrichment_results) <- c("pvalue","estimate","chip_dataset","Event type - Tool")
        return(enrichment_results)
    })
    
    enrichment_plot_ob <- eventReactive(input$run_enr, {
        ht_data <- get_gene_enrichment()
        ht_data$pvalue[as.numeric(ht_data$pvalue) > 0.05] <- NA
        ht_data$pvalue[as.numeric(ht_data$estimate) == "Inf"] <- NA
        ggplot(ht_data, aes(x = chip_dataset, y = `Event type - Tool`, fill = as.numeric(pvalue))) + 
            geom_tile() + 
            labs(fill = "Chi Square pvalue") + 
            ggtitle("Enrichment analysis of ChIP-seq\npeak in spliced genes")+
            geom_text(aes(label=round(as.numeric(estimate),digits = 2))) +
            scale_fill_gradientn(colors= rev(RColorBrewer::brewer.pal(n = 6,name = "BuPu")),na.value = "white")+
            theme(axis.text.x=element_text(angle = 45))
        
    })
    
    pro_enrichment_plot_ob <- eventReactive(input$run_enr, {
        ht_data <- get_promoter_enrichment()
        ht_data$pvalue[as.numeric(ht_data$pvalue) > 0.05] <- NA
        ht_data$pvalue[as.numeric(ht_data$estimate) == "Inf"] <- NA
        ggplot(ht_data, aes(x = chip_dataset, y = `Event type - Tool`, fill = as.numeric(pvalue))) + 
            geom_tile() + 
            labs(fill = "Chi Square pvalue") + 
            ggtitle("Enrichment analysis of ChIP-seq\npeak in spliced genes")+
            geom_text(aes(label=round(as.numeric(estimate),digits = 2))) +
            scale_fill_gradientn(colors= rev(RColorBrewer::brewer.pal(n = 6,name = "BuPu")),na.value = "white")+
            theme(axis.text.x=element_text(angle = 45))
    })
    
    output$pro_enrichment_plot <- renderPlot({ pro_enrichment_plot_ob() })
    output$enrichment_plot <- renderPlot({ enrichment_plot_ob() })
    output$show_fileassembly <- renderText({ input$fileassembly })
    output$show_biosample <- renderText({ input$biosample })
    
    #################### download handler #########################
    output$deseq2_pca_down <- downloadHandler(
        filename = function() {paste0("dasire_", "deseq2_pca", ".", input$deseq2_pca_down_ext,  sep="")},
        content = function(file) {ggsave(file, plot = rnapcaplot(), device = input$deseq2_pca_down_ext, width = 10)})

    output$deseq2_hm_down <- downloadHandler(
        filename = function() {paste0("dasire_", "deseq2_heatmap", ".", input$deseq2_hm_down_ext,  sep="")},
        content = function(file) {ggsave(file, plot = rnaheatmap(), device = input$deseq2_hm_down_ext, width = 10)})
    
    output$sf_deseq2_down <- downloadHandler(
        filename = function() {paste0("dasire_", "sf_deseq2", ".", input$sf_deseq2_down_ext,  sep="")},
        content = function(file) {ggsave(file, plot = sf_deseq2_plot(), device = input$sf_deseq2_down_ext, width = 10)})
    
    output$deseq2_gc_down <- downloadHandler(
        filename = function() {paste0("dasire_", "rna_gc", ".", input$deseq2_gc_down_ext,  sep="")},
        content = function(file) {ggsave(file, plot = rna_genecount_plot(), device = input$deseq2_gc_down_ext, width = 10)})

    output$deseq2_vc_down <- downloadHandler(
        filename = function() {paste0("dasire_", "rna_volcano", ".", input$deseq2_vc_down_ext,  sep="")},
        content = function(file) {ggsave(file, plot = rna_volcano_plot(), device = input$deseq2_vc_down_ext, width = 10)})

    output$star_pp_down <- downloadHandler(
        filename = function() {paste0("dasire_", "star_pp", ".", input$star_pp_down_ext,  sep="")},
        content = function(file) {ggsave(file, plot = star_qc_plot(), device = input$star_pp_down_ext, width = 10)})

    output$star_rp_down <- downloadHandler(
        filename = function() {paste0("dasire_", "star_rp", ".", input$star_rp_down_ext,  sep="")},
        content = function(file) {ggsave(file, plot = star_qc_readlength(), device = input$star_rp_down_ext, width = 10)})

    output$kall_pp_down <- downloadHandler(
        filename = function() {paste0("dasire_", "kallisto_pp", ".", input$kall_pp_down_ext,  sep="")},
        content = function(file) {ggsave(file, plot =pseudo_qc_plot(), device = input$kall_pp_down_ext, width = 10)})
    
    output$iso_switch_sum_down <- downloadHandler(
        filename = function() {paste0("dasire_", "iso_switch_sum", ".", input$iso_switch_sum_down_ext,  sep="")},
        content = function(file) {ggsave(file, plot =iso_switch_sumplot(), device = input$iso_switch_sum_down_ext, width = 10)})
    
    output$iso_switch_plot_down <- downloadHandler(
        filename = function() {paste0("dasire_", "iso_switch_gene_", input$iso_gene_name, ".", input$iso_switch_plot_down_ext,  sep="")},
        content = function(file) {ggsave(file, plot=switch_plotplot(), device = input$iso_switch_plot_down_ext, width = 10)})
    
    output$exons_gene_plot_down <- downloadHandler(
        filename = function() {paste0("dasire_", "dexseq_gene",input$exon_gene_name, ".", input$exons_gene_plot_down_ext,  sep="")},
        content = function(file) {ggsave(file, plot=exons_per_gene_plot(), device = input$exons_gene_plot_down_ext, width = 10)})
    
    output$majiq_event_sum_down <- downloadHandler(
        filename = function() {paste0("dasire_", "majiq_sum", ".", input$majiq_event_sum_down_ext,  sep="")},
        content = function(file) {ggsave(file, plot=event_sum_plot(), device = input$majiq_event_sum_down_ext, width = 10)})
    
    output$exons_hm_down <- downloadHandler(
        filename = function() {paste0("dasire_", "majiq_hm", ".", input$exons_hm_down_ext,  sep="")},
        content = function(file) {ggsave(file, plot=majiq_heatmap_plot(), device = input$exons_hm_down_ext, width = 10)})
    
    output$chip_geneenr_down <-  downloadHandler(
        filename = function() {paste0("dasire_", "chip_gene_en", ".", input$chip_geneenr_down_ext,  sep="")},
        content = function(file) {ggsave(file, plot=enrichment_plot_ob(), device = input$chip_geneenr_down_ext, width = 10)})
    
    output$chip_pro_down <- downloadHandler(
        filename = function() {paste0("dasire_", "chip_pro_en", ".", input$chip_pro_down_ext,  sep="")},
        content = function(file) {ggsave(file, plot=pro_enrichment_plot_ob(), device = input$chip_pro_down_ext, width = 10)})
    
    output$chip_bedpk_down <- downloadHandler(
        filename = function() {paste0("dasire_", "chip_bed_peak", ".", input$chip_bedpk_down_ext,  sep="")},
        content = function(file) {ggsave(file, plot=df_peakanno_plot(), device = input$chip_bedpk_down_ext, width = 10)})
    
    output$chip_genetrack_down <- downloadHandler(
        filename = function() {paste0("dasire_", "chip_gene_track", ".", input$chip_genetrack_down_ext,  sep="")},
        content = function(file) {ggsave(file, plot=plot_gene_track(), device = input$chip_genetrack_down_ext, width = 10)})
    
    output$deseq2_table_csv <- downloadHandler(
        filename=function() {paste0("dasire_","deseq2_table.csv")},
        content=function(file) {write.table(x=rna_deseq_res_table(), file=file, sep = ",")}
    )
    
    output$deseq2_table_excel <- downloadHandler(
        filename = function() {paste0("dasire_","deseq2_table.xlsx")},
        content=function(file) {write.xlsx(rna_deseq_res_table(), file=file)}
    )
    
    output$iso_table_csv <- downloadHandler(
        filename=function() {paste0("dasire_","iso_table.csv")},
        content=function(file) {write.table(x=switch_table_tb(), file=file, sep = ",")}
    )
    
    output$iso_table_excel <- downloadHandler(
        filename=function() {paste0("dasire_","iso_table.xlsx")},
        content=function(file) {write.xlsx(x=switch_table_tb(), file=file)}
    )
    
    output$exons_table_csv <- downloadHandler(
        filename=function() {paste0("dasire_","exons_table.csv")},
        content=function(file) {write.table(x=exons_table_tb(), file=file, sep = ",")}
    )
    
    output$exons_table_excel <- downloadHandler(
        filename=function() {paste0("dasire_","exons_table.xlsx")},
        content=function(file) {write.xlsx(x=exons_table_tb(), file=file)}
    )
    
    output$encode_table_csv <- downloadHandler(
        filename=function() {paste0("dasire_","encode_table.csv")},
        content=function(file) {write.table(x=filter_encode_meta(), file=file, sep = ",")}
    )
    
    output$encode_table_excel <- downloadHandler(
        filename=function() {paste0("dasire_","encode_table.xlsx")},
        content=function(file) {write.xlsx(x=filter_encode_meta(), file=file)}
    )
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

shinyApp(ui=ui, server=server)