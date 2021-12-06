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
        ),
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
            fluidRow(
                id="deseq2_panel",
                uiOutput("deseq2_panels")
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
                              tabName="deseq2",
                              h3("Select parameters:"),
                              selectInput("rna_meta_var1", "Select condition (choose 'time')", choices=c()),
                              selectInput("rna_meta_var2", "Select a variable", choices=c()),
                              checkboxGroupInput("rna_samples", "Select samples for analysis.", choices = c()),
                              bsButton("run_pca_rna", "Run DESeq2", icon=icon('chevron-right'))
                              ),
                     menuItem("Differential exon"),
                     menuItem("Isoform switch"),
                     menuItem("Splicing event")
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
            id="deseq2_panel",
            box(
                withSpinner(
                    plotlyOutput("rna_pca"),
                    type=4
                ),
                width=6
            ), 
            box(
                withSpinner(
                    plotOutput("rna_heatmap"),
                    type=4
                ),
                width=6
            )
        )
        # fluidRow(
        #     id="deseq2_panel",
        #     box(
        #         verbatimTextOutput("check_meta")
        #     )
        # )
    })
    #################### data handler ##########################
    import_data_modal <- function(failed=F){
        if (failed==F){
            if (input$useexamples == TRUE){
                modalDialog(
                    h3("Confirm using examples dataset."),
                    p("If you are importing your own data, please click 'cancel' and uncheck 'Use example data'."),
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
    
    #################### reactive object #####################3
    ##show import dialog
    rnaseq_df <- reactive({
        rnafile <- input$rnaseq_csv
        if (is.null(rnafile)) {
            showModal(import_data_modal(failed = TRUE))
            return(NULL)
        }
        df <- read.csv(rnafile$datapath, sep=input$rnaseq_sep, header=input$header_rna)
        return(df)
    })
    
    # import meta data
    rnameta_df <- reactive({
        metafile <- input$rnaseq_meta
        
        if (is.null(metafile)){
            showModal(import_data_modal(failed = TRUE))
            return(NULL)
        } 
        df <- read.csv(metafile$datapath, sep=input$rnaseq_meta_sep, header=input$header_rnameta)
        removeModal()
        return(df)
    })
    
    ## filter meta data 
    filtered_meta <- reactive({
        if (input$run_pca_rna == 0){return()}
        metadf <- rnameta_df()
        samplename <- input$rnaseq_metacol
        sub <- metadf %>% dplyr::select(samplename) 
        return(metadf[sub[,1] %in% input$rna_samples,])
        
        })
    
    ## generate dds object 
    #### TODO add covariates
    dds_obj <- reactive({
        if (input$run_pca_rna == 0){return()}
        seqdf <- rnaseq_df()
        metadf <- filtered_meta()
        row.names(metadf) <- metadf$`input$rna_meta_samples`
        
        dds <- DESeqDataSetFromMatrix(countData=as.matrix(seqdf[,input$rna_samples]),
                                      colData=metadf,
                                      design=formula(c("~", input$rna_meta_var1)))

        dds
    })
    
    #####output
    check_meta1 <- reactive({
        if (input$run_pca_rna == 0){return()}
        filtered_meta()
    })
    
    #output
    output$check_meta <- renderPrint({input})
    
    ### make reactive pca plot
    rnapcaplot <- reactive({
        if (input$run_pca_rna == 0){return()}
        
        seqdf <- rnaseq_df()
        metadf <- rnameta_df()
        samplenames <- input$rnaseq_metacol
        
        pr <- prcomp(t(seqdf[,input$rna_samples]))
        pca <- data.frame(pr$x)
        pca$samples <- row.names(pca)
        allpca <- inner_join(pca, metadf, by=c("samples"=samplenames))
        
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
    
    ## 
    rnaheatmap <- reactive({
        if (input$run_pca_rna==0){return()}
        dds <- dds_obj()
        vst <- vst(dds, blind=FALSE)
        sampleDists <- dist(t(assay(vst)))
        sampleDistsMatrix <- as.matrix(sampleDists)
        
        colors <- colorRampPalette( rev(brewer.pal(length(dds$sample), "Blues")) )(255)
        
        
        pheatmap(sampleDistsMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors,legend = FALSE,main = "Sample distances",treeheight_row = 10,treeheight_col = 20)
    })
    
    output$rna_heatmap <- renderPlot({rnaheatmap()})
    ######################### obersever #################################
    observeEvent(input$renderimport_rna, {
        if (is.null(rnaseq_df()) | is.null(rnameta_df())){return()}
        updateSelectInput(session, "rnaseq_genecol", "Select the gene column", choices=colnames(rnaseq_df()))
        updateSelectInput(session, "rnaseq_metacol", "Select the sample names column", choices=colnames(rnameta_df()))
        
    })
    
    observeEvent(input$ok_import_rna,{
        removeModal()
        if(is.null(rnaseq_df())){return()}
        
        updateSelectInput(session, "rna_meta_var1", "Select group column", choices=colnames(rnameta_df()))
        updateSelectInput(session, "rna_meta_var2", "Select group column", choices=colnames(rnameta_df()))
        updateCheckboxGroupInput(session, "rna_samples", "Select samples for analysis", choices=colnames(rnaseq_df()), selected=colnames(rnaseq_df()))
    })
    

    output$check <- renderPrint({input$run_pca_rna})
    output$deseq2_panel <- renderUI({
        input$rna_pca_rna
        fluidRow(
            id="deseq2_panel",
            box(
                withSpinner(
                    plotlyOutput("rna_pca"),
                    type=4
                ),
                width=6
            ), 
            box(
                withSpinner(
                    plotOutput("rna_heatmap"),
                    type=4
                ),
                width=6
            )
        )
        # fluidRow(
        #     id="deseq2_panel",
        #     box(
        #         verbatimTextOutput("check_meta")
        #     )
        # )
    })
    
    output$rna_pca <- renderPlotly({
        rnapcaplot()
    })
        
    
    #################### dynamic rendering ##########################
    observeEvent("",{
        show("start_panel")
        hide("import_panel")
        hide("rnaseq")
        hide("deseq2_panel")
    }, once=T)

    observeEvent(input$tabs,{
        if (input$tabs=="import"){
            show("import_panel")
            hide("start_panel")
            hide("deseq2_panel")
        } else if (input$tabs=="start"){
            show("start_panel")
            hide("import_panel")
            hide("deseq2_panel")
        } 
    })
    
    observeEvent(input$run_pca_rna,{
        show("deseq2_panel")
        hide("start_panel")
        hide("import_panel")
    })
}

shinyApp(ui, server)