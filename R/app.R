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
library(shinyBS)


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
                              h4("Select parameters:"),
                              selectInput("rna_meta_var1", "Select condition (choose 'time')", choices=c()),
                              selectizeInput("rna_meta_var2", "Select a variable", choices=c()),
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
            column(
                withSpinner(
                    plotlyOutput("rna_pca"),
                    type=4
                ),
                div(
                    selectInput("deseq_result", label = "Result", choices = c()),
                    withSpinner(
                        plotOutput("rna_volcano"),
                        type=4
                    )
                ),
                width=4
            ),
            column(
                div(
                    selectizeInput("gene_name", label = "Gene", choices = c()),
                    withSpinner(
                        plotOutput("rna_genecount"),
                        type=4
                    )
                ),
                width=4
            ), 
            column(
                withSpinner(
                    plotOutput("rna_heatmap"),
                    type=4
                ),
                width=4
            )
            
        )
        # fluidRow(
        #     id="deseq2_panel",
        #     box(
        #         verbatimTextOutput("check_meta")
        #     )
        # )
    })
    #################### data handler / modal ##########################
    import_data_modal <- function(failed=F){
        if (failed==F){
            if (input$useexamples == TRUE){
                modalDialog(
                    h3("Confirm using examples dataset."),
                    p("If you are importing your own data, please click 'cancel' and uncheck 'Use example data'."),
                    selectInput("rnaseq_metacol", "Automatic sample column:", selected = "Run", choices=c()),
                    selectInput("rnaseq_genecol", "Automatic sample column:", selected = "gene", choices=c()),
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
            df <- read.csv("../examples/sorted_pos_tpm.csv", sep="\t", header=T)
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
 
        if (input$useexamples == TRUE){
            df <- read.csv("../examples/meta.txt", sep="\t", header=T)
        } else {
            if (is.null(metafile)){
                showModal(import_data_modal(failed = TRUE))
                return(NULL)
            } 
            df <- read.csv(metafile$datapath, sep=input$rnaseq_meta_sep, header=input$header_rnameta)
        }
        removeModal()
        return(df)
    })
    
    rnaseq_df_clean <- reactive({
        seqdf <- rnaseq_df()
        req(input$rnaseq_genecol)
        row.names(seqdf) <- seqdf[,input$rnaseq_genecol]
        return(seqdf %>% dplyr::select(-input$rnaseq_genecol))
    })
    
    rnaseq_meta_clean <- reactive({
        req(input$rnaseq_metacol)
        metadf <- rnameta_df()
        row.names(metadf) <- metadf[,input$rnaseq_metacol]
        return(metadf)
    })
    
    ## filter meta data 
    filtered_meta <- reactive({
        if (input$run_pca_rna == 0){return()}
        metadf <- rnaseq_meta_clean()
        samplename <- input$rnaseq_metacol
        sub <- metadf %>% dplyr::select(samplename) 
        return(metadf[sub[,1] %in% input$rna_samples,])
        
    })
    
    ## generate dds object 
    #### TODO add covariates
    dds_obj <- reactive({
        if (input$run_pca_rna == 0){return()}
        seqdf <- rnaseq_df_clean()
        metadf <- filtered_meta()
        
        dds <- DESeqDataSetFromMatrix(countData=as.matrix(seqdf[,input$rna_samples]),
                                      colData=metadf[input$rna_samples,],
                                      design=formula(c("~", input$rna_meta_var1)))
        
        dds <- DESeq(dds)
        dds 
    })
    
    # dds_result <- reactive({
    #     if (input$ok_deseq_level){
    #         
    #         dds <- dds_obj()
    #         dds[,input$rna_meta_var1] <- relevel(dds[,input$rna_meta_var1], input$rna_deseq_level)
    #         dds <- DESeq(dds)
    #         return(dds)
    #     }
    # })
    
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
        
        seqdf <- rnaseq_df_clean()
        metadf <- filtered_meta()
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
    
    ## heatmap rnaseq
    rnaheatmap <- reactive({
        if (input$run_pca_rna==0){return()}
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
    
    check_gene <- reactive({
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
    
    rna_volcano_plot <- reactive({
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
        
    
    
    ######################### obersever #################################
    
    observeEvent(input$renderimport_rna, {
        if (is.null(rnaseq_df()) | is.null(rnameta_df())){return()}
        updateSelectInput(session, "rnaseq_genecol", "Select the gene column", choices=colnames(rnaseq_df()))
        updateSelectInput(session, "rnaseq_metacol", "Select the sample names column", selected = "Run", choices=colnames(rnameta_df()))
        
    })
    
    observeEvent(input$ok_import_rna,{
        removeModal()
        if(is.null(rnaseq_df())){return()}
        
        updateSelectInput(session, "rna_meta_var1", "Select group column", choices=colnames(rnameta_df()))
        updateSelectizeInput(session, "rna_meta_var2", "Select group column", choices=colnames(rnameta_df()))
        updateCheckboxGroupInput(session, "rna_samples", "Select samples for analysis", choices=colnames(rnaseq_df()), selected=colnames(rnaseq_df()))
    })
    
    ## gene selection
    observeEvent(input$rnaseq_genecol,{
        req(input$rnaseq_genecol)
        seqdf <- rnaseq_df()
        updateSelectizeInput(session, "gene_name", "Select gene", choices=seqdf[,input$rnaseq_genecol], server = TRUE)
    })
    
    observeEvent(input$run_pca_rna, {
        dds <- dds_obj()
        updateSelectInput(session, "deseq_result", "Choose a comparison", choices=resultsNames(dds))
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