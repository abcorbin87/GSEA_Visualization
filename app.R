library(readr)
library(dplyr)
library(ggplot2)
library(magrittr)
library(SummarizedExperiment)
library(DESeq2)
library(airway)
library(enrichR)
library(matrixStats)
library(plotly)
library(data.table)
library(DT)
library(shinythemes)

# # functions and needed objects
# BiocManager::install("fgsea", dependencies = TRUE)
# library(fgsea)

GSEA <- function(gene_list, GO_file, pval) {
  set.seed(54321)
  
  library(dplyr)
  library(fgsea)
  
  if (any(duplicated(names(gene_list)))) {
    warning("Duplicates in gene names")
    gene_list <- gene_list[!duplicated(names(gene_list))]
  }
  if (!all(order(gene_list, decreasing = TRUE) == 1:length(gene_list))) {
    warning("Gene list not sorted")
    gene_list <- sort(gene_list, decreasing = TRUE)
  }
  myGO <- fgsea::gmtPathways(GO_file)
  
  fgRes <- fgsea::fgsea(
    pathways = myGO,
    stats = gene_list,
    minSize = 15, ## minimum gene set size
    maxSize = 400, ## maximum gene set size
    nperm = 10000
  ) %>%
    as.data.frame() %>%
    dplyr::filter(padj < !!pval) %>%
    arrange(desc(NES))
  message(paste("Number of significant gene sets =", nrow(fgRes)))
  
  message("Collapsing Pathways -----")
  concise_pathways <- collapsePathways(
    data.table::as.data.table(fgRes),
    pathways = myGO,
    stats = gene_list
  )
  fgRes <- fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  fgRes$Enrichment <- ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes <- rbind(head(fgRes, n = 10),
                   tail(fgRes, n = 10))
  
  total_up <- sum(fgRes$Enrichment == "Up-regulated")
  total_down <- sum(fgRes$Enrichment == "Down-regulated")
  header <- paste0("Top 10 (Total pathways: Up=", total_up, ", Down=", total_down, ")")
  
  colos <- setNames(c("firebrick2", "dodgerblue2"),
                    c("Up-regulated", "Down-regulated"))
  
  rg <- as.numeric(range(fgRes$NES) +
                     as.numeric(range(fgRes$NES) * 0.2))
  
  g1 <- ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point(aes(fill = Enrichment, size = size), shape = 21) +
    scale_fill_manual(values = colos) +
    scale_size_continuous(range = c(2, 10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    scale_x_discrete(expand = expansion()) +
    labs(x = "Pathway", y = "Normalized Enrichment Score",
         title = header) +
    theme() +
    theme(axis.text.y = element_text(size = 6))  # adjust the axis text size
  
  output <- list("Results" = fgRes, "Plot" = g1)
  return(output)
}

# Define UI for application
ui <- fluidPage(
  theme = shinytheme("darkly"),
  title = "Useful RNAseq Visualizations",
  tabsetPanel(
    tabPanel(
      titlePanel("GSEA"),
      # Sidebar
      sidebarLayout(
        sidebarPanel(
          p("Input .gmt file including your desired gene sets."),
          br(),
          p("Commonly used gene sets can be found here: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp"),
          br(),
          fileInput(
            'go_file',
            "Go_file",
            buttonLabel = "Browse",
            placeholder = ".symbol.gmt file not selected"
          ),
          fileInput(
            'gene_list',
            "Gene Lets",
            buttonLabel = "Browse",
            placeholder = ".csv file not selected"
          ),
          sliderInput(
            "pval",
            dragRange = FALSE,
            label = "P value cutoff",
            min = 0.001,
            max = 0.06,
            value = 0.05,
            step = 0.001
          ),
          actionButton(
            "go_button",
            "Generate GSEA Plot",
            class = "btn-info",
            style = "color: #fff; background-color: #bbbbbb; border-color: #2e6da4"
          ), width = 5
        ),
        # Show a plot
        mainPanel(
          plotlyOutput(outputId = "gsea_plot"),
          uiOutput("message"), width = 7
        )
      )
    )
  )
)

# Define server logic required to print plot
options(shiny.maxRequestSize = 30 * 1024^2)
server <- function(input, output, session) {
  
  colos <- setNames(c("firebrick2", "dodgerblue2"),
                    c("Up-regulated", "Down-regulated"))
  
  output$message <- renderUI({
    HTML('<div class="shinyjs-message" style="display: none;">Generating GSEA Results...</div>')
  })
  
  observeEvent(input$go_button, {
    shinyjs::enable("message")
    withProgress(
      message = "Generating GSEA Results...",
      detail = "This may take a moment...",
      value = 0, {
        
        pval <- reactive(input$pval)
        S4table <- read.csv(file = input$gene_list$datapath, header = TRUE, skip = 1) %>%
          filter(Gene.Symbol != "")
        gene_list <- S4table$DESeq2.Log2.Fold.Change
        names(gene_list) <- S4table$Gene.Symbol
        gene_list <- sort(gene_list, decreasing = TRUE)
        gene_list <- gene_list[!duplicated(names(gene_list))]
        
        # perform GSEA
        res <- GSEA(gene_list = gene_list, GO_file = input$go_file$datapath, pval = pval())
        
        # total_up <- sum(res$Results$Enrichment == "Up-regulated")
        # total_down <- sum(res$Results$Enrichment == "Down-regulated")
        # header <- paste0("Top 10 (Total pathways: Up=", total_up, ", Down=", total_down, ")")
        # 
        # # calculate plot width dynamically based on the number of pathways and size of dots
        # num_pathways <- nrow(res$Results)
        # max_width <- 800  # set your desired maximum width
        # plot_width <- min(max_width, num_pathways * 5)  # adjust the multiplier as needed
        # 
        # # calculate x-axis range dynamically based on the number of pathways
        # x_range <- c(0, max(res$Plot$data$pathway, na.rm = TRUE))
        # 
        # # calculate y-axis range dynamically based on the size of the dots with buffer
        # y_range <- c(min(res$Plot$data$NES) - max(res$Plot$data$size) * 1.5,
        #              max(res$Plot$data$NES) + max(res$Plot$data$size) * 1.5)
        # 
        # plot results
        up_reg <- res$Results[res$Results$Enrichment == "Up-regulated",]
        up_reg <- up_reg[order(up_reg$NES,decreasing=TRUE),]
        top_upreg <- head(up_reg, 10)
        down_reg <- res$Results[res$Results$Enrichment == "Down-regulated",]
        down_reg <- down_reg[order(down_reg$NES,decreasing=FALSE),]
        top_downreg <- head(down_reg, 10)
        
        
        df <- rbind(top_upreg, top_downreg, stringsAsFactors = FALSE)
        
        gsea_plot <- ggplot(df, aes(x = NES, y = pathway, color = padj, size = size)) + 
          geom_point(stat = 'identity') + 
          xlab("Ratio") + ylab("Pathway") + ggtitle("GSEA Plot") + 
          theme_bw()
        
        output$gsea_plot <- renderPlotly(gsea_plot)
        shinyjs::disable("message")
      }
    )
  })
}

# Run application 
shinyApp(ui = ui, server = server)

  