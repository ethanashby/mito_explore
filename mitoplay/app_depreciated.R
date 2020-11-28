library(shiny)
library(tidyverse)
library(parallel)
library(seqinr)
library(data.table)
library(ggplot2)
library(stringi)
library(stringr)
library(reshape2)
library(entropy)
library(gggenes)

fasta <- read.fasta(file="mitochondrion_sequence.fasta", seqtype="DNA", as.string=TRUE, seqonly=TRUE)

anno<-data.table(read.csv("mito_coords.csv", header=FALSE))
colnames(anno)<-c("Locus", "Start", "Stop", "Feature", "Description", "Ref", "NA")
anno<-setDT(anno[,2:4])
anno_1<-anno[Start<Stop, with=TRUE]
anno_1<-as.data.frame(anno_1)

for (i in 1:3){
  anno_1<-rbind(anno_1, c("Start"=as.numeric(anno[Start>Stop, with=TRUE][i,1]), "Stop"=16569, "Feature"=as.character(anno[Start>Stop, with=TRUE][i,3])))
  anno_1<-rbind(anno_1, c("Start"=0, "Stop"=as.numeric(anno[Start>Stop, with=TRUE][i,2]), "Feature"=as.character(anno[Start>Stop, with=TRUE][i,3])))
}

anno_1<-setDT(anno_1[,c(3,1,2)])
anno_1$Start<-as.numeric(anno_1$Start)
anno_1$Stop<-as.numeric(anno_1$Stop)

# Define UI for app that draws a histogram ----

ui <- fluidPage(
  # App title ----
  titlePanel("Shannon Entropy of Mitochondrial Genome"),
  sidebarPanel(
    selectizeInput(
      inputId = "featureChoice", 
      label = "Select a Feature", 
      choices = unique(anno_1$Feature), 
      selected = "CR:HVS2",
      multiple = TRUE
    ),
    verbatimTextOutput("rage"),
    # checkboxInput("smoothCheckbox", label = "Smooth", value = FALSE),
    # checkboxInput("boxCheckbox", label = "Box", value = FALSE),
    textInput("breaks", "Enter window size(s) separated by a comma...", value="10"),
    actionButton("plotButton", label = "Plot Graphs"),
    # textOutput("breakresult")
    # dataTableOutput("data")
  ),
  mainPanel(
    plotOutput(outputId = "plot", width = "100%"),
  ),
  verbatimTextOutput("info"),
  plotOutput(outputId = "genomePlot", 
             # width = "100%", 
             dblclick = "plot1_dblclick",
             click = "plot_click",
             hover = "plot_hover",
             brush = brushOpts(
               id = "plot1_brush",
               resetOnNew = TRUE)
  )
)
server <- function(input, output) {
  output$breakresult <- renderPrint({ 
    # nums <- extract(input$breaks)
    # nums
    # print(input$do)
  })
  
  extract <- function(text) {
    text <- gsub(" ", "", text)
    split <- strsplit(text, ",", fixed = FALSE)[[1]]
    as.numeric(split)
  }
  
  get_windowNames <- function(text_w){
    window_choices <- c("index")
    for(i in rev(text_w)){
      # print(paste0("w",i))
      window_choices <- rbind(paste0("w",i), window_choices)
    }
    window_choices[,1]
  }
  
  genomeDataInput <- eventReactive(input$plotButton, {
    nums <- extract(input$breaks)
    entropies <- mclapply(nums, function(w){
      addon <- substring(fasta, 1, w)
      fasta_added <- paste(fasta, addon, sep="")
      string_chunks <- mclapply(1:nchar(fasta), function(i){
        s <- substring(fasta_added, i, i + w - 1) # %>% str_split("") %>% table()
        c(str_count(s, "A"), str_count(s, "C"), str_count(s, "G"), str_count(s, "T"))
      }, mc.cores=4)
      get_entropies <- mclapply(string_chunks, function(i){
        entropy(i)
      }, mc.cores=4)
      unlist(get_entropies)
    })
    entropies <- as.data.frame(do.call(cbind, entropies)) %>% mutate(index=1:16569)
    colnames(entropies) <- get_windowNames(nums)
    entropies_df <- reshape::melt(entropies, id=c("index"))
    plot_var = FALSE
    entropies_df
  })
  
  output$plot <- renderPlot({
    feature_row <- anno_1 %>% filter(Feature==input$featureChoice)
    start <- feature_row$Start
    end <- feature_row$Stop
    filtered_data <- genomeDataInput() %>% filter(start <= index & index <= end)
    label_pos <- summary(filtered_data$value)[3]
    ggplot(anno_1) +
      geom_gene_arrow(mapping=aes(xmin=start, xmax=end, y=label_pos, fill=Feature, label=Feature)) +
      geom_gene_label(mapping=aes(xmin=start, xmax=end, y=label_pos, fill=Feature, label=Feature), 
                      align = "left") + 
      geom_smooth(data=filtered_data, aes(x = index, y = value, color=variable)) +theme(legend.position="none") +ylab("Shannon Entropy")+xlab("MTchr position")
    # geom_smooth(aes(x=index, y=value, color=variable)) 
    # coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) + 
  })
  
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  output$genomePlot <- renderPlot({
    # print("plotting!")
    nums <- extract(input$breaks)
    print(ranges$x)
    # genomeDataInput() %>% filter(variable %in% get_windowNames(nums)) %>% 
    #   ggplot() + geom_smooth(aes(x=index, y=value, color=variable))
    filtered_data <- genomeDataInput() %>% filter(variable %in% get_windowNames(nums))
    label_pos <- summary(filtered_data$value)[3]
    ggplot(anno_1) +
      geom_gene_arrow(mapping=aes(xmin=Start, xmax=Stop, y=label_pos, fill=Feature, label=Feature)) +
      geom_gene_label(mapping=aes(xmin=Start, xmax=Stop, y=label_pos, fill=Feature, label=Feature), 
                      align = "left") + 
      geom_smooth(data=filtered_data, aes(x = index, y = value, color=variable)) +theme(legend.position="none")+
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
  })
  
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  output$info <- renderText({
    xy_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("x=", round(e$x, 1), " y=", round(e$y, 1), "\n")
    }
    xy_range_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1), 
             " ymin=", round(e$ymin, 1), " ymax=", round(e$ymax, 1))
    }
    paste0(
      "click: ", xy_str(input$plot_click),
      "dblclick: ", xy_str(input$plot1_dblclick),
      "hover: ", xy_str(input$plot_hover),
      "brush: ", xy_range_str(input$plot1_brush)
    )
  })
}

shinyApp(ui, server)

