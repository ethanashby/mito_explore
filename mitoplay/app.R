library(shiny)
library(entropy)
library(stringr)
library(seqinr)
library(stringi)
library(future.apply)
library(tidyverse)
library(dygraphs)
library(reshape)
library(ggplot2)
library(data.table)
library(gggenes)


####### Import Genome ###########
fasta <- read.fasta(file="mitochondrion_sequence.fasta", seqtype="DNA", as.string=TRUE, seqonly=TRUE)


##### Get Anno_1 for plotting purposes ########
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
anno_1 <- anno_1[abs(anno_1$Stop-anno_1$Start)>50,]

####### Get Feature Coordinates anno without wrap separated #########
anno<-data.table(read.csv("mito_coords.csv", header=FALSE))
colnames(anno)<-c("Locus", "Start", "Stop", "Feature", "Description", "Ref", "NA")
anno<-setDT(anno[,2:4])
anno <- anno[abs(anno$Stop-anno$Start)>50,]


########Create master entropy data frame ##########
get_windowNames <- function(text_w){
  window_choices <- c("index")
  for(i in rev(text_w)){
    # print(paste0("w",i))
    window_choices <- rbind(paste0("w",i), window_choices)
  }
  window_choices[,1]
}

nums <- seq(10, 510, 40)
entropies <- future_lapply(nums, function(w){
  addon <- substring(fasta, 1, w)
  fasta_added <- paste(fasta, addon, sep="")
  string_chunks <- future_lapply(1:nchar(fasta), function(i){
    s <- substring(fasta_added, i, i + w - 1) # %>% str_split("") %>% table()
    c(str_count(s, "A"), str_count(s, "C"), str_count(s, "G"), str_count(s, "T"))
  })
  get_entropies <- future_lapply(string_chunks, function(i){
    entropy(i)
  })
  unlist(get_entropies)
})
entropies <- as.data.frame(do.call(cbind, entropies)) %>% mutate(index=1:16569)
colnames(entropies) <- get_windowNames(nums)
entropies_df <- reshape::melt(entropies, id=c("index"))
# entropies_df


########## Shiny App ###############
# Define UI for app that draws a histogram ----
ui <- fluidPage(
  # App title ----
  titlePanel("Shannon's Entropy on the Mitochondrial Genome"),
  sidebarPanel(
    selectInput("featureChoice", label = h3("Select Feature"), 
                choices = unique(anno$Feature),
                selected = "CR:HVS2"),
    verbatimTextOutput("rage"),
    selectizeInput(
      inputId = "windowChoice", 
      label = "Select a Window Size", 
      choices = unique(seq(10, 510, 40)), 
      selected = 10,
      multiple = TRUE
    ),
    checkboxInput("pointCheck", label = "points", value = FALSE),
    actionButton("plotButton", label = "Plot Graphs"),
    verbatimTextOutput("featureStats"),
  ),
  mainPanel(
    plotOutput(outputId = "plot", width = "100%"),
  ),
  fluidRow(
    column(12, 
           h4("Entropy of Entire Mitochondrial Genome"),
           plotOutput(outputId = "genomePlot", 
                      # width = "100%", 
                      dblclick = "plot1_dblclick",
                      click = "plot_click",
                      hover = "plot_hover",
                      brush = brushOpts(
                        id = "plot1_brush",
                        resetOnNew = TRUE)
           ),
           verbatimTextOutput("info"),
    )
  ),

)
server <- function(input, output) {
  output$breakresult <- renderPrint({})
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
    nums <- as.numeric(input$windowChoice)
    entropies_df %>% filter(variable %in% get_windowNames(nums))
  })
  output$plot <- renderPlot({
    feature_row <- anno %>% filter(Feature==input$featureChoice)
    start <- feature_row$Start
    end <- feature_row$Stop
    title <- paste0("Entropy of Feature: ", input$featureChoice)
    p <- ggplot(anno)
    if(start > end){
      chunk1 <- genomeDataInput() %>% filter(index >= start)
      chunk2 <- genomeDataInput() %>% filter(index <= end)
      # filtered_data <- 
      label_pos <- summary(chunk1$value)[3]
      p <- p + geom_gene_arrow(mapping=aes(xmin=start, xmax=16596, y=label_pos, fill=Feature, label=Feature)) +
        geom_gene_arrow(mapping=aes(xmin=0, xmax=end, y=label_pos, fill=Feature, label=Feature)) +
        geom_gene_label(mapping=aes(xmin=start, xmax=16596, y=label_pos, fill=Feature, label=Feature), 
                        align = "left") + 
        geom_gene_label(mapping=aes(xmin=0, xmax=end, y=label_pos, fill=Feature, label=Feature), 
                        align = "left") + 
        geom_smooth(data=chunk1, aes(x = index, y = value, color=variable)) +
        geom_smooth(data=chunk2, aes(x = index, y = value, color=variable)) +
        theme(legend.position="none", plot.title = element_text(hjust = 0.5)) + 
        ggtitle(title) +
        xlab("Index") + ylab("Entropy")
      
      if(input$pointCheck){
        p <- p + geom_point(data=chunk1, aes(x = index, y = value, color=variable)) + 
          geom_point(data=chunk2, aes(x = index, y = value, color=variable))
      }
    } else {
      filtered_data <- genomeDataInput() %>% filter(start <= index & index <= end)
      label_pos <- summary(filtered_data$value)[3]
      p <- ggplot(anno_1) + geom_gene_arrow(mapping=aes(xmin=start, xmax=end, y=label_pos, fill=Feature, label=Feature)) +
        geom_gene_label(mapping=aes(xmin=start, xmax=end, y=label_pos, fill=Feature, label=Feature), 
                        align = "left") + 
        geom_smooth(data=filtered_data, aes(x = index, y = value, color=variable))+
        theme(legend.position="none", plot.title = element_text(hjust = 0.5)) + 
        ggtitle(title) +
        xlab("Index") + ylab("Entropy")
        if(input$pointCheck){
          p <- p + geom_point(data=filtered_data, aes(x = index, y = value, color=variable))
        }
    }
    p
  })
  
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  output$genomePlot <- renderPlot({
    nums <- as.numeric(input$windowChoice)
    filtered_data <- genomeDataInput() #%>% filter(variable %in% get_windowNames(nums))
    label_pos <- summary(filtered_data$value)[3]
    p <- ggplot(anno_1) +
      geom_gene_arrow(mapping=aes(xmin=Start, xmax=Stop, y=label_pos, fill=Feature, label=Feature)) +
      geom_gene_label(mapping=aes(xmin=Start, xmax=Stop, y=label_pos, fill=Feature, label=Feature), 
                      align = "left") + 
      geom_smooth(data=filtered_data, aes(x = index, y = value, color=variable)) +
      theme(legend.position="none", plot.title = element_text(hjust = 0.5))+
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +
      ggtitle("Entropy of Mitochondial Genome") +
      xlab("Index") + ylab("Entropy")
    if(input$pointCheck){
      p <- p + geom_point(data=filtered_data, aes(x = index, y = value, color=variable))
    }
    p
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
  
  output$featureStats <- renderText({
    feature_row <- anno %>% filter(Feature==input$featureChoice)
    start <- feature_row$Start
    end <- feature_row$Stop
    filtered_data <- genomeDataInput() %>% filter(start <= index & index <= end)
    maxE <- 0
    avgE <- 0
    minE <- 0
    if(start > end){
      chunk1 <- genomeDataInput() %>% filter(index >= start)
      chunk2 <- genomeDataInput() %>% filter(index <= end)
      if(max(chunk1$value) > max(chunk2$value)){
        maxE <- max(chunk1$value)
      }else {
        maxE <- max(chunk2$value)
      }
      if(mean(chunk1$value) > mean(chunk2$value)){
        avgE <- mean(chunk1$value)
      }else {
        avgE <- mean(chunk2$value)
      }
      if(min(chunk1$value) > min(chunk2$value)){
        minE <- min(chunk1$value)
      }else {
        minE <- min(chunk2$value)
      }
    } else {
      maxE <- max(filtered_data$value)
      AvgE <- mean(filtered_data$value)
      minE <- min(filtered_data$value)
    }
    
    # print(AvgE)
    paste0(
      "Feature: ", paste0(input$featureChoice,"\n"),
      "Max Entropy: ", paste0(maxE,"\n"),
      "Average Entropy: ", paste0(AvgE,"\n"),
      "Minx Enropy: ", paste0(minE,"\n")
    )
  })
}

shinyApp(ui, server)