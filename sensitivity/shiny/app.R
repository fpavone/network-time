library(shiny)
library(visNetwork)

load('ex.Rdata')

ui <- fluidPage(
  
  # App title ----
  titlePanel("Signal-related sensitivity"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      titlePanel("Visualization"),
      textInput(inputId = "seed",
                label = "Network-node seed",
                value = 32423),
      selectInput(inputId = "fill",
                  label = "Coloring attribute",
                  choices = names(meta_data)[-1],
                  selected = 'Degree'),
      
      titlePanel("Data simulation parameters"),
      selectInput(inputId = "tau",
                  label = "Tau",
                  choices = tau_seq,
                  selected = 1),
      
      selectInput(inputId = "kappa",
                  label = "Kappa",
                  choices = kappa_seq,
                  selected = 1),
      
      selectInput(inputId = "alpha",
                  label = "Alpha",
                  choices = alpha_seq,
                  selected = 1),
      
      
      
      #shinythemes::themeSelector(),
      # checkboxGroupInput(inputId = "rows",
      #                    label = "Rows",
      #                    choices = c("Alpha" = "alpha", 
      #                                "Kappa" = "kappa",
      #                                "Tau" = "tau"),
      #                    selected = "alpha"),
      # checkboxGroupInput(inputId = "cols",
      #                    label = "Columns",
      #                    choices = c("Alpha" = "alpha", 
      #                                "Kappa" = "kappa",
      #                                "Tau" = "tau"),
      #                    selected = c("kappa","tau")),
      # Input: Slider for the number of bins ----
      
      ),

    
    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Histogram ----
      visNetworkOutput("network"),
      plotOutput(outputId = "sensitivity",
                 width = "100%",
                 height = "500px")
    )
  ),
  fluidRow(
    includeMarkdown("notation.md")
  )
)


server <- function(input, output) {
  
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$network <- renderVisNetwork({
    color_palette <- col_numeric("viridis", domain = range(pull(meta_data,input$fill)))
    nodes <- data.frame(id = 1:n_nodes,
                        label = as.character(1:n_nodes),
                        color = color_palette(pull(meta_data,input$fill)))
    edges <- data.frame(from = graph_exp$E[,1],
                        to = graph_exp$E[,2],
                        color = 'black')
    visNetwork(nodes, edges) %>%
      visIgraphLayout(randomSeed = as.numeric(input$seed),
                      type = "full") %>%
      visNodes(#shape = "ellipse",
               font = list(color = 'black',
                           size = 50,
                           strokeWidth = 5,
                           strokeColor = "white"))
  })
  output$sensitivity <- renderPlot({
    data_plot %>%
      filter(tau == input$tau,
             kappa == input$kappa,
             alpha == input$alpha) %>%
      ggplot(aes(x = factor(node, levels = 1:n_nodes), 
                 y = abs(value))) +
      geom_boxplot(aes(fill = !!sym(input$fill))) +
      scale_fill_viridis_c() +
      facet_grid(rows = vars(what), scale = 'free') +
      labs(x = 'Node') +
      theme(axis.title.y = element_blank())
    
  })
  
}

shinyApp(ui = ui, server = server)
