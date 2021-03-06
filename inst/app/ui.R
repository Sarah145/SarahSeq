#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinycssloaders)
library(plotly)


# Define UI for application that draws a histogram
shinyUI(fluidPage(
    includeCSS("www/custom.css"), # link to css

    
    # page header                  
    titlePanel(tags$div(class='header-panel',  
                        h1("SarahSeq", align = 'center'), 
                        img(src='pal1.png') 
                        ), windowTitle = 'SarahSeq'),

    # sidebar
    sidebarLayout(
        sidebarPanel(
            fileInput("vcf", h4("Choose file"),
                      multiple = F,
                      accept = c(".vcf", ".vcf.gz")),
            tags$hr(),
            withSpinner(textOutput("success"), color="#4b36c7", size = 0.5),
            width = 3
        ),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
                tabPanel("Pick A Gene", 
                         selectInput('goi', label = '', choices = readLines(system.file('gns.txt', package = 'SarahSeq')),
                                     selected = 'BRCA1'),
                         actionButton('goi_go', label = 'GO'),
                         plotlyOutput('goi_plot'),
                         tags$hr(),
                         conditionalPanel(
                             condition = "input.goi_go > 0",
                             downloadButton("download_goi", "Download"),
                             dataTableOutput('goi_table'))),
                tabPanel("Rare Variants",
                         tags$div(style="margin-top:40px; margin-bottom:30px;",
                                  h5('Compare your genome with the UK BioBank to find rare variants.'),
                                  actionButton('rare_var_go', label = 'GO'),
                                  tags$hr(),
                                  tags$br(),
                                  sliderInput('p_threshold', 'Find variants present in less than this percentage of UK BioBank participants:',
                                              min = 0.01, max = 10, value = 1, width = '100%')
                         ),
                         tags$hr(),
                         conditionalPanel(
                             condition = "input.rare_var_go > 0",
                             downloadButton("download_rare_var", "Download"),
                             withSpinner(dataTableOutput('rare_var_df'), color = "#4b36c7")),
                         tags$hr(),
                         conditionalPanel(
                             condition = "input.rare_var_go > 0",
                             uiOutput('plot_title'),
                             plotOutput('rare_var_plot'),
                             radioButtons('pal', '', inline = T, choiceValues = 1:4, selected = 1, 
                                          choiceNames = paste('Palette', 1:4)))
                ),
                tabPanel("Ancestry",
                         tags$div(style="margin-top:40px; margin-bottom:30px;",
                             h5('Compare your genome with genomes from the 1,000 genomes project.'),
                             actionButton('pca_go', label = 'GO')),
                         tags$hr(),
                         plotlyOutput("genomes_1k_pca"),
                         conditionalPanel(
                             condition = "input.pca_go > 0",
                             img(src='legend.png', align = "center", width = "100%"),
                             tags$hr(),
                             h6('For detailed information about what each population code stands for, see', 
                                tags$a(href="https://www.internationalgenome.org/faq/which-populations-are-part-your-study/", "this website.")))
                ),
                tabPanel("Visualisation",
                         tags$div(style="margin-top:40px; margin-bottom:30px;",
                                  h5('View your genome as a Hilbert curve coloured by the density of variants in a given region.'),
                                  actionButton('hilbert_go', label = 'GO')),
                         tags$hr(),
                         radioButtons('hil_pal', '', inline = T, choiceValues = 1:2, selected = 1, 
                                      choiceNames = paste('Palette', 1:2)),
                         withSpinner(plotOutput('HilbertCurve'), color = "#4b36c7"),
                         #plotOutput('hilbert_legd', height = "30px"),
                         downloadButton('download_hil', 'Download Plot'),
                         tags$br(),
                         tags$br()
                )
            )
        ))
    ))

