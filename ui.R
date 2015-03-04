
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Functional Identification of Target by Expressional Proteomics"),
  
  # Sidebar with a slider input for number of bins
  sidebarPanel(
    textInput('cell_r', 'Cell RegExp', "LFQ.intensity.([^_]+)_([^_]+)_.*"),
    textInput('drug_r', 'Drug RegExp', "LFQ.intensity.[^_]+_([^_]+)_.*"),
    textInput('ctrl', 'Control', 'CTRL'),
    tags$hr(),
#    checkboxInput('header', 'Header', TRUE),
    selectInput('sep', 'Separator',
                 list(Tab='\t',
                   Comma=',',
                   Semicolon=';'),
                 '\t'),
    tags$hr(),
    fileInput('file1', 'Choose CSV File',
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv','.txt'))
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    h3('Result'),
    dataTableOutput('contents')
  )
))
