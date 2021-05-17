library(shinydashboard)
library(DT)
library(mdatools)

#sidebar design
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Data", tabName = "reference", icon = icon("table"))
  )
)

#_____________________________________________________________________________________

#body design
body <- dashboardBody(
  tabItems(
    tabItem(tabName = "reference",
            tabsetPanel(
              tabPanel("Upload",
                       fluidRow(
                         column(3,
                                box(width = 12, title = "Upload data",
                                    # Input: Select a file ----
                                    fileInput("file1", "Choose File",
                                              multiple = FALSE,
                                              accept = c("text/csv",
                                                         "text/comma-separated-values,text/plain",
                                                         ".csv"))
                                    ,
                                    
                                    # Horizontal line ----
                                    tags$hr(),
                                    
                                    # Input: File type ----
                                    radioButtons("selectAPI", "Select API",
                                                 choices = c("Acetaminophen" = 'acetaminophen',
                                                             "Albendazole" = 'albendazole',
                                                             "Amoxicillin" = 'amoxicillin',
                                                             "Ampicillin" = 'ampicillin',
                                                             "Azithromycin" = 'azithromycin',
                                                             "Benzyl Penicillin" = 'benzyl penicillin',
                                                             "Ceftriaxone" = 'ceftriaxone',
                                                             "Chloroquine" = 'chloroquine',
                                                             "Ciprofloxacin" = 'ciprofloxacin',
                                                             "Doxycycline" = 'doxycycline',
                                                             "Epinephrine" = 'epinephrine',
                                                             "Ethambutol" = 'ethambutol',
                                                             "Ferrous Sulfate" = 'ferrous sulfate',
                                                             "Hydroxychloroquine" = 'hydroxychloroquine',
                                                             "Isoniazid" = 'isoniazid',
                                                             "Promethazine Hydrochloride" = 'promethazine hydrochloride',
                                                             "Pyrazinamide" = 'pyrazinamide',
                                                             "Rifampicin" = 'rifampicin',
                                                             "RIPE" = 'ripe',
                                                             "Sulfamethoxazole" = 'sulfamethoxazole',
                                                             "Tetracycline" = 'tetracycline'),
                                                 selected = 'acetaminophen'),
                                    
                                    conditionalPanel(
                                      condition = "input.datatype == 'training'",
                                      checkboxInput("rowselect", "Select Data for Model")
                                      
                                    )
                                    # 
                                    # # Input: Checkbox if file has header ----
                                    # checkboxInput("header", "Header", TRUE),
                                    # 
                                    # # Input: Select separator ----
                                    # radioButtons("sep", "Separator",
                                    #              choices = c(Comma = ",",
                                    #                          Semicolon = ";",
                                    #                          Tab = "\t"),
                                    #              selected = ","),
                                    # 
                                    # # Input: Select quotes ----
                                    # radioButtons("quote", "Quote",
                                    #              choices = c(None = "",
                                    #                          "Double Quote" = '"',
                                    #                          "Single Quote" = "'"),
                                    #              selected = '"'),
                                    # 
                                    # # Input: Skip lines in record ----
                                    # numericInput("begin", "Data begins on line number",
                                    #              value = 0),
                                    # 
                                    # # Input: Checkbox if need to transpose dataset ----
                                    # checkboxInput("transpose", "Transpose", TRUE)
                                    
                                    
                                )#end box
                                
                                
                         ), #end column 1
                         column(8,
                                tags$style(type="text/css",
                                           ".shiny-output-error { visibility: hidden; }",
                                           ".shiny-output-error:before { visibility: hidden; }"
                                ),
                                radioButtons("datatype", "Data",
                                             choices = c(Test = "test",
                                                         Training = "training"),
                                             selected = "training",
                                             inline = TRUE),
                                plotOutput("spectra"),
                                br(),
                                DTOutput("contents"),
                                
                                verbatimTextOutput("info")
                                
                         ) #end column 2
                         
                       ) #end fluidRow
                       
              ), #end tabPanel
              tabPanel("Model",
                       fluidRow(
                         column(3,
                                box(width = 12, title = "PLS Parameters",
                                    radioButtons("images", tags$b("Images in Model"),
                                                 choices = c(First = "First",
                                                             Mid = "Mid",
                                                             Last = "Last"),
                                                 selected = "First"),
                                    checkboxInput("center", tags$b("Center Data"), value = TRUE),
                                    checkboxInput("scale", tags$b("Scale Data"), value = TRUE),
                                    checkboxInput("msc", tags$b("MSC"), value = FALSE),
                                    checkboxInput("snv", tags$b("SNV"), value = FALSE),
                                    #checkboxInput("log", tags$b("Log Transform"), value = FALSE),
                                    radioButtons("algorithm", tags$b("PLS Algorithm"),
                                                 choices = c(Kernel = "kernelpls",
                                                             NIPALS = "oscorespls",
                                                             SIMPLS = "simpls"),
                                                 selected = "oscorespls"),
                                    checkboxInput("validation", tags$b("Cross Validation"), value = FALSE)
                                    # conditionalPanel(
                                    #   condition = "input.validation == '1'",
                                    #   sliderInput("sgolay", "Derivative Order", min = 0, max = 2, value = 0),
                                    #   numericInput("porder", "Polynomial Order", value = 1),
                                    #   numericInput("window", tags$b("Window"), value = 15)
                                    # ),
                                    # checkboxInput("filter", tags$b("Savitzky-Golay Filter"), value = FALSE),
                                    # conditionalPanel(
                                    #   condition = "input.filter == '1'",
                                    #   sliderInput("sgolay", "Derivative Order", min = 0, max = 2, value = 0),
                                    #   numericInput("porder", "Polynomial Order", value = 1),
                                    #   numericInput("window", tags$b("Window"), value = 15)
                                    # ),
                                    # checkboxInput("residuals", tags$b("Use Residuals"), value = TRUE),
                                    # numericInput("components", tags$b("Components"), value = 10)
                                )#end box
                                
                         ), #end column 1
                         column(8,
                                plotOutput("biplot"),
                                downloadButton("plot1"),
                                br(),
                                plotOutput("loadplot"),
                                downloadButton("plot2"),
                                br(),
                                plotOutput("coefplot"),
                                downloadButton("plot3"),
                                br(),
                                verbatimTextOutput("plsSummary"),
                                br(),
                                DTOutput("plsLoadings"),
                                downloadButton("downloadLoadings", "Loadings"),
                                br(),
                                br(),
                                DTOutput("plsScores"),
                                downloadButton("downloadScores", "Scores"),
                                # downloadButton("downloadPCAPlot", "Plot"),
                                # downloadButton("downloadPCAData", "Selected Data"),
                                verbatimTextOutput("modelCoefs"),
                                downloadButton("downloadModel", "Coefficients")
                         )
                         # column(5,
                         #        plotOutput("rmsep"),
                         #        br(),
                         #        # plotOutput("compsuggest"),
                         #        # br(),
                         #        plotOutput("biplot"),
                         #        downloadButton("downloadPCAPlot", "Plot"),
                         #        downloadButton("downloadPCAData", "Selected Data")
                         #        # DTOutput("mdist"),
                         #        # downloadButton("data1Download", "Download Analysis")
                         #        
                         # )
                         # , #end column 2
                         # column(5,
                         #        plotOutput("screeplot"),
                         #        br(),
                         #        plotOutput("residuals")
                         # 
                         #        )
                         
                       ) #end fluidRow
                       
              ) #end tabPanel
              ,
              tabPanel("Predict",
                       fluidRow(
                         column(3,
                                box(width = 12, title = "Model Parameters",
                                    numericInput("components", tags$b("Components"), value = 10)
                                    ,
                                    checkboxGroupInput("laneBands", label = tags$b("Select Lanes"),
                                                       choices = list("A"="A", "B"="B","C"="C","D"="D","E"="E","F"="F", 
                                                                      "G"="G","H"="H","I"="I","J"="J","K"="K","L"="L"),
                                                       selected = c("A","B","C","D","E","F","G","H","I","J","K","L"))
                                    ,
                                    checkboxInput("controls", tags$b("Use Controls"), value = FALSE),
                                    downloadButton("data1Download", "Download"),
                                    actionButton("lowRMSEP","Find Lowest RMSEP"),
                                    verbatimTextOutput("minRMSEP")
                                )#end box
                                
                                
                                
                         ), #end column 1
                         column(8,
                                # radioButtons("datatype2", "Data",
                                #              choices = c(Test = "test",
                                #                          Training = "training"
                                #                          ),
                                #              selected = "test",
                                #              inline = TRUE),
                                # plotOutput("predictplot"),
                                # br(),
                                # plotOutput("residuals"),
                                #verbatimTextOutput("lanes")
                                plotOutput("boxplot"),
                                br(),
                                #verbatimTextOutput("predictvalues"),
                                DTOutput("predictvalues")
                                
                                
                         )
                         # , #end column 2
                         # column(5,
                         #        plotOutput("screeplot"),
                         #        br(),
                         #        plotOutput("residuals")
                         # 
                         #        )
                         
                       ) #end fluidRow
                       
              ) #end tabPanel
            ) #end tabsetPanel
    ) #end tabItem
  ) #end tabItems
) #end dasbhoardBody


#_____________________________________________________________________________________

#App arguments
dashboardPage(
  dashboardHeader(title = "PAD PLSR"),
  sidebar,
  body
)
