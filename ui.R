
require(shiny)
require(shinythemes)
require(knitr)
#require(HMP)
#require(multSampWald)
require(msWaldHMP)
#setwd("C:/Users/Federico/Documents/SharedDocs/Dropbox/Microbiome/Architect/") 
#runApp("shinyMB-0.2")

shinyUI(fluidPage(theme = shinytheme("spacelab"), 
        titlePanel(h1("Microbiome: Power & Sample Sizes Determination")),
        title = "Microbiome Sample Size Determination",
        sidebarLayout(
            sidebarPanel(
                fluidRow(
                    h3("Generate R.A.D.", align = "center"),
                    column(4, align = "center",
                        selectInput(
                            inputId = "entero",
                            label = "Stratification",
                            choices = c("NO", "YES"))
                    ),
                    column(4, align = "center",
                        selectInput(
                            inputId = "kindPiOne",
                            label = "Shape",
                            choices = c("Geometric", "Broken stick", 
                                "Two-pieces Geometric", "Exponential"))
                    ),
                    column(4, align = "center",
                        numericInput("numOTUs", 
                            label = "# of OTUs",
                            value = 50, min = 5, max = 200)
                    ),
                    column(12, 
                        h3("OR", align = "center")
                    ),
                    column(6, 
                        fileInput("userFiles", label = "Load your data", multiple = TRUE), 
                        helpText(paste0(
                                "In case of Enterotype stratification you need to ",
                                "select 3 files." )
                        )
                    ),
#                    column(12, h4("___")), 
                    column(6, align = "center", 
                        actionButton("reset", label = "Remove File"),
                        helpText(paste0(
                                "Click to ignore the uploaded files and go back to ",
                                "the simulated abundance curves (only possible once)."
                            ))
#                        helpText(paste0(
#                                "NOTE: you can upload files only once. ",
#                                "To load other files refresh the browser.")
#                        )
                    ),
#                    column(4, 
##                        h4("___"),
#                        helpText(paste0(
#                                "Click to ignore the uploaded files and go back to ",
#                                "the simulated abundance curves (only possible once)."
#                            ))
#                    ),
                    column(12, 
#                        h4("Differential Abundance", align = "center")
                        h3("Define Alternative", align = "center")
                    ),
                    column(3, 
                        numericInput("diffOTUs1", 
                            label = "The", 
                            value = 2, min = 0),
                        numericInput("diffOTUs2", 
                            label = "The next", 
                            value = 3)
                    ),
                    column(4, 
                        selectInput(
                            inputId = "mostLeastAb1", label = "most/least abundant", 
                            choices = c("most abund.", "least abund.")),
                        selectInput(
                            inputId = "mostLeastAb2", label = "most/least abundant", 
                            choices = c("most abund.", "least abund."))
                    ),
                    column(5,
                        numericInput("relAbund1", 
                            label = "OTUs increased by (%)", 
                            value = 33, min = -300, max = 300),
                        numericInput("relAbund2", 
                            label = "OTUs increased by (%)", 
                            value = 25, min = -300, max = 300)
                    ),
#                    column(12, 
#                        h3("Now choose the theta parameter and library sizes", 
#                            align = "center"), 
#                        helpText(paste0("In case of enterotype stratification, ", 
#                                "sample sizes will respect ", 
#                                "proportions found in real data."))
##                        h3("Dirichlet-Multinomial data generation", align = "center")
#                    ), 
                    column(width = 6, 
                        numericInput(inputId = "theta", 
                            label = "Overdispersion", 
                            value = 0.01, min = 0.005, max = 1, step = 0.005) 
                    ),
                    column(width = 6, 
                        numericInput(inputId = "totCounts", 
                            label = "Library size", 
                            value = 1e4, min = 100)
                    ), 
                    column(12, 
                        helpText(paste0("In case of enterotype stratification, ", 
                                "sample sizes will respect ", 
                                "proportions found in real data."))
#                        h3("Dirichlet-Multinomial data generation", align = "center")
                    ), 
#                    column(width = 12, 
#                        numericInput(inputId = "seed", 
#                            label = "Seed for the RNG", 
#                            value = 12345, min = 100, 
#                            max = round(.Machine$integer.max ^.5))
#                    ), 
                    column(width = 12, 
                        h3("MonteCarlo Simulation Settings", align = "center")
                    ),
#                    column(width = 12, 
#                        helpText(paste0("With these settings you can ",
#                                "now run MonteCarlo simulations"))
#                    ),
                    column(width = 6, 
                        numericInput(inputId = "MC", 
                            label = "# of replications", 
                            value = 100, min = 100, 
                            max = 1e4)
                    ),
                    column(width = 6, 
                        numericInput(inputId = "alpha", 
                            label = "Significance level", 
                            value = 0.1, min = 1e-3, step = 1e-3,
                            max = 0.25)
                    ),
                    column(width = 12, align = "center", 
                        h4("Compute power with these sample sizes")
                    ),
                    column(width = 6, 
                        numericInput(inputId = "n1", 
                            label = "Sample Size 1", 
                            value = 30, min = 3) 
                    ),
                    column(width = 6, 
                        numericInput(inputId = "n2", 
                            label = "Sample Size 2", 
                            value = 30, min = 3)
                    ),
                    column(width = 12, 
                        actionButton("mcStart", label = "Compute Power")
                    ),
                    column(width = 12, align = "center", 
                        h4("Or select a range of sample sizes"),
                        sliderInput(inputId = "sampleSizes",
                            label = "Select Range",
                            min = 5, max = 200, value = c(10, 40), step = 5
                        ), 
                        helpText(paste0("In this case the two datasets have the ",
                                "same sample size"))
                    ),
                    column(width = 12, 
                        actionButton("powerSimStart", label = "Start Simulations")
                    ), 
                    column(width = 12, 
                        h4(paste0("If you click below you generate a PDF report with ",
                                "all inputs and outputs from the analysis."), 
                            align = "center"),
                        downloadButton("downloadReport", label = "Download")
                    )
                )# END - fluidRow
            ),# END - sidebarPanel 
            
            mainPanel(
                tabsetPanel(id = "inTabs", 
                    tabPanel(title = "Settings", 
                        fluidRow(
                            column(8, 
                                h3("Desired R.A.D.s", align = "center"),
#                                plotOutput("piPlot"),
#                                h3("One Generated R.A.D.s", align = "center"),
                                plotOutput("estPiPlot"),
                                plotOutput("estPiPlot2"),
                                plotOutput("estPiPlot3")
                            ), 
                            column(4, 
                                p(img(src = "hmp_logo_NIH_retina.png", height = 200), 
                                    align = "center"),
                                h3("Wald-type test", align = "left"), 
#                                textOutput("hmpDiffTest")
                                verbatimTextOutput("hmpDiffTest"),
#                                h3("One Generated R.A.D.s", align = "center"),
                                h3("Resulting Power"),
                                verbatimTextOutput("singleMcWaldResults") 
                            )
                        )# END - fluidRow
                    ), # END - tabPanel
                    tabPanel(title = "Results",
                        fluidRow(
                            h3("Power vs. Sample Size", align = "left"),
                            plotOutput("resPowSimPlot", width = "90%", height = "550px",
                                clickId = "powPlotClickCoord")#,
#                            column(8,
#                                h3("Rejection proportion", align = "left"),
#                                verbatimTextOutput("mcHmpWaldResPrint")
#                            )
                        )# END - fluidRow
                    )# END - tabPanel
                )# END - tabsetPanel
            )
        )# END - sidebarLayout
    )# END - fluidPage
)# END - shinyUI

