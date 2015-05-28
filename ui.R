#require(HMP)
require(shiny)
require(shinythemes)
require(knitr)
require(rmarkdown)
require(msWaldHMP)

#if (no <- require(msWaldHMP))
#{
##  install.packages("./www/msWaldHMP_0.2.tar.gz", repos = NULL, type = "code")
#  owd <- getwd()
#  setwd(tempdir())
#  devtools::install_github(repo = "mafed/msWaldHMP")
#  setwd(owd)
#} else {}

#setwd("C:/Users/Federico/Documents/SharedDocs/Dropbox/Microbiome/Architect/") 
#runApp("shinyMB")
#runApp("../shinyMB")

shinyUI(fluidPage(theme = shinytheme("spacelab"), 
        titlePanel(h1("Microbiome: Power & Sample Sizes Determination")),
        title = "Microbiome Sample Size Determination",
        sidebarLayout(
            sidebarPanel(
                fluidRow(
                    h3("Generate R.A.D. for controls", align = "center"),
                    column(12, align = "center",
                        selectInput(
                            inputId = "kindPiOne",
                            label = "Shape",
                            choices = c(
                                "Stool",
                                "Anterior_nares",
                                "Attached_Keratinized_gingiva",
                                "Buccal_mucosa",
                                "Hard_palate",
                                "Left_Antecubital_fossa",
                                "Left_Retroauricular_crease",
                                "Mid_vagina",
                                "Palatine_Tonsils",
                                "Posterior_fornix",
                                "Right_Antecubital_fossa",
                                "Right_Retroauricular_crease",
                                "Saliva",
                                "Subgingival_plaque",
                                "Supragingival_plaque",
                                "Throat",
                                "Tongue_dorsum",
                                "Vaginal_introitus",
                                "Geometric (Theoretical)", 
                                "Stick-breaking (Theoretical)", 
                                "Two-pieces Geometric (Theoretical)", 
                                "Exponential (Theoretical)")
                        )
                    ),
                    column(4, align = "center",
                        selectInput(
                            inputId = "strata",
                            label = "Stratification",
                            choices = c("NO", "YES"))
                    ),
                    column(4, align = "center",
                        numericInput("numOTUs", 
                            label = "# of OTUs", 
                            value = 50, min = 5, max = 200)
                    ),
                    column(4, align = "center",
                        numericInput(inputId = "theta", 
                            label = "Overdispersion", 
                            value = 0.01, min = 0.005, max = 1, step = 0.005) 
                    ),
#                    column(width = 4, 
#                        numericInput(inputId = "totCounts", 
#                            label = "Library size", 
#                            value = 1e4, min = 100)
#                    ), 
                    column(12, 
                        h3("OR", align = "center")
                    ),
                    column(7, 
                        fileInput("userFiles", label = "Load your data", multiple = TRUE), 
                        helpText(paste0(
                                "If stratification is 'YES' the strata variable ",
                                "has to be in the second column." )
                        )
                    ),
#                    column(12, h4("___")), 
                    column(5, align = "right", 
                        actionButton("reset", label = "Remove File"),
                        helpText(
                            "Click to ignore the uploaded files (only possible once)."
                        )
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
                    )#,
#                    column(12, 
#                        h3("Now choose the theta parameter and library sizes", 
#                            align = "center"), 
#                        helpText(paste0("In case of enterotype stratification, ", 
#                                "sample sizes will respect ", 
#                                "proportions found in real data."))
                    ##                        h3("Dirichlet-Multinomial data generation", align = "center")
#                    ), 
#                    column(12, 
#                        helpText(paste0("In case of stratification, ", 
#                                "sample sizes will respect ", 
#                                "proportions found in real data."))
##                        h3("Dirichlet-Multinomial data generation", align = "center")
#                    ) 
#                    column(width = 12, 
#                        numericInput(inputId = "seed", 
#                            label = "Seed for the RNG", 
#                            value = 12345, min = 100, 
#                            max = round(.Machine$integer.max ^.5))
#                    ),
                ),            # two divided fluidRows
                fluidRow(
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
                    column(width = 4, 
                        numericInput(inputId = "n1", 
                            label = "Sample Size 1", 
                            value = 30, min = 3) 
                    ),
                    column(width = 4, 
                        numericInput(inputId = "n2", 
                            label = "Sample Size 2", 
                            value = 30, min = 3)
                    ),
                    column(width = 4, 
                        h6("."),
                        actionButton("mcStart", label = "Compute Power")
                    ),
                    column(width = 12, align = "center", 
                        h4("Or select a range of sample sizes")
                    ),
                    column(width = 10, align = "center", 
                        sliderInput(inputId = "sampleSizes",
                            label = "Select Range",
                            min = 5, max = 200, value = c(10, 40), step = 5
                        ), 
                        helpText(paste0("In this case the two datasets have the ",
                                "same sample size"))
                    ),
                    column(width = 2, 
                        h5("."),
                        actionButton("powerSimStart", label = "START"),
                        hr()
                    ), 
                    column(width = 12,
                        h3(paste0("Report Generation"), align = "center")
                        ), 
                    column(width = 9, 
                        h4(paste0("Click on the right to generate a PDF report with ",
                                "all inputs and outputs from the analysis."), 
                            align = "left") 
                    ), 
                    column(width = 3, 
                        downloadButton("downloadReport", label = "Download") 
                    )
                )# END - fluidRow
            ),# END - sidebarPanel 
            
            mainPanel(
                tabsetPanel(id = "inTabs", 
                    tabPanel(title = "Settings", 
                        fluidRow(
                            column(9, 
                                h3("Desired R.A.D.s", align = "center"),
#                                plotOutput("piPlot"),
#                                h3("One Generated R.A.D.s", align = "center"),
                                plotOutput("estPiPlot"),
                                plotOutput("estPiPlot2"),
                                plotOutput("estPiPlot3")
                            ), 
                            column(3, 
#                                p(img(src = "hmp_logo_NIH_retina.png", height = 200), 
#                                    align = "center"),
                                h3("Wald-type test", align = "left"), 
#                                textOutput("hmpDiffTest")
                                verbatimTextOutput("hmpDiffTest"),
#                                h3("One Generated R.A.D.s", align = "center"),
                                h3("Resulting Power"),
                                verbatimTextOutput("singleMcWaldResults") 
                            )
                        )# END - fluidRow
                    ), # END - tabPanel
                    tabPanel(title = "Power vs. Sample Size",
                        fluidRow(
                            h3("Power vs. Sample Size", align = "left"),
                            plotOutput("resPowSimPlot", width = "100%", height = "550px",
                                click = "powPlotClick")#,
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



