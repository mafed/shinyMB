#require(HMP)
require(shiny)
require(shinythemes)
require(knitr)
require(rmarkdown)
# devtools::install_github(repo = "mafed/msWaldHMP")
require(msWaldHMP)

#setwd("C:/Users/Federico/Documents/SharedDocs/Dropbox/Microbiome/Architect/") 
#runApp("shinyMB")
#runApp("../shinyMB")
#options(error = recover)

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
                        h3("OR Upload File", align = "center")
                    ),
#                    column(12, h4("___")), 
                    column(3, align = "left", 
#                        helpText(paste0(
#                                "If stratification is 'YES' the strata variable ",
#                                "has to be in the second column." )
#                        ),
                        radioButtons('sep', 'Separator', inline = FALSE,
                            c(Comma=',', Semicolon = ';', Tab='\t'),
                            selected = ','),
                        checkboxInput("strataRow", "Strata Factor in 2nd row", 
                            value = FALSE)
                    ),
                    column(3, align = "left", 
                        radioButtons('quote', 'Quote', inline = FALSE,
                            c(None='',
                                'Double (\")'='"',
                                'Single (\')'="'"),
                            selected = '"')
                    ),
                    column(6, 
                        fileInput("userFiles", label = "Load your data", 
                            multiple = FALSE), 
                        actionButton("reset", label = "Remove File"),
                        helpText("Ignore uploaded file (only once).")
                    ),
#                    column(6, align = "left",
#                    ), 
#                    column(6, align = "left",
#                        ), 
#                    column(4, 
#                        h4("___"),
#                        helpText(paste0(
#                                "Click to ignore the uploaded files and go back to ",
#                                "the simulated abundance curves (only possible once)."
#                            ))
#                    ),
                    column(12, 
#                        h4("Differential Abundance", align = "center")
                        tags$hr(),
                        h3("Define Alternative", align = "center")
                    ),
                    column(3, 
                        numericInput("diffOTUs1", 
                            label = "The", 
                            value = 2, min = 0),
                        numericInput("diffOTUs2", 
                            label = "The next", 
                            value = 2)
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
                            value = 15, min = -300, max = 300),
                        numericInput("relAbund2", 
                            label = "OTUs increased by (%)", 
                            value = 20, min = -300, max = 300)
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
                    column(width = 4, 
                        numericInput(inputId = "MC", 
                            label = "# replications", 
                            value = 100, min = 50, 
                            max = 1e4)
                    ),
                    column(width = 5, 
                        numericInput(inputId = "alpha", 
                            label = "Significance level", 
                            value = 0.1, min = 1e-3, step = 1e-3,
                            max = 0.25)
                    ),
                    column(width = 3, 
                        checkboxInput(inputId = "wmwTest", 
                                label = "Add WMW test (slow)", 
                                value = FALSE)
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
                            min = 5, max = 200, value = c(20, 60), step = 5
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
##                                p(img(src = "hmp_logo_NIH_retina.png", height = 200), 
##                                    align = "center"),
#                               h3("Single Dataset Wald", align = "left"), 
##                                textOutput("hmpDiffTest")
#                               verbatimTextOutput("hmpDiffTest"),
##                                h3("One Generated R.A.D.s", align = "center"),
                                h3("Resulting Power"),
                                verbatimTextOutput("singleMcWaldResults") 
#                                h4(tableOutput("singleMcWaldResults")) 
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
                    ), # END - tabPanel
                    tabPanel(title = "How to use it",
                        fluidRow(
                            column(8, align = "center", 
                                h2("Instructions", align = "center")),
                            column(10, align = "left", 
                                h3(paste0("Select the shape of the ",
                                        "Ranked Abundance Distribution"))
#                                tags$hr()
                            ),
                            column(3, align = "left", #offset = 1,
                                h4(paste0("Choose if you want to use stratification ", 
                                    "or not ('YES' or 'NO')"))
                            ),
                            column(2, align = "left",
                                h4("Number of most abundant OTUs to consider")
                            ),
                            column(3, align = "left",
                                h4(paste0("If a theoretical shape is selected, ",
                                    "overdispersion parameter can be chosen."))
                            ),
                            column(9, 
                                tags$hr(),
                                h3("OR upload your own file", align = "left"),
                                h4(paste0(
                                        "Required file format: ",
                                        "each row is an OTU and each column a ", 
                                        "sample. Stratifying variable, if present, ",
                                        "has to be in the second row."))
                            ),
                            column(6, align = "left", 
                                h4(paste0(
                                        "Specify before uploading: separator ",
                                        "character used in the file, ",
                                        "character defining quotes, ",
                                        "and presence/absence of stratifying variable."))
                            ),
                            column(3, 
                                h4(paste0(
                                        "Click on 'Remove File' to ignore the ",
                                        "uploaded file and go back to the selected ",
                                        "shape (only possible once)."))
                            ),
                            column(10, 
                                tags$hr(),
                                h3("Now define the alternative", align = "left")
                            ),
#                            column(4, h3(".")),
                            column(4, 
                                h4(paste0(
                                        "How many OTUs are to be changed? Two ",
                                        "separate sets of OTUs can be specified, with ",
                                        "different characteristics (see the following).",
                                        "Defaults are 2 and 3 respectively."))
                            ),
                            column(3, 
                                h4(paste0(
                                        "Are these OTUs the 'most abundant' or the ",
                                        "'least abundant' among the ones displayed?"))
                            ),
                            column(3,
                                h4(paste0(
                                        "How much should these OTUs be increased/decreased, ",
                                        "in percentage, with respect to the null ",
                                        "hypothesis (controls)?"))
                            )#,
                        ),            # two divided fluidRows
                        fluidRow(
                            column(width = 10, 
                                tags$hr(),
                                h3("MonteCarlo Simulation Settings", align = "left")
                            ),
                            column(width = 4, 
                                h4(paste0(
                                        "Number of MonteCarlo replications (simulated ",
                                        "datasets) on which power calculations are based."))
                            ),
                            column(width = 2, 
                                h4(paste0(
                                        "Significance level to consider (type I error)."
                                    ))
                            ),
                            column(width = 4, 
                                h4(paste0(
                                        "Compute also a rank-based test ",
                                        "(Wilcoxon-Mann-Whiney), performed on each OTU ",
                                        "individually (with Benjamini-Hochberg ",
                                        "multiplicity correction). It takes much more ",
                                        "to compute."
                                    ))
                            ),
                            column(width = 10, align = "left", 
                                tags$hr(),
                                h3(paste0("First use of the application: with ",
                                        "specific sample sizes."))
                            ),
                            column(width = 3, 
                                h4(paste0("Sample size for controls"))
                            ),
                            column(width = 3, 
                                h4(paste0("Sample size for cases"))
                            ),
                            column(width = 4, 
                                h4(paste0("Compute the power with specified sample sizes"))
                            ),
                            column(width = 10, align = "left", 
                                tags$hr(),
                                h3(paste0(
                                        "Second use of the application: with a ",
                                        "range of sample sizes."))
                            ),
                            column(width = 6, align = "left", 
                                h4(paste0("In this case the two datasets have the ",
                                        "same sample size"))
                            ),
                            column(width = 4, 
                                h4(paste0("Start the computations, resulting graph is on ", 
                                        "the tab 'Power vs. Sample Size'"))
                            ), 
                            column(width = 10,
                                tags$hr(),
                                h3(paste0(
                                        "At the end you can generate a report ", 
                                        "containing all settings and results."))
                            ) 
                        )# END - fluidRow
                    )# END - tabPanel
                )# END - tabsetPanel
            ) # END - mainPanel
        )# END - sidebarLayout
    )# END - fluidPage
)# END - shinyUI




### Notes: to implement
# count only rejections among the differentially abundant
# 1. individual OTU power averaged
# 2. min(adj.p.val) < alpha
# 3. TPR: avg(#rejections among DA / # true DA)




