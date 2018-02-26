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
        titlePanel(h1("Power & Sample Sizes Tool for Case-Control Microbiome Studies")),
        title = "Power&Sample Sizes Tool",
        withMathJax(),
        sidebarLayout(
            sidebarPanel(
                fluidRow(
                    h3("Choose Relative Abundance Curve for Controls", 
                        "(\\(\\pi_1\\))", align = "center"),
                    column(12, align = "center",
                        selectInput(
                            inputId = "kindPiOne",
                            label = "Shape",
                            choices = c(
                                "StoolEntero",
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
                        h3("OR Estimate It From Your Data", align = "center")
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
                    #column(6, align = "left",
                    #), 
                    #column(6, align = "left",
                    #    ), 
                    #column(4, 
                    #    h4("___"),
                    #    helpText(paste0(
                    #            "Click to ignore the uploaded files and go back to ",
                    #            "the simulated abundance curves (only possible once)."
                    #        ))
                    #),
                    column(12, 
                        #                        h4("Differential Abundance", align = "center")
                        hr(),
                        h3("Define Relative Abundance Curve for Cases",
                            "(\\(\\pi_2\\))", align = "center")
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
                            choices = c("most abundant", "least abundant")),
                        selectInput(
                            inputId = "mostLeastAb2", label = "most/least abundant", 
                            choices = c("most abundant", "least abundant"))
                    ),
                    column(4,
                        numericInput("relAbund1", 
                            label = "OTUs increased by (%)", 
                            value = 20, min = -200, max = 200),
                        numericInput("relAbund2", 
                            label = "OTUs increased by (%)", 
                            value = 40, min = -200, max = 200)
                    )#,
                # column(12, 
                #     h3("Now choose the theta parameter and library sizes", 
                #         align = "center"), 
                #     helpText(paste0("In case of enterotype stratification, ", 
                #             "sample sizes will respect ", 
                #             "proportions found in real data."))
                ##     h3("Dirichlet-Multinomial data generation", align = "center")
                # ), 
                # column(12, 
                #     helpText(paste0("In case of stratification, ", 
                #             "sample sizes will respect ", 
                #             "proportions found in real data."))
                ##     h3("Dirichlet-Multinomial data generation", align = "center")
                # ) 
                # column(width = 12, 
                #     numericInput(inputId = "seed", 
                #         label = "Seed for the RNG", 
                #         value = 12345, min = 100, 
                #         max = round(.Machine$integer.max ^.5))
                # ),
                ),            # two divided fluidRows
                fluidRow(
                    hr(),
                    column(width = 12, 
                        h3("Monte Carlo Simulation Settings", align = "center")
                    ),
                    #column(width = 12, 
                    #    helpText(paste0("With these settings you can ",
                    #            "now run MonteCarlo simulations"))
                    #),
                    column(width = 3, 
                        numericInput(inputId = "MC", 
                            label = "# of Replications", 
                            value = 100, min = 50, 
                            max = 1e4)
                    ),
                    column(width = 3, 
                        numericInput(inputId = "alpha", 
                            label = "Significance", 
                            value = 0.1, min = 1e-3, step = 1e-3,
                            max = 0.25)
                    ),
                    column(width = 3, 
                        numericInput(inputId = "seed", 
                            label = "RNG Seed", 
                            value = 999, min = 100, 
                            max = round(.Machine$integer.max ^.5))
                    ),
                    column(width = 3, 
                        checkboxInput(inputId = "wmwTest", 
                            label = h4("Add WMW test (slow)"), 
                            value = FALSE)
                    ),
                    column(width = 12, align = "center", 
                        h4("Compute power with these sample sizes (option 1)")
                    ),
                    column(width = 4, 
                        numericInput(inputId = "n1", 
                            label = "  Sample Size Controls", 
                            value = 40, min = 3) 
                    ),
                    column(width = 4, 
                        numericInput(inputId = "n2", 
                            label = "  Sample Size Cases", 
                            value = 40, min = 3)
                    ),
                    column(width = 4, align = "center",
                        h6("."),
                        actionButton("mcStart", label = "Compute")
                    ),
                    column(width = 12, align = "center", 
                        h4("Or select a range of sample sizes (Option 2)")
                    ),
                    column(width = 9, align = "center", 
                        sliderInput(inputId = "sampleSizes",
                            label = "Select Range",
                            min = 5, max = 200, value = c(30, 70), step = 5
                        ), 
                        helpText(paste0("In this case controls and cases have the ",
                                "same sample size."))
                    ),
                    column(width = 3, align = "center",
                        h5("."),
                        actionButton("powerSimStart", label = "Compute"),
                        hr()
                    ), 
                    column(width = 12,
                        h3(paste0("Report Generation"), align = "center")
                    ), 
                    column(width = 9, 
                        h4(paste0(
                                "Click on \"Download\" to generate a PDF report with ",
                                "all inputs and outputs from the analysis."), 
                            align = "left"), 
                        helpText(paste0("Wait the end of one computation to click for ",
                            "another one, look at the loading bar at the top of the page."))
                    ), 
                    column(width = 3, 
                        downloadButton("downloadReport", label = "Download") 
                    ),
                    column(width = 12, hr())
                ),# END - fluidRow
                tags$footer(
#                     HTML(paste0(
#                         "<p>For details see the paper at ",
#                         "<a href=\"http://bioinformatics.oxfordjournals.org/",
#                             "content/early/2016/02/19/bioinformatics.btw099.abstract\">",
#                             " [Link to connected paper]</a>.</p>")))
                    HTML(paste0(
                        "<p>For details see the paper at ",
                        "<a href=\"http://doi.org/10.1093/bioinformatics/btw099\", ",
                        "target=\"_blank\">",
                        "[Link to connected Bioinformatics paper]</a>.</p>")))
            ),# END - sidebarPanel 
            mainPanel(
                tabsetPanel(id = "inTabs", 
                    tabPanel(title = "Abundance Curves",  #"Settings", 
                        fluidRow(
                            column(9, 
                                h3("Desired R.A.D.s", align = "center"),
                                # plotOutput("piPlot"),
                                # h3("One Generated R.A.D.s", align = "center"),
                                plotOutput("estPiPlot"),
                                plotOutput("estPiPlot2"),
                                plotOutput("estPiPlot3")
                            ), 
                            column(3, 
                                ##  p(img(src = "hmp_logo_NIH_retina.png", height = 200), 
                                ##      align = "center"),
                                # h3("Single Dataset Wald", align = "left"), 
                                ##  textOutput("hmpDiffTest")
                                # verbatimTextOutput("hmpDiffTest"),
                                ##  h3("One Generated R.A.D.s", align = "center"),
                                h3("Power Estimate (Option 1)"),
                                verbatimTextOutput("singleMcWaldResults") 
                            #                                h4(tableOutput("singleMcWaldResults")) 
                            )
                        )# END - fluidRow
                    ), # END - tabPanel: Settings
                    tabPanel(title = "Power vs. Sample Size",
                        fluidRow(
                            h3("Power vs. Sample Size (Option 2)", align = "left"),
                            plotOutput("resPowSimPlot", width = "100%", height = "550px",
                                click = "powPlotClick")#,
                        #column(8, 
                        #    h3("Rejection proportion", align = "left"),
                        #    verbatimTextOutput("mcHmpWaldResPrint")
                        #)
                        )# END - fluidRow
                    ), # END - tabPanel: Power vs. Sample Size
                    tabPanel(title = "How to use it",
                        fluidRow(
                            column(8, align = "center", 
                                h2("Side-by-Side Instructions", align = "center")),
                            column(10, align = "left", 
                                h3(paste0("Select the shape of the Relative Abundance ",
                                        "Curve for control samples  (\\(\\pi_1\\))"))
                            #                                tags$hr()
                            ),
                            column(3, align = "left", #offset = 1,
                                h4(paste0("Choose if you want to use stratification ", 
                                        "or not ('YES' or 'NO')"))
                            ),
                            column(3, align = "left",
                                h4("Number of most abundant OTUs to consider")
                            ),
                            column(4, align = "left",
                                h4(paste0("Overdispersion parameter can be set only if ",
                                        "\"Shape\" is one of the theoretical ones."))
                            ),
                            column(10, 
                                tags$hr(),
                                h3(paste0("OR upload your own file and estimate ",
                                        "\\(\\pi_1\\) from it"), 
                                    align = "left"),
                                h4(paste0(
                                        "Required file format: ",
                                        "each row is an OTU and each column a ", 
                                        "sample. Stratifying variable, if present, ",
                                        "has to be in the second row. ", 
                                        "Note: the application considers only the ", 
                                        "\"# of OTUs\" most abundant OTUs, so it has ",
                                        "to be changed after loading in case a ",
                                        "different number of OTUs must be considered."))
                            ),
                            column(6, align = "left", 
                                h4(paste0(
                                        "Specify before uploading: separator ",
                                        "character used in the file, ",
                                        "character defining quotes, ",
                                        "and presence/absence of stratifying variable."))
                            ),
                            column(4, 
                                h4(paste0(
                                        "Click on 'Remove File' to ignore the ",
                                        "uploaded file and go back to the selected ",
                                        "shape (only possible once)."))
                            ),
                            column(10, 
                                tags$hr(),
                                h3("Define the Relative Abundance Curve for the ",
                                    "Case samples (\\(\\pi_2\\))", 
                                    align = "left")
                            ),
                            #column(4, h3(".")),
                            column(4, 
                                h4(paste0(
                                        "How many OTUs are to be changed? Two ",
                                        "separate sets of OTUs can be specified, with ",
                                        "different relative changes. ",
                                        "Defaults are 2 and 2 respectively."))
                            ),
                            column(3, 
                                h4(paste0(
                                        "Are these OTUs the 'most abundant' or the ",
                                        "'least abundant' among the ones displayed?"))
                            ),
                            column(3,
                                h4(paste0(
                                        "How much should these OTUs be increased/decreased, ",
                                        "in percentage, with respect to the R.A.D. of ",
                                        "the control samples (null hypothesis)?"))
                            ),
                            column(10, 
                                helpText(paste0(
                                        "Remark: to ensure that \\(\\boldsymbol\\pi_2\\) ",
                                        "sums up to 1, (as it is a vector of proportions) ",
                                        "modified OTUs need to be compensated. The ", 
                                        "compensation strategy is described at the tab ",
                                        "\"Methodology\"."))
                            )
                        ),            # two divided fluidRows
                        fluidRow(
                            column(width = 10, 
                                hr(),
                                h3("Monte Carlo Simulation Settings", align = "left")
                            ),
                            column(width = 3, 
                                h4(paste0(
                                        "Number of Monte Carlo replications (simulated ",
                                        "datasets) on which power calculations are based."))
                            ),
                            column(width = 2, 
                                h4(paste0(
                                        "Significance level to be considered."
                                    ))
                            ),
                            column(width = 2, 
                                h4(paste0(
                                        "Seed for the Random Number Generator."
                                    ))
                            ),
                            column(width = 4, 
                                h4(paste0(
                                        "Add computation of Wilcoxon-Mann-Whitney test, ",
                                        "performed on each OTU separately.",
                                        "Benjamini-Hochberg multiplicity correction. ",
                                        "It takes much more time to finish."
                                    ))
                            ),
                            column(width = 10, align = "left", 
                                hr(),
                                h3("Option 1: with specific sample sizes.")
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
                                hr(),
                                h3(paste0(
                                        "Option 2: with a ",
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
                                hr(),
                                h3(paste0(
                                        "At the end you can generate a report ", 
                                        "containing all settings and results."))
                            ) 
                        )# END - fluidRow
                    ),# END - tabPanel: How to use it
                    tabPanel(title = "Methodology",
                        fluidRow(
                            column(8, align = "center", 
                                h2("Methodological Details", align = "center")),
                            column(10, align = "left", 
                                h3("The Statistical Model"),
                                h4("The Dirichlet-Multinomial distribution is employed ",
                                    "here to fit and to generate data. ", 
                                    "It is described by two parameters:"), 
                                h4("\\(\\theta\\): the overdispersion parameter, which ", 
                                    "measures the within-sample excess of variability ",
                                    "w.r.t. a multinomial distribution;"),
                                h4("\\(\\boldsymbol\\pi\\): the vector of relative ",
                                    "abundances."),
                                h4("We refer to La Rosa", em("et.al."), 
                                    " (2012) for details ",
                                HTML(paste0("<a target=\"_blank\", ",
                                    "href= \"http://journals.plos.org/plosone/article?id=10.1371/",
                                    "journal.pone.0052078\">(Direct link)</a>.")))
                            ), 
                            column(10, align = "left",
                                h3("Compensation Strategy"),
                                h4(
                                    "When OTUs of the \"cases\" group are modified through ",
                                    "the user interface, \\(\\boldsymbol\\pi_2\\) does ",
                                    "not sum to 1 anymore, therefore other elements of ",
                                    "this vector must be changed to compensate this. ", 
                                    "We opted for a simple strategy: the probability mass ",
                                    "to compensate is dissipated among all OTUs left ",
                                    "unchanged by the user, named \"unrequested\" OTUs."),
                                h4(strong("Algorithm Description:"), 
                                    "Let \\(\\boldsymbol\\pi_2^1\\) and \\(\\boldsymbol\\pi_2^0\\)",
                                    "be the vectors of relative abundances of requested and ",
                                    "unrequested OTUs, respectively (changed or not changed ",
                                    "by the user); let also \\(m_1\\) and \\(m_0\\) be ", 
                                    "the total sums of these vectors, and \\(d = m_0 + m_1 - 1\\)"), 
                                h4("In case \\(d > m_0\\) then all unrequested OTUs are ",
                                    "set to 0 and requested OTUs are rescaled by multiplying ",
                                    "each element of \\(\\boldsymbol\\pi_2^1\\) by ",
                                    "\\([m_1 - (d - m_0)] / m_1\\), this because it means the ",
                                    "probability mass to compensate is bigger than the total ",
                                    "mass available for compensation,", em("i.e. "), "the ",
                                    "requested effect size is too big ", em("w.r.t. "),
                                    "the chosen settings. ", br(),
                                    "In case \\(d \\leq m_0\\) then only \\(\\boldsymbol\\pi_2^0\\) ",
                                    "elements have to be rescaled by multiplying them by ",
                                    "\\((m_0 - d)/m_0\\).")
                            ),
                            column(10, align = "left",
                                h3("Stratified tests descriptions"),
                                h4(strong("Wald: "), 
                                    "Let \\(W_i\\) denote the Wald test statistic in stratum ",
                                    "\\(i=1,\\ldots, q\\). The asymptotic null distribution ",
                                    "of \\(W_i\\) is chi-squared with \\(d_i = K-1\\) ",
                                    "degrees of freedom. ",
                                    "Since sample observations from different strata ",
                                    "are mutually independent, the asymptotic null distribution ",
                                    "of the overall test statistic \\(\\sum_{i=1}^q W_i\\) ",
                                    "is chi-squared with \\(\\sum_{i=1}^q d_i\\) degrees of freedom."),
                                h4(strong("WMW: "), 
                                    "We adopted a similar strategy as for the Wald test. ",
                                    "If \\(S_i\\) denotes the standardised WMW statistic ",
                                    "then the asymptotic null distribution of \\(S_i\\) is ",
                                    "standard normal and hence the asymptotic null distribution ",
                                    "of \\(S_i^2\\) is chi-squared with one degree of freedom. ",
                                    "Since sample observations from different strata are ",
                                    "mutually independent, the overall test statistic given by ",
                                    "\\(\\sum_{i=1}^q S_i^2\\) asymptotically has a chi-squared ",
                                    "distribution with \\(q\\) degrees of freedom. By squaring ",
                                    "the \\(S_i\\) statistics before summing, the overall  test ",
                                    "does not rely on the assumption that the sign of the effect ",
                                    "size is the same in all strata. An overall test could also be based ",
                                    "\\(\\sum_{i=1}^q S_i\\) (similar to the construction of a ",
                                    "Cochran-Mantel-Haenszel test), but such a test may show small power ",
                                    "if the signs of the effect sizes are not the same in all strata."),
                                helpText(
                                    "Note that we do not include weights in the construction ",
                                    "of the overall test statistics. The reason is that the Shiny ",
                                    "applet is restricted to equal stratum sample sizes. Weights ",
                                    "(e.g. proportional to the stratum sample sizes) may be included ",
                                    "in the construction of the tests, but this will affect the ",
                                    "asymptotic null distributions.")
                            ),
                            column(10, align = "left",
                                h3("Power definitions"),
                                h4("Let \\(R\\) be the number of Monte Carlo data generations, ",
                                    "\\(K = K_0 + K_1\\) the number of OTUs considered among which ",
                                    "\\(K_0\\) are unrequested (not requested to be differentially ",
                                    "abundant between the two groups by the user) and ",
                                    "\\(K_1\\) are requested by the user (", em("i.e. "),
                                    "requested to be differentially abundant between the two groups), and ",
                                    "\\(\\alpha\\) the significance level."), 
                                #br(),
                                h4(strong("Wald: "), "power is simply given by ",
                                    "\\( \\# (\\mathbf p^W < \\alpha) / R \\), where ",
                                    "\\(\\mathbf p^W\\) is the vector of length \\(R\\) ",
                                    "containing the p-values of the test. It is the ",
                                    "average number of times the test was rejected among ",
                                    "the Monte Carlo replications."),
                                #br(),
                                h4(strong("WMW: "), "power is defined as ",
                                    "\\(w^{MW}=\\# [(\\sum_{k=1}^{K_1} p_k^{adj} < \\alpha) > 0] / R \\), ",
                                    "where \\( p_k^{adj}\\) is the (multiplicity adjusted) ",
                                    "p-value of the test for a requested OTU. It is the ",
                                    "average number of times ", em("at least "), "one ",
                                    "requested OTU was rejected by the test."),
                                #br(),
                                h4(strong("WMW avg: "), "power is defined individually ",
                                    "for each requested OTU and then averaged. In formula: ",
                                    "\\(w_k = \\# (p_k^{adj} < \\alpha) / R \\), \\(k=1, \\ldots, K_1\\), ",
                                    "\\(w_{avg}^{MW} = \\sum_{k=1}^{K_1} w_k / K_1\\). ",
                                    "It is the expected power for a generic requested OTU."),
                                #br(),
                                helpText("Note that often \"WMW avg\" will likely be smaller than ",
                                    "\"WMW\" as for the second, to be called \"rejected\", ",
                                    "is enough to have one rejection among requested OTUs.")
                            ),
                            column(10, align = "left",
                                h3("Enterotype stratification"),
                                h4("We refer to the following link for details about ",
                                    "the enterotype clustering algorithm."), 
                                HTML("<a target=\"_blank\", href = ",
                                    "\"http://enterotype.embl.de/enterotypes.html\">",
                                    "http://enterotype.embl.de/enterotypes.html\ </a>"),
                                h4("In a nutshell:"), 
                                h4(strong("(i)"), " build distance matrix between ", 
                                   "samples using Jensen-Shannon Divergence applied on ",
                                   "genus relative abundances;"), 
                                h4(strong("(ii)"), " compute optimal number ", 
                                   "of clusters with Calinski-Harabasz coefficient ",
                                   "(Calinski and Harabasz, 1974);"),
                                h4(strong("(iii)"), " calculate statistical significance of optimal ",
                                   "clustering with Silhouette coefficient (Rousseeuw, 1984).")
                            ),
                            column(10, align = "left",
                                h3("Examples"),
                                h4("We added two examples in the namesake tab (on the ",
                                    "right of Methodology), to compare stratified and ",
                                    "non-stratified analyses. The first example is with ",
                                    "Enterotype stratification, the second one with ",
                                    "stratification according to sex.")
                            )
                        )# END - fluidRow
                    ), # END - tabPanel: Methodology 
                    tabPanel(title = "Examples",
                        fluidRow(
                            column(9, align = "left",
                                h3(HTML(
                                    "<a href=\"#PooledEnterotype\">",
                                    "Pooled Enterotype Data</a>")),
                                h3(HTML(
                                    "<a href=\"#StrataEnterotype\">",
                                    "Enterotype Data Stratified by Enterotype</a>")),
                                h3(HTML(
                                    "<a href=\"#PooledStoolHMP\">",
                                    "Pooled HMP Stool Data</a>")),
                                h3(HTML(
                                    "<a href=\"#StrataStoolHMP\">",
                                    "HMP Stool Data Stratified by Sex</a>"))),
                                column(11, tags$hr())
                        ),
                        fluidRow(
                            column(9, align = "left", 
#                                 h1("Pooled or Stratified by Enterotype"),
                                HTML("<a name=\"PooledEnterotype\"></a>"),
                                includeHTML("./www/stoolData-noStrata.html"),
                                br(),
                                HTML("<a name=\"StrataEnterotype\"></a>"),
                                includeHTML("./www/stoolData-enteroStrata.html")),
                            column(11, tags$hr()), 
                            column(9, align = "left", 
#                                 h1("Pooled or Stratified by Sex"),
                                HTML("<a name=\"PooledStoolHMP\"></a>"),
                                includeHTML("./www/stoolHMP-noStrata.html"),
                                br(),
                                HTML("<a name=\"StrataStoolHMP\"></a>"),
                                includeHTML("./www/stoolHMP-sexStrata.html"))
                        )# END - fluidRow
                    )# END - tabPanel: Examples
                )# END - tabsetPanel
            )# END - mainPanel
        )# END - sidebarLayout
    )# END - fluidPage
)# END - shinyUI




### Notes: to implement
# count only rejections among the differentially abundant
# 1. individual OTU power averaged
# 2. min(adj.p.val) < alpha
# 3. TPR: avg(#rejections among DA / # true DA)





