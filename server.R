################################################################################
### Title: server.R                                              ###
###                                                                          ###
### Project: MicrobiomeTry                                                ###
###                                                                          ###
### Version: 0.2 - 18/mar/2015                                    ###
###                                                                          ###
### Description: short description                                           ###
###                                                                          ###
### Authors: Federico Mattiello <Federico.Mattiello@UGent.be>                ###
###                                                                          ###
### Maintainer: Federico Mattiello <Federico.Mattiello@UGent.be>             ###
###                                                                          ###
### Versions:                                                                ###
###     > 0.2 - 18/mar/2015 - 14:36:33:                                      ###
###         creation                                                         ###
###                                                                          ###
################################################################################
####
#
#setwd("../")
### super slow!!!
#dat <- t(read.table("v13_psn_otu.genus.fixed.txt", header = TRUE, skip = 1L))
#str(dat)
#
#metaData <- read.table("v13_map_uniquebyPSN.txt", HMPbodysubsite)

#rm(list = ls()); gc(reset = TRUE)


#setwd("./shinyMB/")
source("./www/auxCode.R")


#### manual input for testing purposes
### H_0 to test FDR
#input <- list(strata = "YES", MC = 50, "wmwTest" = TRUE,
#    kindPiOne = "StoolEntero", numOTUs = 30, simulatedPiOne = FALSE,
#    mostLeastAb1 = "most abundant", mostLeastAb2 = "most abundant", 
#    diffOTUs1 = 2, diffOTUs2 = 3,
#    relAbund1 = 0, relAbund2 = 0, 
#    n1 = 30, n2 = 30, theta = 0.01, totCounts = 1e4, 
#    seed = 12345, alpha = .1,
#    sampleSizes = c(20, 60), userFiles = list(datapath = "./shinyMB/data/entTot.csv"))
#
### effect size at minimum
# input <- list(strata = "YES", MC = 100, "wmwTest" = TRUE,
# #input <- list(strata = "YES", MC = 50, "wmwTest" = TRUE,
#    kindPiOne = "StoolEntero", numOTUs = 50, simulatedPiOne = FALSE,
#    mostLeastAb1 = "most abundant", mostLeastAb2 = "most abundant", 
#    diffOTUs1 = 2, diffOTUs2 = 2,
#    relAbund1 = 15, relAbund2 = 20, 
#    n1 = 30, n2 = 30, theta = 0.01, totCounts = 1e4, 
#    seed = 12345, alpha = .1,
#    sampleSizes = c(20, 70), userFiles = list(datapath = "./shinyMB/data/entTot.csv"))
#
### effect size with over AND under abundance
#input <- list(strata = "YES", MC = 50, "wmwTest" = TRUE,
#    kindPiOne = "StoolEntero", numOTUs = 30, simulatedPiOne = FALSE,
#    mostLeastAb1 = "most abundant", mostLeastAb2 = "most abundant", 
#    diffOTUs1 = 2, diffOTUs2 = 3,
#    relAbund1 = -30, relAbund2 = 20, 
#    n1 = 30, n2 = 30, theta = 0.01, totCounts = 1e4, 
#    seed = 12345, alpha = .1,
#    sampleSizes = c(10, 40), userFiles = list(datapath = "./shinyMB/data/entTot.csv"))
#
### effect size like in the paper
# input <- list(strata = "YES", MC = 200, "wmwTest" = TRUE,
#     kindPiOne = "StoolEntero", numOTUs = 50, simulatedPiOne = FALSE,
#     mostLeastAb1 = "most abundant", mostLeastAb2 = "most abundant", 
#     diffOTUs1 = 2, diffOTUs2 = 2,
#     relAbund1 = 20, relAbund2 = 40, 
#     n1 = 40, n2 = 40, theta = 0.01, totCounts = 1e4, 
#     seed = 999, alpha = .1,
#     sampleSizes = c(30, 70), userFiles = list(datapath = "./shinyMB/data/stoolHMPsexStrata.csv"))
#
# input <- list(strata = "NO", MC = 1000, "wmwTest" = TRUE,
#    kindPiOne = "StoolEntero", numOTUs = 50, simulatedPiOne = FALSE,
#    mostLeastAb1 = "most abundant", mostLeastAb2 = "most abundant", 
#    diffOTUs1 = 4, diffOTUs2 = 2,
#    relAbund1 = 25, relAbund2 = 30, 
#    n1 = 30, n2 = 30, theta = 0.01, totCounts = 1e4, 
#    seed = 12345, alpha = .1,
#    sampleSizes = c(10, 40), userFiles = list(datapath = "./shinyMB/data/entTot.csv"))



### 
shinyServer(
    function(input, output, session) {
      
      ### first RAD generation, OR estimated from uploaded dataset
      ## FIRST way of generating *piOne*
      piOne <- reactive({
            ## check if  *piOne* needs to be simulated or not
            type <- input$kindPiOne
            simulatedPiOne <- type %in% c(
                "Geometric (Theoretical)", 
                "Stick-breaking (Theoretical)", 
                "Two-pieces Geometric (Theoretical)", 
                "Exponential (Theoretical)")
            nReadsFromData <- NULL
            
            if (input$reset > 0 | is.null(input$userFiles))
            {
              if(input$strata == "NO")
              {
                updateNumericInput(session, "kindPiOne", value = type)
                
                if (simulatedPiOne)
                {
                  maxNumOTU <- input$numOTUs
                  piVec <- list(piOneGen(K = input$numOTUs, 
                          kind = switch(
                              input$kindPiOne, 
                              "Exponential (Theoretical)"  = {"exp"}, 
                              "Geometric (Theoretical)"    = {"geom"},
                              "Stick-breaking (Theoretical)" = {"bstick"}, 
                              "Two-pieces Geometric (Theoretical)" = {"brGeom"}
                          )
                      ))
                  theta <- list(input$theta)
                } else
                {
                  if (type == "StoolEntero")
                  {
                    load("./data/stoolMostAbJointEnt.RData")
                    aux <- as.matrix(entTot[, -1L])
                    ord <- order(colSums(aux), decreasing = TRUE)
                    maxNumOTU <- min(length(ord), input$numOTUs)
                    aux <- aux[, ord[seq_len(maxNumOTU)]]
                    piVec <- list(msWaldHMP:::piMoM4Wald(aux))
                    theta <- list(msWaldHMP:::weirMoM4Wald(aux))
                  } else
                  {
                    load("./data/paramsDM.RData")
                    piVec <- piEstMoM[type]
#                     maxNumOTU <- min(length(piVec[[1L]]), input$numOTUs)
                    
                    maxNumOTU <- min(round(length(piVec[[1L]]) * .2), 1000L)
                    if (maxNumOTU < input$numOTUs)
                    {
                      input$numOTUs <- maxNumOTU
                      updateNumericInput(session, "numOTUs", value = maxNumOTU)
                    } else {}
                    
                    ### truncation of RAD curve
                    piVec[[type]] <- piVec[[type]][seq_len(maxNumOTU)]
                    piVec[[type]] <- piVec[[type]] / sum(piVec[[type]])
                    theta <- list(thEstMoM[type])
                  }# END - ifelse: Stool default
                }# END - ifelse: simulatedPiOne
              } else        # stratified _piOne_
              {
                type <- "StoolEntero"
                
                updateNumericInput(session, "kindPiOne", value = type)
                
                load("./data/stoolMostAbJointEnt.RData")
                strata <- as.factor(entTot[, "entero"])
                piVec <- vector("list", nlevels(strata))
                theta <- vector("list", nlevels(strata))
                names(piVec) <- names(theta) <- levels(strata)
                
                maxNumOTU <- ncol(entTot)
                if (maxNumOTU < input$numOTUs)
                {
                  input$numOTUs <- maxNumOTU
                  updateNumericInput(session, "numOTUs", value = maxNumOTU)
                } else {}
                
                for (stratumRun in levels(strata))
                {
                  aux <- as.matrix(entTot[entTot$"entero" == stratumRun, -1L])
                  ord <- order(colSums(aux), decreasing = TRUE)
                  aux <- aux[, ord[seq_len(maxNumOTU)]]
                  piVec[[stratumRun]] <- msWaldHMP:::piMoM4Wald(aux)
                  theta[[stratumRun]] <- msWaldHMP:::weirMoM4Wald(aux)
                }# END - for: strata
              }# END - ifelse: enterotype stratification
            } else      ## user uploaded file or after pressing *reset* button
            {
              ## set first column as rownames then remove it, first column must 
              ## contain observations IDs
              inFile <- input$userFiles
              
              inData <- read.table(inFile$datapath[1L], header = TRUE, row.names = 1,
                  sep = input$sep, quote = input$quote, check.names = FALSE)
              inData <- t(as.matrix(inData))
              
              if (input$strataRow)
              {
                strata <- as.factor(inData[, 1L])
                inData <- inData[, -1]
              } else {}
              
              dimNames <- dimnames(inData)
              inData <- apply(inData, 2, as.integer)
              dimnames(inData) <- dimNames
              
              if(input$strata == "NO")
              {
                ## standardisation if the dataset is in the format "obs x OTUs"
                ord <- order(colSums(inData), decreasing = TRUE)
                maxNumOTU <- length(ord)
                if (maxNumOTU < input$numOTUs)
                {
                  input$numOTUs <- maxNumOTU
                  updateNumericInput(session, "numOTUs", value = maxNumOTU)
                } else {}
                
                aux <- inData[, ord[seq_len(maxNumOTU)]]
                piVec <- list(msWaldHMP:::piMoM4Wald(aux))
                theta <- list(msWaldHMP:::weirMoM4Wald(aux))
                nReadsFromData <- list(rowSums(inData), rowSums(inData))
              } else
              {
                piVec <- vector("list", nlevels(strata))
                theta <- vector("list", nlevels(strata))
                names(piVec) <- names(theta) <- levels(strata)
                nReadsFromData <- piVec
                nReadsTmp <- rowSums(inData)
                
                ord <- order(colSums(inData), decreasing = TRUE)
                maxNumOTU <- length(ord)
                if (maxNumOTU < input$numOTUs)
                {
                  input$numOTUs <- maxNumOTU
                  updateNumericInput(session, "numOTUs", value = maxNumOTU)
                } else {}
                
                for (stratumRun in levels(strata))
                {
                  aux <- as.matrix(inData[strata == stratumRun, ])
                  ord <- order(colSums(aux), decreasing = TRUE)
                  aux <- aux[, ord[seq_len(maxNumOTU)]]
                  piVec[[stratumRun]] <- msWaldHMP:::piMoM4Wald(aux)
                  theta[[stratumRun]] <- msWaldHMP:::weirMoM4Wald(aux)
                  nReadsFromData[[stratumRun]] <- list(
                      nReadsTmp[strata == stratumRun], 
                      nReadsTmp[strata == stratumRun])
                }
                
              }# END - ifelse: enterotype stratification
            }# END - ifelse: userFiles or generated data
            
#            piOne <- 
            list("piOne" = piVec, "theta" = theta, "simulatedPiOne" = simulatedPiOne,
                "nReadsFromData" = nReadsFromData, "maxNumOTU" = maxNumOTU)
          })
      
      
      ### rest of the code
#      changedOTUs <- function() 1L:4
#      changedOTUs <- function() 1L:6
      changedOTUs <- reactive({
            otus2Change <- seq_len(input$diffOTUs1 + input$diffOTUs2)
            nOTUs <- length(piOne()$"piOne"[[1L]])
            
            if(input$mostLeastAb1 == "least abundant")
            {
              otus2Change[seq_len(input$diffOTUs1)] <- 
                  (nOTUs - input$diffOTUs1 - input$diffOTUs2 + 1):(
                    nOTUs - input$diffOTUs2)
            } else {}
            
            if(input$mostLeastAb2 == "least abundant")
            {
              otus2Change[seq_len(input$diffOTUs2) + input$diffOTUs1] <- 
                  (nOTUs - input$diffOTUs2 + 1):nOTUs
            } else {}
            
            otus2Change
          })
      
      rho <- reactive({
#            rho <- 
            rhoGenSpecOtus(
                K = input$numOTUs, 
                m1 = input$diffOTUs1, m2 = input$diffOTUs2, 
                relAbund1 = (input$relAbund1 + 100) / 100, 
                relAbund2 = (input$relAbund2 + 100) / 100, 
                otus2Change = changedOTUs())
          })
      
      
      piTwo <- reactive({
            tmp <- lapply(piOne()$"piOne", FUN = piTwoGenSimple, rho = rho())
#            tmp <- lapply(piOne$"piOne", FUN = piTwoGenSimple, rho = rho)
#            piTwo <- 
            list(
                "piTwo" = lapply(tmp, elNamed, name = "piTwo"),
                "otuType" = lapply(tmp, elNamed, name = "otuType"))
          })
      
      
      ### select tab with simulation results
      observe({
            input$powerSimStart
            updateTabsetPanel(session, inputId = "inTabs", 
                selected = "Power vs. Sample Size")
          })
      
      ### select tab with settings by default
      observe({
            input$mcStart
            updateTabsetPanel(session, inputId = "inTabs", 
                selected = "Abundance Curves")
          })
      
      
      ### plot for generated pi RADs
      piPlot1 <- reactive({
            main <- ifelse(test = input$strata == "YES", 
#                yes = names(piOne$"piOne")[1L], 
                yes = names(piOne()$"piOne")[1L], 
                no = "")
#                no = "Ranked Abundance Distribution")
#            pdf(file = "plotRadsStrata.pdf", width = 10, height = 6)
#            main <- names(piOne$piOne)
#            pdf(file = "plotRads.pdf", width = 10, height = 6)
#            main <- ""
#            jpeg(file = "plotRadsStrata.jpg", width = 800, height = 600, 
#                pointsize = 16, quality = 90)
#            png(file = "plotRad.png", width = 1000, height = 750, pointsize = 24)
#            layout(t(1:3))
#            pdf(file = "Bacteroides.pdf", width = 8, height = 6)
#            pdf(file = "Prevotella.pdf", width = 8, height = 6)
#            pdf(file = "Ruminococcus.pdf", width = 8, height = 6)
#            for (i in 1:3)
            
#             numOTUs2Plot <- min(piOne$maxNumOTU, input$numOTUs)
            numOTUs2Plot <- min(piOne()$maxNumOTU, input$numOTUs)
            drawPiPlot(
#                countsData = generatedCounts()$piDirList[[1L]], 
                countsData = NULL, 
                ylab = "Abundance Proportions", main = main, 
                piOneObj = piOne()$"piOne"[[1L]][seq_len(numOTUs2Plot)],
                piTwoObj = piTwo()$"piTwo"[[1L]][seq_len(numOTUs2Plot)], 
                theta = round(piOne()$"theta"[[1L]], 3L))
#                piOneObj = piOne$"piOne"[[i]],
#                piTwoObj = piTwo$"piTwo"[[i]], main = main[i], 
#                theta = round(piOne$"theta"[[i]], 3L))
#            dev.off()
          })
      
      output$estPiPlot <- renderPlot({
            piPlot1()
          })
      
      ## second enterotype plot
      piPlot2 <- reactive({
            if (input$strata == "YES" && length(piOne()$"piOne") > 1L)
            {
              numOTUs2Plot <- min(piOne()$maxNumOTU, input$numOTUs)
              drawPiPlot(
#                  countsData = generatedCounts()$piDirList[[2L]], 
                  countsData = NULL, 
                  piOneObj = piOne()$"piOne"[[2L]][seq_len(numOTUs2Plot)],
                  piTwoObj = piTwo()$"piTwo"[[2L]][seq_len(numOTUs2Plot)], 
#                  main = "Enterotype: 'Prevotella'", 
                  main = names(piOne()$"piOne")[2L],  
                  ylab = "Abundancy Proportions", 
                  theta = round(piOne()$"theta"[[2L]], 3L))
            } else {}
          })
      output$estPiPlot2 <- renderPlot({
            piPlot2()
          })
      
      ## third enterotype plot
      piPlot3 <- reactive({
            if (input$strata == "YES" && length(piOne()$"piOne") > 2L)
            {
              numOTUs2Plot <- min(piOne()$maxNumOTU, input$numOTUs)
              drawPiPlot(
#                  countsData = generatedCounts()$piDirList[[3L]], 
                  countsData = NULL, 
                  piOneObj = piOne()$"piOne"[[3L]][seq_len(numOTUs2Plot)],
                  piTwoObj = piTwo()$"piTwo"[[3L]][seq_len(numOTUs2Plot)], 
#                  main = "Enterotype: 'Ruminococcus'", 
                  main = names(piOne()$"piOne")[3L], 
                  ylab = "Abundancy Proportions", 
                  theta = round(piOne()$"theta"[[3L]], 3L))
            } else {}
          })
      
      output$estPiPlot3 <- renderPlot({
            piPlot3()
          })
      
      
      ### generation of all MC counts datasets
      totCountsGen <- reactive({
            
            ## take dependence only on the simulation start buttons
            input$powerSimStart
            input$mcStart
            
            isolate({
                set.seed(input$seed)
                invisible(a <- rnorm(1e3))
                type <- input$kindPiOne
                maxSampSize1 <- max(input$n1, input$sampleSizes)
                maxSampSize2 <- max(input$n2, input$sampleSizes)
                  
                  ## check if random library size or not
#                if (piOne()$"simulatedPiOne")
#                {
#                  nReads <- list(
#                      round(
#                          rgamma(maxSampSize1, shape = 3.5, rate = 7e-4), 
#                          digits = 0L), 
#                      round(
#                          rgamma(maxSampSize2, shape = 3.5, rate = 7e-4), 
#                          digits = 0L) 
#                  )
#                  tmpLibSize <- list(0)
#                  names(tmpLibSize) <- type
#                } else
#                {
#              if (is.null(piOne()$"nReadsFromData"))
#              {
#              } else
#              {
#                tmpLibSize <- piOne()$"nReadsFromData"
#              }# END - ifelse: library sizes taken from user-uploaded data
                if (piOne()$"simulatedPiOne")
                {
                  type <- "StoolEntero"
                } else {}
                
                load("./data/librarySizes.RData")
                
                if (type == "StoolEntero")
                {
                  load("./data/stoolMostAbJointEnt.RData")
                  entLibSizes <- rowSums(as.matrix(entTot[, -1]), na.rm = TRUE)
                  libSizesOrig$"StoolEntero" <- entLibSizes
                  libSizesOrigRaref$"StoolEntero" <- entLibSizes
                } else {}
                
                # XXXX
                tmpLibSize <- libSizesOrigRaref[[type]]
#                 inds <- tmpLibSize > input$numOTUs
#                 tmpLibSize <- libSizesOrig[[type]]
#                 inds <- tmpLibSize > piOne$maxNumOTU
                inds <- tmpLibSize > piOne()$maxNumOTU
                tmpLibSize <- tmpLibSize[inds]
                nReads <- list(
                    sample(tmpLibSize, size = maxSampSize1, replace = TRUE), 
                    sample(tmpLibSize, size = maxSampSize2, replace = TRUE))
#                }#END - ifelse: simulatedPiOne libSizes
                
                
                # nOTUs   <- length(piOne$"piOne"[[1L]])
                # nOTUs   <- length(piOne()$"piOne"[[1L]])
                # nOTUs <- min(input$numOTUs, piOne$maxNumOTU)
                # nStrata <- length(piOne$"piOne")
                # nGroups <- 2L
                nOTUs <- min(input$numOTUs, piOne()$maxNumOTU)
                nStrata <- length(piOne()$"piOne")
                nGroups <- 2L        # length(generatedCounts()$"dmDataList"[[1L]])
                
#                   obj1 <-  piOne$"piOne"
#                   obj2 <-  piTwo$"piTwo"
#                   theta <- piOne$"theta"
                obj1 <-  piOne()$"piOne"
                obj2 <-  piTwo()$"piTwo"
                theta <- piOne()$"theta"
                
                withProgress(message = "Generating data...", 
                    value = 0, min = 0, max = 1,
                    expr = {
                      setProgress(0)
                      totData <- lapply(seq_len(input$MC), 
                          FUN = function(x)
                          {
                            incProgress(round(1/input$MC, 3L))
                            tmp <- vector("list", length(obj1))
                            names(tmp) <- names(obj1)
                            for (stratumRun in seq_along(obj1))
                            {
                              aux <- countsGen(
#                              sampleSizes = c(input$n1, input$n2),
                                  sampleSizes = c(maxSampSize1, maxSampSize2),
                                  alphas = cbind(
                                      obj1[[stratumRun]], obj2[[stratumRun]]),
                                  theta = theta[[stratumRun]],
                                  K = length(obj1[[stratumRun]]),
#                                      N = input$totCounts,
                                  seed = NA, libSizes = nReads)[["dmDataList"]]
                              tmp[[stratumRun]] <- list(
                                  aux[[1L]][, seq_len(nOTUs)], 
                                  aux[[2L]][, seq_len(nOTUs)])
                            }# END - for: strata
                            # lapply(tmp, elNamed, name = "dmDataList")
                            tmp
                          }# END - function: lapply generation
                      )# END - lapply: totData
                    })# END - withProgress: end simulations
#            totCountsGen <- 
                c(totData, "nReads" = list(nReads))
              })
          })
      
      
      
      ### power calculation among MC replications with settings already defined
      auxSingleWald <- reactive({
            input$mcStart
            
            isolate({
                  dataList <- totCountsGen()
#                 dataList <- totCountsGen
                  wmwTest <- input$"wmwTest"
                  
#                   set.seed(12345)
                  set.seed(input$seed)
                  invisible(a <- rnorm(1e3))
                  
                  nSubsets <- 1L
                  n1 <- input$n1
                  n2 <- input$n2
                  mcReps <- length(dataList) - 1L
                  nStrata <- length(dataList[[1L]])
                  nGroups <- length(dataList[[c(1L, 1L)]])
                  nOTUs <- ncol(dataList[[c(1L, 1L, 1L)]])
#                  h1OTUs <- changedOTUs()
#                  h1OTUs <- which(el(piTwo$otuType) %in% 1L:2)
                  h1OTUs <- which(el(piTwo()$otuType) %in% 1L:2)
                  if (length(h1OTUs) == 0L)
                  {
                    h1OTUs <- seq_len(input$numOTUs)
                    minNumRej <- input$alpha * nOTUs
                  } else
                  {
                    minNumRej <- 0L
                  }
#                  nReads <- list(
#                      rep.int(input$totCounts, input$n1), 
#                      rep.int(input$totCounts, input$n2))
#                  nReads <- generatedCounts()$"libSizes"
#                  ## create *alphaDM* for the function
#                  alphaDM <- lapply(seq_len(nStrata), 
#                      FUN = function(i, one, two)
#                      {
#                        cbind(one[[i]], two[[i]])
#                      }, one = piOne()$"piOne", two = piTwo()$piTwo)
#                  names(alphaDM) <- names(piOne()$"piOne")
#                  
#                  ## create *thetaDM* for the function
#                  thetaDM <- matrix(rep.int(unlist(piOne()$"theta"), nGroups), 
#                      nrow = nGroups, ncol = nStrata, byrow = TRUE, 
#                      dimnames = list(paste0("group", seq_len(nGroups)), names(alphaDM)))
#                  thetaDM <- unlist(piOne()$"theta")
#                  nReads <- lapply(seq_along(nReads), FUN = function(i, nums) 
#                      {
#                        nReads[[i]][seq_len(nums[i])]
#                      }, nums = c(input$n1, input$n2))
#                  tmpReads <- totCountsGen$"nReads"
                  tmpReads <- dataList$"nReads"
                  nReads <- list(
                      tmpReads[[1L]][seq_len(n1)], 
                      tmpReads[[2L]][seq_len(n2)])
                  
#                  tmpWald <- rep.int(NA, nStrata)
#                  names(tmpWald) <- names(dataList[[1L]])
                  
                  
                  withProgress(message = "Computing...", 
                      value = 0, min = 0, max = 1,
                      expr = {
                        setProgress(0)
                        resTot <- sapply(seq_len(mcReps), 
                            FUN = function(mcRun)
                            {
                              incProgress(round(1/mcReps, 3L))
                              tmpWald <- rep.int(NA, nStrata)
                              for (strRun in seq_len(nStrata))
                              {
                                aux <- list(
                                    dataList[[c(mcRun, strRun, 1L)]][seq_len(n1), ], 
                                    dataList[[c(mcRun, strRun, 2L)]][seq_len(n2), ])
                                thAux <- sapply(aux, msWaldHMP:::weirMoM4Wald)
                                piAux <- sapply(aux, msWaldHMP:::piMoM4Wald)
                                tmpWald[strRun] <- 
                                    msWaldHMP:::msWaldStat(nReads, piAux, thAux)
                              }
                              tmpWald
                            })# END - lapply: resTot
                      })# END - withProgress
                  
                  resTot <- as.matrix(resTot)
                  
                  
                  ### in case Wilcoxon-Mann-Whitney is desired, 
                  ### if more than one stratum, use normal approximation
                  if (wmwTest)
                  {
                    pValsWMW <- array(NA, dim = c(nOTUs, nStrata, mcReps))
                    if (nStrata > 1L)
                    {
                      statWMW <- array(NA, dim = c(nOTUs, nStrata, mcReps))
                    } else {}
                    
                    withProgress(message = "Computing WMW...", 
                        value = 0, min = 0, max = 1,
                        expr = {
                          setProgress(0)
                          auxWMW <- sapply(seq_len(mcReps), 
                              FUN = function(mcRun)
                              {
                                incProgress(round(1/mcReps, 3L))
                                for (strRun in seq_len(nStrata))
                                {
                                  tmp <- msWaldHMP:::wrapperWMW(
                                      x = dataList[[c(mcRun, strRun, 1L)]][seq_len(n1), ],
                                      y = dataList[[c(mcRun, strRun, 2L)]][seq_len(n2), ],
                                      adjMethod = "fdr")
                                  pValsWMW[, strRun, mcRun] <<- tmp$p.value
                                  ## if more than one strata, obtain test-statistic
                                  if (nStrata > 1L)
                                  {
                                    statWMW[, strRun, mcRun] <<- tmp$statistic
                                  } else {}
                                }# END - for: strata
                              })# END - lapply: auxWMW
                        })# END - withProgress
                  } else {}# END - if: WMW test
                  
                  
                  
                  ### organise results
                  if (NCOL(resTot) > 1 && is.null(rownames(resTot)))
                  {
                    rownames(resTot) <- names(dataList[[1L]])
                  } else {}
                  
                  ### quantiles of the reference distribution
                  ## Degrees of Freedom
                  DoF <- (nGroups - 1) * (nOTUs - 1)
                  globDoF <- (nGroups - 1) * (nOTUs * nStrata - 1)
                  
                  pValTot <- pchisq(q = resTot, df = DoF, lower.tail = FALSE)
                  
                  ## multiplicity correction and sum up together strata
                  if (nStrata > 1L)
                  {
                    pValTot <- t(apply(pValTot, MARGIN = 2, FUN = p.adjust, 
                            method = "fdr"))
                    pValGlob <- pchisq(q = colSums(resTot), df = globDoF, 
                        lower.tail = FALSE)
                    powTotWald <- c(
                        colMeans(pValTot < input$alpha), 
                        "Global" = mean(pValGlob < input$alpha))
                    powTotWald <- as.matrix(powTotWald)
                    colnames(powTotWald) <- "Power"
                  } else
                  {
                    powTotWald <- colMeans(pValTot < input$alpha)
                  }
                  
                  if (length(powTotWald) == 1L)
                  {
                    names(powTotWald) <- "Value"
                  } else {}
                  
                  otuPowTable <- matrix(0, nrow = 2, ncol = 1)
                  
                  ## in case Wilcoxon-Mann-Whitney is desired
                  if (wmwTest)
                  {
                    ## select only Differentially Abundant OTUs
                    pValsWMWorig <- pValsWMW
                    pValsWMW <- pValsWMW[h1OTUs, , , drop = FALSE]
                    ## power for each truly DA OTU, averaged across all DA OTUs
                    singleOtuRej <- rowMeans(pValsWMW < input$alpha, na.rm = TRUE, 
                        dims = 2L)
                    powAvgWMW <- colMeans(singleOtuRej)
                    ## at least one rejection among DA OTUs, unless all are under H0
                    tmpRej <- colSums(pValsWMW < input$alpha, na.rm = TRUE, dims = 1L)
                    powTotWMW <- rowMeans(tmpRej > minNumRej)
                    
#              tmpOtuWMW <- apply(resWMW < input$alpha, MARGIN = c(1L, 3), sum, 
#                  na.rm = TRUE)
#              rejOtuWMW <- rowMeans(tmpOtuWMW > 0)
                    
                    ### data.frame containing power for individual OTUs, both truly DA and
                    ### used for compensation
#                     h1TotOTU <- which(el(piTwo$otuType) != 0)
                    h1TotOTU <- which(el(piTwo()$otuType) != 0)
                    allDaOtuRej <- rowMeans(
                        pValsWMWorig[h1TotOTU, , , drop = FALSE] < input$alpha, 
                        na.rm = TRUE, dims = 2L)
                    auxOtuPow <- vector("list", 2 * nStrata)
                    names(auxOtuPow) <- seq_len(2 * nStrata)
                    count <- 1L
                    for (stRun in seq_len(nStrata))
                    {
#                       tmp <- (piOne$piOne[[stRun]])[h1TotOTU]
                      tmp <- (piOne()$piOne[[stRun]])[h1TotOTU]
                      auxOtuPow[[count]] <- names(tmp)
                      names(tmp) <- NULL
                      auxOtuPow[[count + 1L]] <- round(tmp, 3L)
                      count <- count + 2L
                    }# END - for: arranging in a list
                    otuPowTable <- as.data.frame(auxOtuPow, stringsAsFactors = FALSE)
                    colnames(otuPowTable) <- c("Name", "Power")
                    
                    if (nStrata > 1L)
                    {
#                auxRej <- apply(pValsWMW, c(1L, 3L), p.adjust, method = "fdr")
#                globRej <- colSums(auxRej < input$alpha, na.rm = TRUE, dims = 2L)
#                powTotWMW <- c(powTotWMW, mean(globRej > 0))
#                globAvgRej <- colSums(auxRej < input$alpha, na.rm = TRUE, dims = 1L)
#                globAvgRej <- mean(rowMeans(globAvgRej > 0))
#                powAvgWMW <- c(powAvgWMW, globAvgRejjWMW
                      
                      ## select only Differentially Abundant OTUs
                      statWMW <- statWMW[h1OTUs, , , drop = FALSE]
                      ## normalise WMW test-statistics
                      ## (U - n1 n2/2) sqrt(n1 n2 (n1 + n2 + 1)/12
                      muU <-  n1 * n2/ 2
                      sig2U <- n1 * n2 * (n1 + n2 + 1)/12
                      statWMW <- (statWMW - muU)^2 / sig2U
                      ## sum together statistics for each stratum and compute p.values
                      auxStatWMW <- apply(statWMW, c(1L, 3L), sum)
                      auxPValWMW <- pchisq(q = auxStatWMW, df = nStrata, lower.tail = FALSE)
                      
                      ## power to reject *at least one* DA
                      globRej <- colSums(auxPValWMW < input$alpha, na.rm = TRUE, dims = 1L)
                      powTotWMW <- c(powTotWMW, mean(globRej > minNumRej))
                      ## average of individual DAs rejection rate
                      globAvgRej <- rowMeans(auxPValWMW < input$alpha, na.rm = TRUE)
                      powAvgWMW <- c(powAvgWMW, mean(globAvgRej))
                    } else {}
                    
                    totOut <- cbind(round(powTotWald, 3L), round(powTotWMW, 3L), 
                        round(powAvgWMW, 3L))
                    colnames(totOut) <- c("Wald", "WMW", "WMW avg")
                  } else
                  {
                    totOut <- round(powTotWald, 3L)
                  }# END - if: WMW test
                  
#                   auxSingleWald <- 
                  list(
                      "table2Display" = totOut, 
                      "otuPowTable" =  otuPowTable)
                })# END - isolate
          })# END - auxSingleWald
      
      output$singleMcWaldResults <- renderPrint({
#      output$singleMcWaldResults <- renderTable({
#            input$mcStart
#            which(el(piTwo()$otuType) %in% 1L:2)
            powTab <- auxSingleWald()$"table2Display"
            nStrata <- NROW(powTab)
            sampSizes <- matrix(NA_integer_, nrow = nStrata, ncol = 2L)
            colnames(sampSizes) <- c("n1", "n2")
            sampSizes[nStrata, ] <- c(as.integer(input$n1), as.integer(input$n2))
            powTab <- cbind("power" = powTab, sampSizes)
            powTab
          })
      
      
      ### power calculation among MC replications with different sample sizes
      ## TODO: modify *myFunMC* s.t. takes in input a list of lists as *nReads* to 
      ## deal with multiple sample sizes on the same data (reduces computation time)
      mcHmpWaldResults <- reactive({
            ## Take a dependency on input$powerSimStart
            input$powerSimStart
            
            isolate({
                  #set.seed(12345)
                  set.seed(input$seed)
                  invisible(a <- rnorm(1e3))
                  
                  dataList <- totCountsGen()
#                  dataList <- totCountsGen
                  nSubsets <- 7L
                  n1 <- input$n1
                  n2 <- input$n2
                  mcReps <- length(dataList) - 1L
                  nStrata <- length(dataList[[1L]])
                  nGroups <- length(dataList[[c(1L, 1L)]])
                  nOTUs <- ncol(dataList[[c(1L, 1L, 1L)]])
                  sampleSizes <- round(
                      seq(from = input$sampleSizes[1L], to = input$sampleSizes[2L], 
                          length = nSubsets))
                  
                  ### focus on the OTUs that are truly Differentially Abundant,
                  ### unless there is no true DA OTU, in that case rejecting
                  ### at most FDR = alpha * nOTUs is accetable because all discoveries
                  ### are false discoveries and alpha percent is still acceptable
#                  h1OTUs <- which(el(piTwo$otuType) %in% 1L:2)
                  h1OTUs <- which(el(piTwo()$otuType) %in% 1L:2)
                  if (length(h1OTUs) == 0L)
                  {
                    h1OTUs <- seq_len(nOTUs)
                    minNumRej <- input$alpha * nOTUs
                  } else
                  {
                    minNumRej <- 0L
                  }
                  
                  tmpReads <- dataList$"nReads"
                  nReads <- lapply(seq_along(sampleSizes), 
                      FUN = function(i)
                      {
                        list(
                            tmpReads[[1L]][seq_len(sampleSizes[i])], 
                            tmpReads[[2L]][seq_len(sampleSizes[i])]) 
                      })
                  names(nReads) <- sampleSizes
                  
                  resArr <- array(NA, dim = c(nSubsets, nStrata, mcReps), 
                      dimnames = list(sampleSizes, names(dataList[[1L]]), 
                          paste0("MC", seq_len(mcReps))))
                  
                  withProgress(message = "Computing...", 
                      value = 0, min = 0, max = 1,
                      expr = {
                        setProgress(0)
#                        resTot <- lapply(seq_len(mcReps), 
                        auxTot <- lapply(seq_len(mcReps), 
                            FUN = function(mcRun)
                            {
#                              setProgress(round(mcRun/mcReps, 2L))
                              incProgress(round(1/mcReps, 3L))
                              for (ssRun in seq_along(sampleSizes))
                              {
                                resArr[ssRun, , mcRun] <<- sapply(dataList[[mcRun]], 
                                    FUN = function(dat)
                                    {
                                      tmp <- list(
                                          dat[[1L]][seq_len(sampleSizes[ssRun]), ], 
                                          dat[[2L]][seq_len(sampleSizes[ssRun]), ])
                                      thAux <- sapply(tmp, msWaldHMP:::weirMoM4Wald)
                                      piAux <- sapply(tmp, msWaldHMP:::piMoM4Wald)
                                      msWaldHMP:::msWaldStat(nReads[[ssRun]], piAux, thAux)
                                    })
                              }# END - for: sampleSizes
                            })# END - lapply: mcRun
                      })# END - withProgress
                  
                  
                  ## Degrees of Freedom
                  DoF <- (nGroups - 1) * (nOTUs - 1)
                  globDoF <- (nGroups - 1) * (nOTUs * nStrata - 1)
                  pValTot <- pchisq(q = resArr, df = DoF, lower.tail = FALSE)
                  
                  ## multiplicity correction
                  if (nStrata > 1L)
                  {
                    pValTot <- apply(pValTot, MARGIN = c(1, 3), FUN = p.adjust, 
                        method = "fdr")
                    pValTot <- aperm(pValTot, c(3, 1, 2))
                    pValGlob <- pchisq(q = apply(resArr, c(1, 3), sum), df = globDoF, 
                        lower.tail = FALSE)
                    powTotWald <- rbind(
                        colMeans(pValTot < input$alpha),
                        "Global" = rowMeans(pValGlob < input$alpha))
                  } else
                  {
                    powTotWald <- rowMeans(drop(pValTot) < input$alpha)
                  }# END - ifelse: nStrata
                  
                  
                  ## in case Wilcoxon-Mann-Whitney is desired
                  if (input$"wmwTest")
                  {
                    pValArrWMW <- array(NA, 
                        dim = c(nOTUs, nStrata, nSubsets, mcReps), 
                        dimnames = list(
                            colnames(dataList[[c(1L, 1L, 1L)]]), 
                            names(dataList[[1L]]), 
                            sampleSizes, 
                            paste0("MC", seq_len(mcReps))))
                    if (nStrata > 1L)
                    {
                      statArrWMW <- pValArrWMW
                    } else {}
                    
                    withProgress(message = "Computing power WMW...", 
                        value = 0, min = 0, max = 1,
                        expr = {
                          setProgress(0)
#                        resTot <- lapply(seq_len(mcReps), 
                          auxTot <- lapply(seq_len(mcReps), 
                              FUN = function(mcRun)
                              {
#                              setProgress(round(mcRun/mcReps, 2L))
                                incProgress(round(1/mcReps, 3L))
                                for (ssRun in seq_along(sampleSizes))
                                {
                                  tmp <- lapply(dataList[[mcRun]], 
                                      FUN = function(dat)
                                      {
                                        msWaldHMP:::wrapperWMW(
                                            x = dat[[1L]][seq_len(sampleSizes[ssRun]), ],
                                            y = dat[[2L]][seq_len(sampleSizes[ssRun]), ],
                                            adjMethod = "fdr")
                                      })
                                  pValArrWMW[, , ssRun, mcRun] <<- 
                                      sapply(tmp, elNamed, name = "p.value")
                                  if (nStrata > 1L)
                                  {
                                    statArrWMW[, , ssRun, mcRun] <<- 
                                        sapply(tmp, elNamed, name = "statistic")
                                  } else {}#END - if stratified
                                }# END - for: sampleSizes
                              })# END - lapply: mcRun
                        })# END - withProgress
                    
                    
                    ## select only Differentially Abundant OTUs
                    pValArrWMW <- pValArrWMW[h1OTUs, , , , drop = FALSE]
                    ## power for each DA OTU, averaged across all DA OTUs
                    singleOtuRej <- rowMeans(pValArrWMW < input$alpha, na.rm = TRUE, 
                        dims = 3L)
                    powAvgWMW <- colMeans(singleOtuRej)
                    
                    ## at least one rejection among DA OTUs
                    tmpRej <- colSums(pValArrWMW < input$alpha, na.rm = TRUE, dims = 1L)
                    powTotWMW <- rowMeans(tmpRej > minNumRej, dims = 2L)
                    
                    if (nStrata > 1L)
                    {
#                      auxRej <- apply(pValArrWMW, c(1L, 3L, 4L), p.adjust, method = "fdr")
#                      globRej <- colSums(auxRej < input$alpha, na.rm = TRUE, dims = 2L)
#                      powTotWMW <- rbind(powTotWMW, "Global" = rowMeans(globRej > 0))
#                      
#                      globAvgRej <- colSums(auxRej < input$alpha, na.rm = TRUE, dims = 1L)
#                      globAvgRej <- rowMeans(globAvgRej > 0, dims = 2L)
#                      powAvgWMW <- rbind(powAvgWMW, "Global" = colMeans(globAvgRej))
#                      
#                      powOut <- cbind(
#                          powTotWald["Global", ], 
#                          powTotWMW["Global", ], powAvgWMW["Global", ])
#                      colnames(powOut) <- c("Wald", "WMW", "WMW avg")
                      
                      ## select only Differentially Abundant OTUs
                      statArrWMW <- statArrWMW[h1OTUs, , , , drop = FALSE]
                      
                      ## normalise WMW test-statistics, exploit array indices permutation
                      ## (U - n1 n2/2) sqrt(n1 n2 (n1 + n2 + 1)/12
                      muU <-  sampleSizes^2 / 2
                      sig2U <- sampleSizes^2 * (2 * sampleSizes + 1)/12
                      tmp <- aperm(statArrWMW, c(3, 1, 2, 4))
                      stdStatArrWMW <- (tmp - muU)^2 / sig2U
                      stdStatArrWMW <- aperm(stdStatArrWMW, c(2, 3, 1, 4))
                      
                      ## sum together statistics for each stratum and compute p.values
                      auxStatArrWMW <- apply(stdStatArrWMW, c(1L, 3L, 4L), sum, 
                          na.rm = TRUE)
                      auxPValArrWMW <- pchisq(q = auxStatArrWMW, df = nStrata, 
                          lower.tail = FALSE)
                      
                      ## power to reject *at least one* DA
                      globRej <- colSums(auxPValArrWMW < input$alpha, na.rm = TRUE, 
                          dims = 1L)
                      powTotWMW <- rbind(powTotWMW, "Global" = rowMeans(globRej > minNumRej))
                      ## average of individual DAs rejection rate
                      globAvgRej <- rowMeans(auxPValArrWMW < input$alpha, na.rm = TRUE, 
                          dims = 2L)
                      powAvgWMW <- rbind(powAvgWMW, "Global" = colMeans(globAvgRej))
                      
                      powOut <- cbind(
                          powTotWald["Global", ], 
                          powTotWMW["Global", ], powAvgWMW["Global", ])
                    } else
                    {
                      powOut <- cbind(
                          powTotWald, drop(powTotWMW), drop(powAvgWMW))
                    }# END - ifelse: more than one stratum
                    
                    colnames(powOut) <- c("Wald", "WMW", "WMW avg")
                    
                  } else        # *** only Wald test ***
                  {
                    if (nStrata > 1L)
                    {
                      powTmp <- t(powTotWald["Global", , drop = FALSE])
                    } else
                    {
                      powTmp <- as.matrix(powTotWald)
                    }# END - ifelse: more than one stratum
                    
                    powOut <- matrix(NA, nrow = NROW(powTmp), ncol = 3L, 
                        dimnames = list(
                            rownames(powTmp), 
                            c("Wald", "WMW", "WMW avg")))
                    powOut[, "Wald"] <- powTmp
                  }# END - ifelse: only Wald or also WMW test?
                  
#                    save.image(file = "debugSessionWMW.RData")
#                    load("debugSessionWMW.RData")
#                  ## create *alphaDM* for the function
#                  alphaDM <- lapply(seq_along(piOne()$"piOne"), 
#                      FUN = function(i, one, two)
#                      {
#                        cbind(one[[i]], two[[i]])
#                      }, one = piOne()$"piOne", two = piTwo()$piTwo)
#                  names(alphaDM) <- names(piOne()$"piOne")
#                  
#                  ## create *thetaDM* for the function
#                  thetaDM <- matrix(rep.int(unlist(piOne()$"theta"), nGroups), 
#                      nrow = nGroups, ncol = nStrata, byrow = TRUE, 
#                      dimnames = list(paste0("group", seq_len(nGroups)), names(alphaDM)))
                })# END - isolate
            
#            mcHmpWaldResults <- 
            list("seqSizes" = sampleSizes, 
                "pow" = powOut)
          })# END - reactive: mcHmpWaldResults
      
      
      ### print results
      output$mcHmpWaldResPrint <- renderPrint({
            if(input$powerSimStart > 0)
            {
#              input$"powPlotClick"
              round(mcHmpWaldResults()$"pow", 4L)
            } else
            {
              "Not yet started"
            }
          })
      
      
      powPlotCoords <- reactiveValues(x = NULL, y = NULL)
      observe({
            if (is.null(input$"powPlotClick"))
            {
              return()
            } else {}
            ### show selected point
#            isolate({
            powPlotCoords$x <- input$"powPlotClick"$x
            powPlotCoords$y <- input$"powPlotClick"$y
#                })
          })
      
      
      ### plot results
      powSimPlot <- reactive({
            input$powerSimStart
            input$"powPlotClick"
            
#            seqSizes <- mcHmpWaldResults$"seqSizes"
#            powerData <- mcHmpWaldResults$"pow"
            seqSizes <- mcHmpWaldResults()$"seqSizes"
            powerData <- mcHmpWaldResults()$"pow"
            
            ### interpolate points with natural splines
            ## interpolated Wald data
            waldLine <- spline(
                x = seqSizes, 
                y = powerData[, "Wald"], method = "natural")
            waldLine$y[waldLine$y > 1] <- 1
            
            
            ### in case "wmwTest" is selected from the interface but 
            ### simulations are not yet performed
            notYetComputedWMW <- all(is.na(powerData[, "WMW"]))
            
            if (input$"wmwTest" && !notYetComputedWMW)
            {
              ## interpolated WMW data, total power: at least one DA detected
              wmwLineTot <- spline(
                  x = seqSizes, 
                  y = powerData[, "WMW"], method = "natural")
              wmwLineTot$y[wmwLineTot$y > 1] <- 1
              ## interpolated WMW data, single OTU power average across DA OTUs
              wmwLineAvg <- spline(
                  x = seqSizes, 
                  y = powerData[, "WMW avg"], method = "natural")
              wmwLineAvg$y[wmwLineAvg$y > 1] <- 1
            } else {}
            
#            sampleVec <- seq(
#                from = input$sampleSizes[1L], to = input$sampleSizes[2L],
#                length = 100)
#                findSpecificPowFun <- approxfun(
#                    x = seqSizes, 
#                    y = powerData)
            
            par(mar = c(4, 4, 1, 1))
            if(input$powerSimStart == 0)
            {
              plot(0, 0 , type = "n", xlim = input$sampleSizes, ylim = c(0, 1), 
                  xlab = "sample size", ylab = "power")
              text(mean(input$sampleSizes), .5, labels = "NOT YET STARTED", 
                  cex = 2)
            } else {}
            
            plot(seqSizes, powerData[, "Wald"], col = "black", bg = "red2",
                xlim = input$sampleSizes + c(-2L, 2L), 
                ylim = c(0, 1.1), pch = 21, lwd = 2, 
                type = "p", #main = "Power vs. Sample Size",
                xlab = "sample size", ylab = "power", cex = 2)
            lines(waldLine, col = "red2", lty = 1, lwd = 2)
            
            if (input$"wmwTest" && !notYetComputedWMW)
            {
              points(seqSizes, powerData[, "WMW"], 
                  col = "black", bg = "blue2", pch = 22, cex = 2, lwd = 2)
              lines(wmwLineTot, col = "blue2", lty = 1, lwd = 2)
              points(seqSizes, powerData[, "WMW avg"], 
                  col = "black", bg = "blue2", pch = 23, cex = 2, lwd = 2)
              lines(wmwLineAvg, col = "blue2", lty = 2, lwd = 2)
              legend(x = "left", pch = 21:23, 
                  pt.bg = c("red2", "blue2", "blue2"), col = "black", 
                  lty = c(1, 1, 2), lwd = 2, pt.cex = 1.5, 
                  legend = c("Wald", "WMW", "WMW avg"))
            } else {}
            
            
            abline(h = c(0, input$alpha, 1), lty = 4, col = "gray70", lwd = 2)
            text(x = .8 * diff(input$sampleSizes), y = input$alpha, pos = 3, 
                labels = paste0("alpha = ", input$alpha), cex = 1.3)
          })
      
#      ## saving plot on file in a temporary directory 
#      savePdfPlot <- reactive({
#            fileName <- paste0(tempdir(), .Platform[["file.sep"]], "copyPowPlot.pdf")
#            pdf(file = fileName)
#            powSimPlot()
#            dev.off()
#            fileName
#          })
      
      powPlotPoint <- reactive({
            input$"powPlotClick"
#            dev.copy2pdf(file = "copyPowPlot.pdf")
            powSimPlot()
#            seqSizes <- mcHmpWaldResults$"seqSizes"
#            powerData <- mcHmpWaldResults$"pow"
            seqSizes <- mcHmpWaldResults()$"seqSizes"
            powerData <- mcHmpWaldResults()$"pow"
            
            ### in case "wmwTest" is selected from the interface but 
            ### simulations are not yet performed
            notYetComputedWMW <- all(is.na(powerData[, "WMW"]))
            
            ## function creating horizontal or vertical line
            lineCoords <- function(dat, direction = c("horizontal", "vertical"))
            {
              if (direction == "horizontal")
              {
                list("x" = c(0, dat["x"]), "y" = rep.int(dat["y"], 2L))
              } else
              {
                list("x" = rep.int(dat["x"], 2L), "y" = c(0, dat["y"]))
              }
            }# END: function - lineCoords
            
            ## function that interpolates data
            specificPowFunWald <- splinefun(
                x = seqSizes, 
                y = powerData[, 1L])
            
            if (input$"wmwTest" && !notYetComputedWMW)
            {
              specificPowFunWMW <- splinefun(
                  x = seqSizes, 
                  y = powerData[, 2L])
            } else {}
            
            
            
            if (!is.null(powPlotCoords$x))
            {
              coordsPointWald <- c(
                  "x" = powPlotCoords$x, 
                  "y" = min(1, specificPowFunWald(powPlotCoords$x)))
              
              if (input$"wmwTest" && !notYetComputedWMW)
              {
                coordsPointWMW <- c(
                    "x" = powPlotCoords$x, 
                    "y" = min(1, specificPowFunWMW(powPlotCoords$x)))
              } else {}
              
              ## cross corresponding to click
              points(powPlotCoords, pch = 3, lwd = 3, cex = 2.5, col = "black")
              
              ### Wald power line
              points(as.list(coordsPointWald), pch = 3, lwd = 3, cex = 1.5, col = "red2")
              ## horizontal line
              lines(lineCoords(coordsPointWald, direction = "horizontal"), 
                  lty = 4, col = "red2", lwd = 2)
              ## vertical line
#              lines(x = c(resCoordsWald["x"], resCoordsWald["x"]), 
#                  y = c(0, resCoordsWald["y"]), 
#                  lty = 4, col = "red2", lwd = 2)
              lines(lineCoords(coordsPointWald, direction = "vertical"), 
                  lty = 4, col = "red2", lwd = 2)
              ## text showing coordinates values
#              text(x = powPlotCoords$x, y = powPlotCoords$y, 
              text(as.list(coordsPointWald), 
                  labels = paste0(
                      "Size=", round(powPlotCoords$x, 2L),
                      "\n Power=", 
                      round(coordsPointWald["y"], 3L)), 
                  pos = 3, offset = 1, cex = 1.5)
              
              ### Wilcoxon-Mann-Whitney power line
              if (input$"wmwTest" && !notYetComputedWMW)
              {
                points(as.list(coordsPointWMW), pch = 3, lwd = 3, cex = 1.5, col = "blue2")
                ## horizontal line
                lines(lineCoords(coordsPointWMW, direction = "horizontal"), 
                    lty = 4, col = "blue2", lwd = 2)
                ## vertical line
                lines(lineCoords(coordsPointWMW, direction = "vertical"), 
                    lty = 4, col = "blue2", lwd = 2)
                ## text showing coordinates values
#                text(x = powPlotCoords$x, y = powPlotCoords$y, 
                text(as.list(coordsPointWMW), 
                    labels = paste0(
                        "Size=", round(powPlotCoords$x, 2L),
                        "\n Power=", 
                        round(coordsPointWMW["y"], 3L)), 
                    pos = 3, offset = 1, cex = 1.5)
              }# END: if - WMW test
            } else {}# END: if - clicked
          })
      
      output$resPowSimPlot <- renderPlot({
#            validate(
#                need(!is.null(input$"powPlotHover"), "not yet clicked")
#            )
#            powSimPlot()
            input$powPlotClick
            input$powerSimStart
            powPlotPoint()
          })
      
      
      ### download PDF report
#      inputPars <- reactive({
#            tmp <- c(
#                "Enterotypes Stratification"     = input$strata, 
#                "MonteCarlo Replications"           = input$MC, 
#                "Pi One type"    = input$kindPiOne, 
#                "Number of OTUs"      = input$numOTUs, 
#                "Most/Least Abundant 1" = input$mostLeastAb1, 
#                "Most/Least Abundant 2" = input$mostLeastAb2, 
#                "Changed OTUs 1"    = input$diffOTUs1, 
#                "Changed OTUs 2"    = input$diffOTUs2, 
#                "Rel. Abund. Diff. 1"    = paste(input$relAbund1, "%"), 
#                "Rel. Abund. Diff. 2"    = paste(input$relAbund2, "%"), 
#                "Theta"        = input$theta, 
#                "Total Reads"    = input$totCounts, 
#                "Significance Level"       = input$alpha,
#                "Sample 1 Size"           = input$n1, 
#                "Sample 2 Size"           = input$n2, 
#                "Min Sample Size"      = input$sampleSizes[1L], 
#                "Max Sample Size"      = input$sampleSizes[2L], 
#                " " = " ")
#            aux <- matrix(" ", ncol = 4, nrow = length(tmp)/2)
#            aux[, 1L] <- names(tmp)[seq_len(length(tmp)/2)]
#            aux[, 2L] <- tmp[seq_len(length(tmp)/2)]
#            aux[, 3L] <- names(tmp)[-seq_len(length(tmp)/2)]
#            aux[, 4L] <- tmp[-seq_len(length(tmp)/2)]
#            aux
#          })
      
      
      outputPars <- reactive({
            
          })
      
      output$downloadReport <- downloadHandler(
          filename = function() {
            aux <- format(Sys.time(), "%Y-%m-%d_%X")
            for (i in 1L:3)
              aux <- sub(pattern = ":", replacement = ".", x = aux, fixed = TRUE)
            paste0('ReportMB_', aux, '.pdf')
          },
          
          content = function(file) {
            src <- normalizePath('./template.Rmd')
            
            # temporarily switch to the temp dir, in case you do not have write
            # permission to the current working directory
            owd <- setwd(tempdir())
            on.exit(setwd(owd))
            file.copy(src, 'template.Rmd')
            
#            require(rmarkdown)
            out <- render('template.Rmd', pdf_document(highlight = "haddock"))
#             out <- render('template.Rmd', NULL)
            file.rename(out, file)
          }
      )
      
    })# END - shinyServer

### rarefication strategy based on order statistics (max 10 obs removed)
# sortLibSizes <- lapply(libSizesOrig, function(x) log(sort(x, decreasing = FALSE)))
# diffLibSizes <- lapply(sortLibSizes, diff)
# sapply(diffLibSizes, function(x) which.max(x[1L:10]))
# newLibSizeRaref <- vector("list", length(sortLibSizes))
# names(newLibSizeRaref) <- names(sortLibSizes)
# aux <- lapply(seq_along(sortLibSizes), 
#     function(i)
#     {
#       seqSubjRmv <- seq_len(which.max(diffLibSizes[[i]][1L:10]))
#       tmpNames <- names(sortLibSizes[[i]][-seqSubjRmv])
#       newLibSizeRaref[[i]] <<- libSizesOrig[[i]][tmpNames]
#     })
# sapply(newLibSizeRaref, range)
# 
# plot(sortLibSizes$Subgingival_plaque)
# plot(diffLibSizes$Right_Antecubital_fossa)

