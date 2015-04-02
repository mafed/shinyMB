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
###

#rm(list = ls()); gc(reset = TRUE)
#source("./shinyMB-0.2/www/auxCode.R")


source("./www/auxCode.R")


### manual input for testing purposes
#input <- list(entero = "YES", MC = 100,
#    kindPiOne = "Geometric", numOTUs = 100, 
#    mostLeastAb1 = "most abundant", mostLeastAb2 = "most abundant", 
#    diffOTUs1 = 2, diffOTUs2 = 3,
#    relAbund1 = 33, relAbund2 = 25, 
#    n1 = 30, n2 = 30, theta = 0.01, totCounts = 1e4, 
#    seed = 12345, alpha = .1,
#    sampleSizes = c(10, 40))


### 
shinyServer(
    function(input, output, session) {
      
      ### first RAD generation, OR estimated from uploaded dataset
      ## FIRST way of generating *piOne*
      piOne <- reactive({
            if (input$reset > 0 | is.null(input$userFiles))
            {
              if(input$entero == "NO")
              {
                piVec <- list(piOneGen(K = input$numOTUs, 
                        kind = switch(
                            input$kindPiOne, 
                            "Exponential"  = {"exp"}, 
                            "Geometric"    = {"geom"},
                            "Broken stick" = {"bstick"}, 
                            "Two-pieces Geometric" = {"brGeom"}
                        )
                    ))
                theta <- list(input$theta)
              } else
              {
                ent1 <- readRDS("./data/stoolMostAbEnt1_3.rds")
                ent2 <- readRDS("./data/stoolMostAbEnt2_3.rds")
                ent3 <- readRDS("./data/stoolMostAbEnt3_3.rds")
                piVec <- list(
                    "ent1" = msWaldHMP:::piMoM4Wald(ent1), 
                    "ent2" = msWaldHMP:::piMoM4Wald(ent2),
                    "ent3" = msWaldHMP:::piMoM4Wald(ent3))
                theta <- list(
                    "ent1" = msWaldHMP:::weirMoM4Wald(ent1), 
                    "ent2" = msWaldHMP:::weirMoM4Wald(ent2),
                    "ent3" = msWaldHMP:::weirMoM4Wald(ent3))
              }# END - ifelse: enterotype stratification
            } else
            {
              if(input$entero == "NO")
              {
                inFile <- input$userFiles
                #              if (length(inFile$datapath) > 1L)
                inData <- as.matrix(read.csv(inFile$datapath[1L], header = TRUE))
                rownames(inData) <- inData[, 1L]
                inData <- inData[, -1]
                ## standardisation if the dataset is in the format "obs x OTUs"
                piVec <- list(msWaldHMP:::piMoM4Wald(inData))
                theta <- list(msWaldHMP:::weirMoM4Wald(inData))
                
              } else
              {
                inFile <- input$userFiles
#                if (length(inFile$datapath) != 3)
                ## remove the first column as it usually contains rownames 
                inData1 <- as.matrix(read.csv(inFile$datapath[1L], header = TRUE))
                inData2 <- as.matrix(read.csv(inFile$datapath[2L], header = TRUE))
                inData3 <- as.matrix(read.csv(inFile$datapath[3L], header = TRUE))
                rownames(inData1) <- inData1[, 1L]
                rownames(inData2) <- inData2[, 1L]
                rownames(inData3) <- inData3[, 1L]
                inData1 <- inData1[, -1]
                inData2 <- inData2[, -1]
                inData3 <- inData3[, -1]
                
                ## standardisation if the dataset is in the format "obs x OTUs"
                piVec <- list(
                    "ent1" = msWaldHMP:::piMoM4Wald(inData1), 
                    "ent2" = msWaldHMP:::piMoM4Wald(inData2),
                    "ent3" = msWaldHMP:::piMoM4Wald(inData3))
                theta <- list(
                    "ent1" = msWaldHMP:::weirMoM4Wald(inData1), 
                    "ent2" = msWaldHMP:::weirMoM4Wald(inData2),
                    "ent3" = msWaldHMP:::weirMoM4Wald(inData3))
              }# END - ifelse: enterotype stratification
            }# END - ifelse: userFiles or generated data
            
#            piOne <- 
            list("piVec" = piVec, "theta" = theta)
          })
      
#      ## SECOND way for generating *piOne*
#      auxOneGen <- reactive({
#            piOneGen(K = length(piOne()$"piVec"), 
#                kind = switch(
#                    input$kindPiOne, 
#                    "Exponential"  = {"exp"}, 
#                    "Geometric"    = {"geom"},
#                    "Broken stick" = {"bstick"}, 
#                    "Two-pieces Geometric" = {"brGeom"}
#                )
#            )
#          })
#      
#      auxOneDat <- reactive({
#            ### look at functioning of "observe"
      ##              isolate({
#            inFile <- input$userFiles
#            inData <- read.csv(inFile$datapath, header = TRUE)[, -1]
#            multSampWald:::piMoM4Wald(inData)
      ##              })
#          })
#      
#      piOne <- eventReactive(input$reset, {
#            if (is.null(input$userFiles)) 
#            {
#              auxOneGen()
#            } else
#            {
#              auxOneDat()
#            }
#          }
#      )
      
      
      
      ### rest of the code
      changedOTUs <- reactive({
            otus2Change <- seq_len(input$diffOTUs1 + input$diffOTUs2)
            nOtus <- length(piOne()$"piVec"[[1L]])
            
            if(input$mostLeastAb1 == "least abund.")
            {
              otus2Change[seq_len(input$diffOTUs1)] <- 
                  (nOtus - input$diffOTUs1 - input$diffOTUs2 + 1):(
                    nOtus - input$diffOTUs2)
            } else {}
            
            if(input$mostLeastAb2 == "least abund.")
            {
              otus2Change[seq_len(input$diffOTUs2) + input$diffOTUs1] <- 
                  (nOtus - input$diffOTUs2 + 1):nOtus
            } else {}
            
            otus2Change
          })
      
      rho <- reactive({
            nOtus <- length(piOne()$"piVec"[[1L]])
            rhoGenSpecOtus(
                K = nOtus, 
                m1 = input$diffOTUs1, m2 = input$diffOTUs2, 
                relAbund1 = (input$relAbund1 + 100) / 100, 
                relAbund2 = (input$relAbund2 + 100) / 100, 
                otus2Change = changedOTUs())
          })
      
      
      piTwo <- reactive({
#            if (!is.list(piOne()$"piVec"))
#            {
#              piTwoGen(piOne()$"piVec", rho(), compensation = "relDiff")
#            } else
#            {
#              list(
#                  "ent1" = piTwoGen(piOne()$"piVec"[["ent1"]], rho()),
#                  "ent2" = piTwoGen(piOne()$"piVec"[["ent2"]], rho()),
#                  "ent3" = piTwoGen(piOne()$"piVec"[["ent3"]], rho()),
#              )
#            }# END - ifelse: enterotype stratification
            aux <- lapply(piOne()$"piVec", FUN = piTwoGen, rho = rho())
#            piTwo <-
            list(
                "piTwo" = lapply(aux, elNamed, name = "piTwo"),
                "otuType" = lapply(aux, elNamed, name = "otuType")
            )
          })
      
      
      output$piPlot <-  renderPlot({
            obj1 <- piOne()$"piVec"[[1L]]
            obj2 <- piTwo()$"piTwo"[[1L]]
            
            par(mar = c(4, 4, 1, 0) + .1)
            plot(obj1, pch = 1, col = 1, lwd = 2, 
                ylim = c(0, max(obj1, obj2)),
                xlab = "OTU Index", ylab = "Abundances", 
#                main = expression(paste("Ranked Abundances for ", pi[1])),
                cex.lab = 1.2, lty = 1:2)
            lines(obj1, lwd = 1, col = 1)
            points(obj2, pch = 3, lwd = 2, col = 2)
            lines(obj2, lwd = 1, col = 2, lty = 2)
            abline(h = 0, lty = 4, col = "gray60", lwd = 2)
            legend(x = "topright", pch = c(1, 3), col = c(1, 2), pt.lwd = 2,
                lwd = 2, lty = 1:2,
                legend = c(expression(pi[1]), 
                    expression(paste(pi[2], ""))
                ))
          })
      
      
      ### Dirichlet-Multinomial counts generation
      generatedCounts <- reactive({
#            countsGen(
#                sampleSizes = c(input$n1, input$n2), 
#                alphas = cbind(piOne()$"piVec"[[1L]], piTwo()[[1L]]$"piTwo"), 
#                theta = piOne()$"theta"[[1L]], 
#                K = length(piOne()$"piVec"[[1L]]), 
#                N = input$totCounts, 
#                seed = 12345)    # input$seed
            
            obj1 <- piOne()$"piVec"
            aux <- lapply(seq_along(obj1), FUN = function(i, obj1, obj2)
                {
                  countsGen(
                      sampleSizes = c(input$n1, input$n2),
                      alphas = cbind(obj1[[i]], obj2[[i]]),
                      theta = piOne()$"theta"[[i]],
                      K = length(obj1[[i]]),
                      N = input$totCounts,
                      seed = 12345 
                  )
                }, obj1 = obj1, obj2 = piTwo()$"piTwo")
            
#            generatedCounts <- 
            list(
                "dmDataList" = lapply(aux, elNamed, name = "dmDataList"),
                "piDirList" = lapply(aux, elNamed, name = "piDirList")
            )
          })
      
      ### plot for generated pi RADs
      output$estPiPlot <- renderPlot({
            main <- ifelse(test = input$entero == "YES", 
                yes = "Enterotype: 'Bacteroides'", # "Enterotype 1", 
                no = "")
#                no = "Ranked Abundance Distribution")
            
#            pdf(file = "plotRad.pdf", width = 8, height = 6)
#            svg(file = "plotRad.svg", width = 8, height = 6)
#            png(file = "plotRad.png", width = 1000, height = 750, pointsize = 24)
            drawPiPlot(
                countsData = generatedCounts()$piDirList[[1L]], 
                piOneObj = piOne()$"piVec"[[1L]],
                piTwoObj = piTwo()$"piTwo"[[1L]], main = main, 
                ylab = "Abundance Proportions", 
                theta = round(piOne()$"theta"[[1L]], 3L))
#            dev.off()
          })
      
      ## second enterotype plot
      output$estPiPlot2 <- renderPlot({
            if (input$entero == "YES")
            {
              drawPiPlot(
                  countsData = generatedCounts()$piDirList[[2L]], 
                  piOneObj = piOne()$"piVec"[[2L]],
                  piTwoObj = piTwo()$"piTwo"[[2L]], 
                  main = "Enterotype: 'Prevotella'", # "Enterotype 2", 
                  ylab = "Abundancy Proportions", 
                  theta = round(piOne()$"theta"[[2L]], 3L))
            } else {}
          })
      
      ## third enterotype plot
      output$estPiPlot3 <- renderPlot({
            if (input$entero == "YES")
            {
              drawPiPlot(
                  countsData = generatedCounts()$piDirList[[3L]], 
                  piOneObj = piOne()$"piVec"[[3L]],
                  piTwoObj = piTwo()$"piTwo"[[3L]], 
                  main = "Enterotype: 'Ruminococcus'", # "Enterotype 3", 
                  ylab = "Abundancy Proportions", 
                  theta = round(piOne()$"theta"[[3L]], 3L))
            } else {}
          })
      
      
      
      ### plot for generated counts data
      output$estCountsPlot <-  renderPlot({
            avgEstCountData <- sapply(generatedCounts()$dmDataList[[1]], FUN = colMeans)
            activeOTUs <- piTwo()$otuType[[1]] %in% 1:2
#            drawPiPlot(
#                countsData = generatedCounts($dmDataList[[1L]], 
#                piOneObj = piOne()$"piVec"[[1L]] * input$totCounts,
#                piTwoObj = piTwo()$"piTwo"[[1L]] * input$totCounts, 
#                main = "Counts", ylab = "Raw Abundances")
            plot(avgEstCountData[, 2] - avgEstCountData[, 1], type = "o", 
                xlab = "OTU index", ylab = expression(paste(bar(x[2]) - bar(x[1]))),
                main = "Avg diff. of generated Counts", cex = 1.5, pch = 20)
            abline(h = 0, lty = 4, col = "gray60", lwd = 2)
            points(avgEstCountData[activeOTUs, 2] - avgEstCountData[activeOTUs, 1], 
                pch = 4, col = "red", lwd = 2, cex = 2)
          })
      
      
      
      ### Single run differential abundance test
      auxHmpTest <- reactive({
            nOtus <- sapply(piOne()$"piVec", length)
            nReads <- lapply(generatedCounts()$"dmDataList"[[1L]], rowSums)
            thetaMoM <- sapply(generatedCounts()$"dmDataList", 
                FUN = function(x) 
                  sapply(x, FUN = msWaldHMP:::weirMoM4Wald)
            )
            piMoM <- lapply(generatedCounts()$"dmDataList", 
                FUN = function(x) sapply(x, FUN = msWaldHMP:::piMoM4Wald)
            )
            names(piMoM) <- names(piOne()$piVec)
            
            wald <- drop(msWald(nReads, alphaDM = piMoM, thetaDM = thetaMoM))
            
            DoFs <- (length(nReads) - 1) * (nOtus - 1)
            pVal <- pchisq(q = wald, df = DoFs, ncp = 0, lower.tail = FALSE)
            
            aux <- c(wald, pVal)
            
            tmp <- matrix(aux, nrow = 2, byrow = TRUE,
                dimnames = list(c("Test Stat.", "p.value"), NULL))
            
            
            ## global test, just sum up the Chi-squares values and DoFs
            if (length(wald) > 1L)
            {
              pValTot <- pchisq(q = sum(wald), df = length(wald) * DoFs[1L], 
                  lower.tail = FALSE)
              tmp <- cbind(tmp, c(sum(wald), pValTot))
              colnames(tmp) <- c(paste0("Ent.", seq_along(wald)), "Global")
#              c("Global Test" = sum(wald), "Global p.value" = pValTot)
              round(t(tmp), 3L)
              
            } else
            {
              round(tmp[, 1], 3L)
            }
          })
      
      output$hmpDiffTest <- renderPrint({
            auxHmpTest()
          })
      
      
      
      ### select tab with simulation results
      observe({
            input$powerSimStart
            updateTabsetPanel(session, inputId = "inTabs", selected = "Results")
          })
      
      ### select tab with settings by default
      observe({
            input$mcStart
            updateTabsetPanel(session, inputId = "inTabs", selected = "Settings")
          })
      
      
      ### power calculation among MC replications with settings already defined
      auxSingleWald <- reactive({
            input$mcStart
            
            set.seed(12345)
            isolate({
                  nSubsets <- 1L
                  nStrata <- length(piOne()$"piVec")
                  nGroups <- length(generatedCounts()$"dmDataList"[[1L]])
                  nOtus <- length(piOne()$"piVec"[[1L]])
                  nReads <- list(
                      rep.int(input$totCounts, input$n1), 
                      rep.int(input$totCounts, input$n2))
                  
                  ## create *alphaDM* for the function
                  alphaDM <- lapply(seq_len(nStrata), 
                      FUN = function(i, one, two)
                      {
                        cbind(one[[i]], two[[i]])
                      }, one = piOne()$piVec, two = piTwo()$piTwo)
                  names(alphaDM) <- names(piOne()$piVec)
                  
                  ## create *thetaDM* for the function
                  thetaDM <- matrix(rep.int(unlist(piOne()$"theta"), nGroups), 
                      nrow = nGroups, ncol = nStrata, byrow = TRUE, 
                      dimnames = list(paste0("group", seq_len(nGroups)), names(alphaDM)))
                  
                  ### quantiles of the reference distribution
                  ## Degrees of Freedom
                  DoF <- (nGroups - 1) * (nOtus - 1)
                  globDoF <- (nGroups - 1) * (nOtus * nStrata - 1)
                  ## for individual tests, simple Bonferroni correction
                  qAlpha <- qchisq(p = 1 - input$alpha/nStrata, df = DoF, 
                      lower.tail = TRUE)
                  qAlphaGlob <- qchisq(p = 1 - input$alpha, df = globDoF, 
                      lower.tail = TRUE)
                  
                  
#            res <- isolate({
                  ##                  res <- 
#                  myFunMC(MC = input$MC, 
#                      nReads =  nReads, 
#                      alphaDM = alphaDM, thetaDM = thetaDM, 
#                      sigLev = input$alpha)
#                })
                  
                  withProgress(message = "Computing...", value = 0, min = 0, max = 1,
                      expr = {
                        setProgress(0)
                        ### MonteCarlo simulations
                        tmp <- lapply(seq_len(input$MC), 
                            FUN = function(x)
                            {
                              setProgress(x/input$MC)
                              msWaldHMP:::msWald(
                                  nReads = nReads, alphaDM = alphaDM, thetaDM = thetaDM)
                            })
                        
                        res <- array(unlist(tmp), dim = c(nStrata, nSubsets, input$MC), 
                            dimnames = list(
                                rownames(el(tmp)), colnames(el(tmp)), 
                                paste0("MC", seq_len(input$MC))
                            ))
                      })# END - withProgress: end simulations
                  
                  rejRes <- rowSums(res > qAlpha, na.rm = TRUE, dims = 2) / input$MC
                  
                  ## if multilple strata, sum up together
                  if (nStrata > 1L)
                  {
                    rejResGlob <- rowSums(colSums(res, na.rm = TRUE) > qAlphaGlob) / input$MC
                    rejRes <- rbind(rejRes, "Global" = rejResGlob)
                  } else {}
                  
                  
                  
                  rejRes <- drop(rejRes)
                  
                  if (length(rejRes) == 1L)
                  {
                    names(rejRes) <- "Value"
                  } else {}
                })
            
            round(rejRes, 3L)
          })
          
      output$singleMcWaldResults <- renderPrint({
            auxSingleWald()
          })
      
      
      ### power calculation among MC replications with different sample sizes
      ## TODO: modify *myFunMC* s.t. takes in input a list of lists as *nReads* to 
      ## deal with multiple sample sizes on the same data (reduces computation time)
      mcHmpWaldResults <- reactive({
            ## Take a dependency on input$powerSimStart
            input$powerSimStart
            
            set.seed(12345)
            isolate({
                  nSubsets <- 5L
                  sampleSizes <- round(
                      seq(from = input$sampleSizes[1L], 
                          to = input$sampleSizes[2L], 
                          length = nSubsets))
                  nGroups <- length(generatedCounts()$"dmDataList"[[1L]])
                  nOtus <- length(piOne()$"piVec"[[1L]])
                  nStrata <- length(piOne()$"piVec")
                  
                  ## put the *nReads* object in the right format
                  tmpReads <- matrix(input$totCounts, 
                      nrow = max(sampleSizes), ncol = nGroups)
                  nReads <- vector("list", length(sampleSizes))
                  
                  for (sampRun in seq_along(sampleSizes))
                  {
                    nReads[[sampRun]] <- lapply(seq_len(nGroups),
                        FUN = function(i) 
                          tmpReads[seq_len(sampleSizes[sampRun]), i])
                  }
                  
                  ## create *alphaDM* for the function
                  alphaDM <- lapply(seq_along(piOne()$"piVec"), 
                      FUN = function(i, one, two)
                      {
                        cbind(one[[i]], two[[i]])
                      }, one = piOne()$piVec, two = piTwo()$piTwo)
                  names(alphaDM) <- names(piOne()$piVec)
                  
                  ## create *thetaDM* for the function
                  thetaDM <- matrix(rep.int(unlist(piOne()$"theta"), nGroups), 
                      nrow = nGroups, ncol = nStrata, byrow = TRUE, 
                      dimnames = list(paste0("group", seq_len(nGroups)), names(alphaDM)))
                  
                  ### quantiles of the reference distribution
                  ## Degrees of Freedom
                  DoF <- (nGroups - 1) * (nOtus - 1)
                  globDoF <- (nGroups - 1) * (nOtus * nStrata - 1)
                  ## for individual tests, simple Bonferroni correction
                  qAlpha <- qchisq(p = 1 - input$alpha/nStrata, df = DoF, 
                      lower.tail = TRUE)
                  qAlphaGlob <- qchisq(p = 1 - input$alpha, df = globDoF, 
                      lower.tail = TRUE)
                  
                  
                  ## NEW code for several sample sizes
                  withProgress(message = "Computing...", value = 0, min = 0, max = 1,
                      expr = {
                        setProgress(0)
#                  rejRes <- isolate({
#                       msWaldMC(MC = input$MC, nReads = nReads, alphaDM = alphaDM, 
#                            thetaDM = thetaDM)
#                      })
#                  setProgress(1)
                        ### MonteCarlo simulations
                        tmp <- lapply(seq_len(input$MC), 
                            FUN = function(x)
                            {
                              setProgress(x/input$MC)
                              msWaldHMP:::msWald(
                                  nReads = nReads, alphaDM = alphaDM, thetaDM = thetaDM)
                            })
                        
                        res <- array(unlist(tmp), dim = c(nStrata, nSubsets, input$MC), 
                            dimnames = list(
                                rownames(el(tmp)), colnames(el(tmp)), 
                                paste0("MC", seq_len(input$MC))
                            ))
                      })# END - withProgress: end simulations
                  
                  
                  rejRes <- rowSums(res > qAlpha, na.rm = TRUE, dims = 2) / input$MC
                  
                  ## if multilple strata, sum up together
                  if (nStrata > 1L)
                  {
                    rejResGlob <- rowSums(colSums(res, na.rm = TRUE) > qAlphaGlob) / input$MC
                    rejRes <- rbind(rejRes, "Global" = rejResGlob)
                  } else {}
                })
            
#            mcHmpWaldResults <- 
            list("pow" = tail(rejRes, n = 1), 
                "seqSizes" = sampleSizes)
          })
      
      
      ### print results
      output$mcHmpWaldResPrint <- renderPrint({
            if(input$powerSimStart > 0)
            {
              input$"powPlotClickCoord"
              round(mcHmpWaldResults()$"pow", 4L)
            } else
            {
              "Not yet started"
            }
          })
      
      ### plot results
      output$resPowSimPlot <- renderPlot({
            if(input$powerSimStart > 0)
            {
              par(mar = c(4, 4, 1, 1))
              plot(mcHmpWaldResults()$"seqSizes", mcHmpWaldResults()$"pow", 
                  xlim = input$sampleSizes + c(-2L, 2L), 
                  ylim = c(0, 1.1), pch = 19, lwd = 1,
                  type = "o", #main = "Power vs. Sample Size",
                  xlab = "sample size", ylab = "power")
              abline(h = c(0, input$alpha, 1), lty = 4, col = "gray70", lwd = 2)
              text(x =  min(input$sampleSizes), y = input$alpha, pos = 3, 
                  labels = paste0("alpha = ", input$alpha), cex = 1.3)
              
              if(!is.null(input$"powPlotClickCoord"))
              {
                coords <- input$"powPlotClickCoord"
                points(coords, pch = 3, lwd = 3, cex = 2.5, col = "red2")
                sampleVec <- seq(
                    from = input$sampleSizes[1L], to = input$sampleSizes[2L],
                    length = 100)
                findSpecificPowFun <- approxfun(
                    x = mcHmpWaldResults()$"seqSizes", 
                    y = mcHmpWaldResults()$"pow")
                
                resCoords <- c(
                    "x" = coords[["x"]], 
                    "y" = findSpecificPowFun(coords[["x"]])
                )
                lines(x = c(resCoords["x"], resCoords["x"]), 
                    y = c(0, resCoords["y"]), 
                    lty = 4, col = "red2", lwd = 2)
                lines(x = c(0, resCoords["x"]), 
                    y = c(resCoords["y"], resCoords["y"]), 
                    lty = 4, col = "red2", lwd = 2)
                
                text(x = coords[["x"]], y = coords[["y"]], 
#                text(x = resCoords["x"], y = resCoords["y"], 
                    labels = paste0(
                        "Size=", round(coords[["x"]], 2L),
#                        "Size=", round(resCoords["x"]),
                        "\n Power=", 
                        round(resCoords["y"], 3L)), 
                    pos = 3, offset = 1, cex = 1.5)
              } else {}
#            desPower <- 0.65
#            indOk <- max(which(findSpecificPowFun(sampleVec) <= desPower))
#            sampleVec[indOk]    # desired power needs at least this sample size
#            
            } else
            {
              par(mar = c(4, 4, 1, 1))
              plot(0, 0 , type = "n", xlim = input$sampleSizes, ylim = c(0, 1), 
                  xlab = "sample size", ylab = "power")
              text(mean(input$sampleSizes), .5, labels = "NOT YET INITIATED", 
                  cex = 2)
            }
            
          })
      
      
      ### download PDF report
      inputPars <- reactive({
            tmp <- c(
                "Enterotypes Stratification"     = input$entero, 
                "MonteCarlo Replications"           = input$MC, 
                "Pi One type"    = input$kindPiOne, 
                "Number of OTUs"      = input$numOTUs, 
                "Most/Least Abundant 1" = input$mostLeastAb1, 
                "Most/Least Abundant 2" = input$mostLeastAb2, 
                "Changed OTUs 1"    = input$diffOTUs1, 
                "Changed OTUs 2"    = input$diffOTUs2, 
                "Rel. Abund. Diff. 1"    = paste(input$relAbund1, "%"), 
                "Rel. Abund. Diff. 2"    = paste(input$relAbund2, "%"), 
                "Theta"        = input$theta, 
                "Total Reads"    = input$totCounts, 
                "Significance Level"       = input$alpha,
                "Sample 1 Size"           = input$n1, 
                "Sample 2 Size"           = input$n2, 
                "Min Sample Size"      = input$sampleSizes[1L], 
                "Min Sample Size"      = input$sampleSizes[2L], 
                " " = " ")
            aux <- matrix(" ", ncol = 4, nrow = length(tmp)/2)
            aux[, 1L] <- names(tmp)[seq_len(length(tmp)/2)]
            aux[, 2L] <- tmp[seq_len(length(tmp)/2)]
            aux[, 3L] <- names(tmp)[-seq_len(length(tmp)/2)]
            aux[, 4L] <- tmp[-seq_len(length(tmp)/2)]
            aux
          })
      
      output$downloadReport <- downloadHandler(
          filename = function() {
            'my-report.pdf'
          },
          
          content = function(file) {
            src <- normalizePath('template.Rmd')
            
            # temporarily switch to the temp dir, in case you do not have write
            # permission to the current working directory
            owd <- setwd(tempdir())
            on.exit(setwd(owd))
            file.copy(src, 'template.Rmd')
            
            require(rmarkdown)
            out <- render('template.Rmd', pdf_document())
            file.rename(out, file)
          }
      )
      
    })# END - shinyServer




#piOne <- piOneGen(kind = "geom")
#rho <- rhoGenSpecOtus(
#    K = 100, 
#    m1 = 3, m2 = 2, 
#    relAbund1 = 1.33, 
#    relAbund2 = 1.25, 
#    otus2Change = 1:5)
#piTwo <- piTwoGen(piOne, rho, compensation = "relDiff")$piTwo
#
#
#MC <- 1000
#Nrs <- list(rep.int(1e4, 30), rep.int(1e4, 30)) 
#piMat0 <- cbind(piOne, piOne)
#piMat1 <- cbind(piOne, piTwo)
#theta <- c(.02, .02) 
#alpha <- .1
#
#library(multSampWald)
#### under H0
#myFunMC(MC = MC, nReads = Nrs, alphaMat = piMat0, 
#    thetaVec = theta, sigLev = alpha)
#### under H1
#myFunMC(MC = MC, nReads = Nrs, alphaMat = piMat1, 
#    thetaVec = theta, sigLev = alpha)
#
#### under H0
#MC.Xmcupo.statistics(Nrs, MC, piOne, t(piMat0), theta, "ha", alpha)
#### under H1
#MC.Xmcupo.statistics(Nrs, MC, piOne, t(piMat1), theta, "ha", alpha)


