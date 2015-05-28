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
#source("./shinyMB/www/auxCode.R")


source("./www/auxCode.R")


### manual input for testing purposes
#rm(list = ls()); gc(reset = TRUE)
# input <- list(strata = "YES", MC = 100,
#    kindPiOne = "Stool", numOTUs = 50, 
#    mostLeastAb1 = "most abundant", mostLeastAb2 = "most abundant", 
#    diffOTUs1 = 2, diffOTUs2 = 3,
#    relAbund1 = 33, relAbund2 = 25, 
#    n1 = 41, n2 = 45, theta = 0.01, totCounts = 1e4, 
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
                if (simulatedPiOne)
                {
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
                  if (type == "Stool")
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
                    maxNumOTU <- min(length(piVec[[1L]]), input$numOTUs)
                    piVec[[type]] <- piVec[[type]][seq_len(maxNumOTU)]
                    piVec[[type]] <- piVec[[type]] / sum(piVec[[type]])
                    theta <- list(thEstMoM[type])
                  }# END - ifelse: Stool default
                }# END - ifelse: simulatedPiOne
              } else
              {
                load("./data/stoolMostAbJointEnt.RData")
                strata <- as.factor(entTot[, "entero"])
                piVec <- vector("list", nlevels(strata))
                theta <- vector("list", nlevels(strata))
                names(piVec) <- names(theta) <- levels(strata)
                
                for (stratumRun in levels(strata))
                {
                  aux <- as.matrix(entTot[entTot$"entero" == stratumRun, -1L])
                  ord <- order(colSums(aux), decreasing = TRUE)
                  maxNumOTU <- min(length(ord), input$numOTUs)
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
              inData <- read.csv(inFile$datapath[1L], header = TRUE)
              rownames(inData) <- inData[, 1L]
              inData <- inData[, -1]
              
              if(input$strata == "NO")
              {
                ## standardisation if the dataset is in the format "obs x OTUs"
                aux <- as.matrix(inData[, -1L])
                ord <- order(colSums(aux), decreasing = TRUE)
                maxNumOTU <- min(length(ord), input$numOTUs)
                aux <- aux[, ord[seq_len(maxNumOTU)]]
                piVec <- list(msWaldHMP:::piMoM4Wald(aux))
                theta <- list(msWaldHMP:::weirMoM4Wald(aux))
                nReadsFromData <- list(rowSums(inData), rowSums(inData))
              } else
              {
                strata <- as.factor(inData[, 1L])
                piVec <- vector("list", nlevels(strata))
                theta <- vector("list", nlevels(strata))
                names(piVec) <- names(theta) <- levels(strata)
                nReadsFromData <- piVec
                inData <- inData[, -1]
                nReadsTmp <- rowSums(inData)
                
                for (stratumRun in levels(strata))
                {
                  aux <- as.matrix(inData[strata == stratumRun, ])
                  ord <- order(colSums(aux), decreasing = TRUE)
                  maxNumOTU <- min(length(ord), input$numOTUs)
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
                "nReadsFromData" = nReadsFromData)
          })
      
      
      ### rest of the code
#      changedOTUs <- 1L:5
      changedOTUs <- reactive({
            otus2Change <- seq_len(input$diffOTUs1 + input$diffOTUs2)
            nOtus <- length(piOne()$"piOne"[[1L]])
            
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
            rhoGenSpecOtus(
                K = input$numOTUs, 
                m1 = input$diffOTUs1, m2 = input$diffOTUs2, 
                relAbund1 = (input$relAbund1 + 100) / 100, 
                relAbund2 = (input$relAbund2 + 100) / 100, 
                otus2Change = changedOTUs())
          })
      
      
      piTwo <- reactive({
#            tmp <- lapply(piOne$"piOne", FUN = piTwoGen, rho = rho)
            tmp <- lapply(piOne()$"piOne", FUN = piTwoGen, rho = rho())
#            piTwo <-
            list(
                "piTwo" = lapply(tmp, elNamed, name = "piTwo"),
                "otuType" = lapply(tmp, elNamed, name = "otuType"))
          })
      
      
      output$piPlot <-  renderPlot({
            obj1 <- piOne()$"piOne"[[1L]]
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
                selected = "Settings")
          })
      
      
      ### Dirichlet-Multinomial counts generation
      generatedCounts <- reactive({
#            countsGen(
#                sampleSizes = c(input$n1, input$n2), 
#                alphas = cbind(piOne()$"piOne"[[1L]], piTwo()[[1L]]$"piTwo"), 
#                theta = piOne()$"theta"[[1L]], 
#                K = length(piOne()$"piOne"[[1L]]), 
#                N = input$totCounts, 
#                seed = 12345)    # input$seed
            
#            obj1 <- piOne$"piOne"
            obj1 <- piOne()$"piOne"
            type <- input$kindPiOne
            
            ## check if random library size or not
            if (piOne()$"simulatedPiOne")
            {
              libSizes <- list(
                  round(
                      rgamma(input$n1, shape = 3.5, rate = 7e-4), 
                      digits = 0L), 
                  round(
                      rgamma(input$n2, shape = 3.5, rate = 7e-4), 
                      digits = 0L) 
              )
              tmpLibSize <- list(0)
              names(tmpLibSize) <- type
            } else
            {
              load("./data/librarySizes.RData")
              tmpLibSize <- libSizesOrigRaref[[type]]
              inds <- tmpLibSize > 2000
              tmpLibSize <- tmpLibSize[inds]
              libSizes <- list(
                  sample(tmpLibSize, size = input$n1, replace = TRUE), 
                  sample(tmpLibSize, size = input$n2, replace = TRUE))
            }#END - ifelse: libSizes
            
            aux <- lapply(seq_along(obj1), FUN = function(i, obj1, obj2)
                {
                  countsGen(
                      sampleSizes = c(input$n1, input$n2),
                      alphas = cbind(obj1[[i]], obj2[[i]]),
#                      theta = piOne$"theta"[[i]],
                      theta = piOne()$"theta"[[i]],
                      K = length(obj1[[i]]),
                      N = input$totCounts,
                      seed = 12345, libSizes = libSizes)
#                }, obj1 = obj1, obj2 = piTwo$"piTwo")
                }, obj1 = obj1, obj2 = piTwo()$"piTwo")
            
#            generatedCounts <- 
            list(
                "dmDataList" = lapply(aux, elNamed, name = "dmDataList"),
                "piDirList" = lapply(aux, elNamed, name = "piDirList"),
                "libSizes" = list(libSizes),
                "libSizesOrig" = tmpLibSize)
          })
      
      ### plot for generated pi RADs
      piPlot1 <- reactive({
            main <- ifelse(test = input$strata == "YES", 
                yes = names(piOne()$"piOne")[1L], 
                no = "")
#                no = "Ranked Abundance Distribution")
            
#            pdf(file = "plotRad.pdf", width = 8, height = 6)
#            svg(file = "plotRad.svg", width = 8, height = 6)
#            png(file = "plotRad.png", width = 1000, height = 750, pointsize = 24)
            drawPiPlot(
                countsData = generatedCounts()$piDirList[[1L]], 
                piOneObj = piOne()$"piOne"[[1L]],
                piTwoObj = piTwo()$"piTwo"[[1L]], main = main, 
                ylab = "Abundance Proportions", 
                theta = round(piOne()$"theta"[[1L]], 3L))
#            dev.off()
            
          })
      output$estPiPlot <- renderPlot({
            piPlot1()
          })
      
      ## second enterotype plot
      piPlot2 <- reactive({
            if (input$strata == "YES" && length(piOne()$"piOne") > 1L)
            {
              drawPiPlot(
                  countsData = generatedCounts()$piDirList[[2L]], 
                  piOneObj = piOne()$"piOne"[[2L]],
                  piTwoObj = piTwo()$"piTwo"[[2L]], 
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
              drawPiPlot(
                  countsData = generatedCounts()$piDirList[[3L]], 
                  piOneObj = piOne()$"piOne"[[3L]],
                  piTwoObj = piTwo()$"piTwo"[[3L]], 
#                  main = "Enterotype: 'Ruminococcus'", 
                  main = names(piOne()$"piOne")[3L], 
                  ylab = "Abundancy Proportions", 
                  theta = round(piOne()$"theta"[[3L]], 3L))
            } else {}
          })
      output$estPiPlot3 <- renderPlot({
            piPlot3()
          })
      
      
      
      ### Single run differential abundance test
      auxHmpTest <- reactive({
#            nOtus <- sapply(piOne$"piOne", length)
            nOtus <- sapply(piOne()$"piOne", length)
#            nReads <- generatedCounts$"libSizes"
            nReads <- generatedCounts()$"libSizes"
            thetaMoM <- sapply(generatedCounts()$"dmDataList", 
                FUN = function(x) sapply(x, FUN = msWaldHMP:::weirMoM4Wald))
            piMoM <- lapply(generatedCounts()$"dmDataList", 
                FUN = function(x) sapply(x, FUN = msWaldHMP:::piMoM4Wald))
            names(piMoM) <- names(piOne()$"piOne")
            
            wald <- drop(msWaldHMP:::msWald(nReads, alphaDM = piMoM, thetaDM = thetaMoM))
            
            DoFs <- (length(el(nReads)) - 1) * (nOtus[1L] - 1)
            pVal <- pchisq(q = wald, df = DoFs, ncp = 0, lower.tail = FALSE)
            
            aux <- c(wald, pVal)
            
            tmp <- matrix(aux, nrow = 2, byrow = TRUE,
                dimnames = list(c("Test Stat.", "p.value"), names(wald)))
            
            
            ## global test, just sum up the Chi-squares values and DoFs
            if (length(wald) > 1L)
            {
              pValTot <- pchisq(q = sum(wald), df = length(wald) * DoFs[1L], 
                  lower.tail = FALSE)
              tmp <- cbind(tmp, "Global" = c(sum(wald), pValTot))
#              colnames(tmp) <- c(paste0("Ent.", seq_along(wald)), "Global")
#              c("Global Test" = sum(wald), "Global p.value" = pValTot)
              round(t(tmp), 4L)
              
            } else
            {
              round(tmp[, 1], 3L)
            }
          })
      
      output$hmpDiffTest <- renderPrint({
            auxHmpTest()
          })
      
      
      totCountsGen <- reactive({
            set.seed(12345)
#            input$mcStart
#            input$powerSimStart
            input$n1
            input$n1
            input$sampleSizes
            type <- input$kindPiOne
            maxSampSize1 <- max(input$n1, input$sampleSizes)
            maxSampSize2 <- max(input$n2, input$sampleSizes)
            
            ## check if random library size or not
            if (piOne()$"simulatedPiOne")
            {
              nReads <- list(
                  round(
                      rgamma(maxSampSize1, shape = 3.5, rate = 7e-4), 
                      digits = 0L), 
                  round(
                      rgamma(maxSampSize2, shape = 3.5, rate = 7e-4), 
                      digits = 0L) 
              )
              tmpLibSize <- list(0)
              names(tmpLibSize) <- type
            } else
            {
#              if (is.null(piOne()$"nReadsFromData"))
#              {
#              } else
#              {
#                tmpLibSize <- piOne()$"nReadsFromData"
#              }# END - ifelse: library sizes taken from user-uploaded data
              load("./data/librarySizes.RData")
              tmpLibSize <- libSizesOrigRaref[[type]]
              inds <- tmpLibSize > 2000
              tmpLibSize <- tmpLibSize[inds]
              nReads <- list(
                  sample(tmpLibSize, size = maxSampSize1, replace = TRUE), 
                  sample(tmpLibSize, size = maxSampSize2, replace = TRUE))
            }#END - ifelse: libSizes
            
            
#            nOtus   <- length(piOne$"piOne"[[1L]])
#            nStrata <- length(piOne$"piOne")
#            nGroups <- length(generatedCounts$"dmDataList"[[1L]])
            nOtus   <- length(piOne()$"piOne"[[1L]])
            nStrata <- length(piOne()$"piOne")
            nGroups <- 2L        # length(generatedCounts()$"dmDataList"[[1L]])
            
#            obj1 <-  piOne$"piOne"
#            obj2 <-  piTwo$"piTwo"
#            theta <- piOne$"theta"
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
#                        setProgress(round(x/input$MC, 2L))
#                        incProgress(round(1/input$MC, 2L))
                        tmp <- vector("list", length(obj1))
                        names(tmp) <- names(obj1)
                        for (stratumRun in seq_along(obj1))
                        {
                          tmp[[stratumRun]] <- countsGen(
#                              sampleSizes = c(input$n1, input$n2),
                              sampleSizes = c(maxSampSize1, maxSampSize2),
                              alphas = cbind(
                                  obj1[[stratumRun]], obj2[[stratumRun]]),
                              theta = theta[[stratumRun]],
                              K = nOtus,
#                                      N = input$totCounts,
                              seed = NA, libSizes = nReads)
                        }# END - for: strata
                        lapply(tmp, elNamed, name = "dmDataList")
                      }# END - function: lapply generation
                  )# END - lapply: totData
                })# END - withProgress: end simulations
#            totCountsGen <- 
            c(totData, "nReads" = list(nReads))
          })
      
      ### power calculation among MC replications with settings already defined
      auxSingleWald <- reactive({
            input$mcStart
            set.seed(12345)
            
            isolate({
                  nSubsets <- 1L
                  dataList <- totCountsGen()
                  nStrata <- length(dataList[[1L]])
                  nGroups <- length(dataList[[c(1L, 1L)]])
                  nOtus <- ncol(dataList[[c(1L, 1L, 1L)]])
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
                      tmpReads[[1L]][seq_len(input$n1)], 
                      tmpReads[[2L]][seq_len(input$n2)])
                  
                  tmpWald <- rep.int(NA, nStrata)
                  names(tmpWald) <- names(dataList[[1L]]) 
                  
                  withProgress(message = "Computing...", 
                      value = 0, min = 0, max = 1,
                      expr = {
                        setProgress(0)
                        resTot <- sapply(seq_len(input$MC), 
                            FUN = function(mcRun)
                            {
#                              setProgress(round(mcRun/input$MC, 2L))
                              incProgress(round(1/input$MC, 2L))
                              tmpWald[] <- NA
                              for (strRun in seq_len(nStrata))
                              {
                                aux <- list(
                                    dataList[[c(mcRun, strRun, 1L)]][seq_len(input$n1), ], 
                                    dataList[[c(mcRun, strRun, 2L)]][seq_len(input$n2), ])
                                thAux <- sapply(aux, msWaldHMP:::weirMoM4Wald)
                                piAux <- sapply(aux, msWaldHMP:::piMoM4Wald)
                                tmpWald[strRun] <- 
                                    msWaldHMP:::msWaldStat(nReads, piAux, thAux)
                              }
                              tmpWald
                            })# END - lapply: resTot
                      })# END - withProgress
                  
                  resTot <- as.matrix(resTot)
                  if (ncol(resTot) > 1 && is.null(rownames(resTot)))
                  {
                    rownames(resTot) <- names(dataList[[1L]])
                  } else {}
                  
                  ### quantiles of the reference distribution
                  ## Degrees of Freedom
                  DoF <- (nGroups - 1) * (nOtus - 1)
                  globDoF <- (nGroups - 1) * (nOtus * nStrata - 1)
                  
                  pValTot <- pchisq(q = resTot, df = DoF, lower.tail = FALSE)
                  
                  ## multiplicity correction and sum up together strata
                  if (nStrata > 1L)
                  {
                    pValTot <- t(apply(pValTot, MARGIN = 2, FUN = p.adjust, 
                            method = "holm"))
                    pValGlob <- pchisq(q = colSums(resTot), df = globDoF, 
                        lower.tail = FALSE)
                    rejTot <- c(
                        colMeans(pValTot <= input$alpha), 
                        "Global" = mean(pValGlob <= input$alpha))
                    rejTot <- as.matrix(rejTot)
                    colnames(rejTot) <- "Power"
                  } else
                  {
                    rejTot <- colMeans(pValTot <= input$alpha)
                  }
                  
                  if (length(rejTot) == 1L)
                  {
                    names(rejTot) <- "Value"
                  } else {}
                })# END - isolate
            
            round(rejTot, 3L)
          })# END - auxSingleWald
      
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
                  nSubsets <- 7L
                  dataList <- totCountsGen()
                  nStrata <- length(dataList[[1L]])
                  nGroups <- length(dataList[[c(1L, 1L)]])
                  nOtus <- ncol(dataList[[c(1L, 1L, 1L)]])
                  sampleSizes <- round(
                      seq(from = input$sampleSizes[1L], to = input$sampleSizes[2L], 
                          length = nSubsets))
                  
                  ## check if random library size or not
                  if (piOne()$"simulatedPiOne")
                  {
                    nReads <- lapply(seq_along(sampleSizes), FUN = function(i)
                        {
                          list(
                              round(
                                  rgamma(sampleSizes[i], shape = 3.5, rate = 7e-4), 
                                  digits = 0L), 
                              round(
                                  rgamma(sampleSizes[i], shape = 3.5, rate = 7e-4), 
                                  digits = 0L)
                          )
                        })
                  } else
                  {
                    tmpReads <- dataList$"nReads"
                    nReads <- lapply(seq_along(sampleSizes), 
                        FUN = function(i)
                        {
#                          list(
#                              sample(generatedCounts()$"libSizesOrig", 
#                                  size = sampleSizes[i], 
#                                  replace = TRUE), 
#                              sample(generatedCounts()$"libSizesOrig", 
#                                  size = sampleSizes[i], 
#                                  replace = TRUE))
                          list(
                              tmpReads[[1L]][seq_len(sampleSizes[i])], 
                              tmpReads[[2L]][seq_len(sampleSizes[i])]) 
                        })
                    names(nReads) <- sampleSizes
                  }#END - ifelse: simulatedPiOne, libSizes
                  
                  tmpWald <- matrix(NA, nrow = nSubsets, ncol = nStrata)
                  rownames(tmpWald) <- sampleSizes
                  colnames(tmpWald) <- names(dataList[[1L]])
                  
                  withProgress(message = "Computing...", 
                      value = 0, min = 0, max = 1,
                      expr = {
                        setProgress(0)
                        resTot <- lapply(seq_len(input$MC), 
                            FUN = function(mcRun)
                            {
#                              setProgress(round(mcRun/input$MC, 2L))
                              incProgress(round(1/input$MC, 2L))
                              for (ssRun in seq_along(sampleSizes))
                              {
                                tmpWald[ssRun, ] <- sapply(dataList[[mcRun]], 
                                    FUN = function(dat)
                                    {
                                      tmp <- list(
                                          dat[[1L]][seq_len(sampleSizes[ssRun]), ], 
                                          dat[[2L]][seq_len(sampleSizes[ssRun]), ])
                                      thAux <- sapply(tmp, msWaldHMP:::weirMoM4Wald)
                                      piAux <- sapply(tmp, msWaldHMP:::piMoM4Wald)
                                      msWaldHMP:::msWaldStat(nReads[[ssRun]], piAux, thAux)
                                    })
                              }# END - for: strata
                              tmpWald
                            })# END - lapply: mcRun
                      })# END - withProgress
                  
                  resArr <- array(NA, dim = c(nSubsets, nStrata, input$MC), 
                      dimnames = list(rownames(resTot[[1L]]), colnames(resTot[[1L]]), 
                          paste0("MC", seq_len(input$MC))))
                  for (mcRun in seq_along(resTot))
                  {
                    resArr[, , mcRun] <- resTot[[mcRun]]
                  }
                  
                  ## Degrees of Freedom
                  DoF <- (nGroups - 1) * (nOtus - 1)
                  globDoF <- (nGroups - 1) * (nOtus * nStrata - 1)
                  pValTot <- pchisq(q = resArr, df = DoF, lower.tail = FALSE)
                  
                  ## multiplicity correction
                  if (nStrata > 1L)
                  {
                    pValTot <- apply(pValTot, MARGIN = c(1, 3), FUN = p.adjust, 
                        method = "holm")
                    pValTot <- aperm(pValTot, c(3, 1, 2))
                    pValGlob <- pchisq(q = apply(resArr, c(1, 3), sum), df = globDoF, 
                        lower.tail = FALSE)
                    rejTot <- rbind(
                        colMeans(pValTot <= input$alpha),
                        "Global" = rowMeans(pValGlob <= input$alpha))
                  } else
                  {
                    rejTot <- rowMeans(drop(pValTot) <= input$alpha)
                  }# END - ifelse: nStrata
                  
                  
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
                "pow" = if(nStrata > 1L)
                    {
                      tail(rejTot, n = 1)
                    } else
                    {
                      rejTot
                    })
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
            
            ### interpolate points with natural splines
            ## interpolated data
            lineData <- spline(
                x = mcHmpWaldResults()$"seqSizes", 
                y = mcHmpWaldResults()$"pow", method = "natural")
            lineData$y[lineData$y > 1] <- 1
#            sampleVec <- seq(
#                from = input$sampleSizes[1L], to = input$sampleSizes[2L],
#                length = 100)
#                findSpecificPowFun <- approxfun(
#                    x = mcHmpWaldResults()$"seqSizes", 
#                    y = mcHmpWaldResults()$"pow")
            
            par(mar = c(4, 4, 1, 1))
            if(input$powerSimStart == 0)
            {
              plot(0, 0 , type = "n", xlim = input$sampleSizes, ylim = c(0, 1), 
                  xlab = "sample size", ylab = "power")
              text(mean(input$sampleSizes), .5, labels = "NOT YET INITIATED", 
                  cex = 2)
            } else {}
            
            plot(mcHmpWaldResults()$"seqSizes", mcHmpWaldResults()$"pow", 
                xlim = input$sampleSizes + c(-2L, 2L), 
                ylim = c(0, 1.1), pch = 19, lwd = 1,
                type = "p", #main = "Power vs. Sample Size",
                xlab = "sample size", ylab = "power")
            lines(lineData, col = 1, lty = 1)
            
            abline(h = c(0, input$alpha, 1), lty = 4, col = "gray70", lwd = 2)
            text(x =  min(input$sampleSizes), y = input$alpha, pos = 3, 
                labels = paste0("alpha = ", input$alpha), cex = 1.3)
          })
      
      ## saving plot on file in a temporary directory 
      savePdfPlot <- reactive({
            fileName <- paste0(tempdir(), .Platform[["file.sep"]], "copyPowPlot.pdf")
            pdf(file = fileName)
            powSimPlot()
            dev.off()
            fileName
          })
      
      powPlotPoint <- reactive({
            input$"powPlotClick"
#            dev.copy2pdf(file = "copyPowPlot.pdf")
            powSimPlot()
            
            ## function that interpolates data
            findSpecificPowFun <- splinefun(
                x = mcHmpWaldResults()$"seqSizes", 
                y = mcHmpWaldResults()$"pow")
            
            if (!is.null(powPlotCoords$x))
            {
              resCoords <- c(
                  "x" = powPlotCoords$x, 
                  "y" = min(1, findSpecificPowFun(powPlotCoords$x))
              )
              
              points(powPlotCoords, pch = 3, lwd = 3, cex = 2.5, col = "red2")
              lines(x = c(resCoords["x"], resCoords["x"]), 
                  y = c(0, resCoords["y"]), 
                  lty = 4, col = "red2", lwd = 2)
              lines(x = c(0, resCoords["x"]), 
                  y = c(resCoords["y"], resCoords["y"]), 
                  lty = 4, col = "red2", lwd = 2)
              
              text(x = powPlotCoords$x, y = powPlotCoords$y, 
                  labels = paste0(
                      "Size=", round(powPlotCoords$x, 2L),
                      "\n Power=", 
                      round(resCoords["y"], 3L)), 
                  pos = 3, offset = 1, cex = 1.5)
            } else {} 
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
      inputPars <- reactive({
            tmp <- c(
                "Enterotypes Stratification"     = input$strata, 
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
                "Max Sample Size"      = input$sampleSizes[2L], 
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
            aux <- format(Sys.time(), "%Y-%m-%d_%X")
            for (i in 1L:3)
              aux <- sub(pattern = ":", replacement = ".", x = aux, fixed = TRUE)
            paste0('ReportMB_', aux, '.pdf')
          },
          
          content = function(file) {
            src <- normalizePath('template.Rmd')
            
            # temporarily switch to the temp dir, in case you do not have write
            # permission to the current working directory
            owd <- setwd(tempdir())
            on.exit(setwd(owd))
            file.copy(src, 'template.Rmd')
            
#            require(rmarkdown)
            out <- render('template.Rmd', pdf_document(highlight = "haddock"))
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
#library(msWaldHMP)
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


