################################################################################
### Title: script4Graphs.R                                              ###
###                                                                          ###
### Project: MicrobiomeTry                                                ###
###                                                                          ###
### Version: 0.1 - 30/mar/2015                                    ###
###                                                                          ###
### Description: short description                                           ###
###                                                                          ###
### Authors: Federico Mattiello <Federico.Mattiello@UGent.be>                ###
###                                                                          ###
### Maintainer: Federico Mattiello <Federico.Mattiello@UGent.be>             ###
###                                                                          ###
### Versions:                                                                ###
###     > 0.1 - 30/mar/2015 - 16:13:06:                                      ###
###         creation                                                         ###
###                                                                          ###
################################################################################
###

#rm(list = ls()); gc(reset = TRUE)

setwd("C:/Users/Federico/Documents/SharedDocs/Dropbox/Microbiome/Architect/shinyMB-0.2")
source("./www/auxCode.R")

input <- list(entero = "YES", MC = 100,
    kindPiOne = "Geometric", nOTUs = 50, 
    mostLeastAb1 = "most abundant", mostLeastAb2 = "most abundant", 
    diffOTUs1 = 1, diffOTUs2 = 2,
    relAbund1 = 20, relAbund2 = 25, 
    n1 = 30, n2 = 30, theta = 0.01, totCounts = 1e4, 
    alpha = .1,
    sampleSizes = c(10, 40))


selMostAb <- 1L:25
ent1 <- readRDS("./data/stoolMostAbEnt1_3.rds")[, selMostAb]
ent2 <- readRDS("./data/stoolMostAbEnt2_3.rds")[, selMostAb]
ent3 <- readRDS("./data/stoolMostAbEnt3_3.rds")[, selMostAb]
piVec <- list(
    "ent1" = multSampWald:::piMoM4Wald(ent1), 
    "ent2" = multSampWald:::piMoM4Wald(ent2),
    "ent3" = multSampWald:::piMoM4Wald(ent3))
theta <- list(
    "ent1" = multSampWald:::weirMoM4Wald(ent1), 
    "ent2" = multSampWald:::weirMoM4Wald(ent2),
    "ent3" = multSampWald:::weirMoM4Wald(ent3))

piOne <- list("piVec" = piVec, "theta" = theta)

changedOTUs <- 1L:5
nOtus <- length(piOne$"piVec"[[1L]])
rho <- rhoGenSpecOtus(
    K = nOtus, 
    m1 = input$diffOTUs1, m2 = input$diffOTUs2, 
    relAbund1 = (input$relAbund1 + 100) / 100, 
    relAbund2 = (input$relAbund2 + 100) / 100, 
    otus2Change = changedOTUs)

aux <- lapply(piOne$"piVec", FUN = piTwoGen, rho = rho)
piTwo <- list(
    "piTwo" = lapply(aux, elNamed, name = "piTwo"),
    "otuType" = lapply(aux, elNamed, name = "otuType"))


aux <- lapply(seq_along(piOne$"piVec"), FUN = function(i, obj1, obj2)
    {
      countsGen(
          sampleSizes = c(input$n1, input$n2),
          alphas = cbind(obj1[[i]], obj2[[i]]),
          theta = piOne$"theta"[[i]],
          K = length(obj1[[i]]),
          N = input$totCounts,
          seed = 12345)
    }, obj1 = piOne$"piVec", obj2 = piTwo$"piTwo")
generatedCounts <- list(
    "dmDataList" = lapply(aux, elNamed, name = "dmDataList"),
    "piDirList" = lapply(aux, elNamed, name = "piDirList"))


mains <- c(
    "Enterotype: 'Bacteroides'", 
    "Enterotype: 'Prevotella'", 
    "Enterotype: 'Ruminococcus'")

#png(file = "plotRads.png", width = 750, height = 1200, pointsize = 24)
#jpeg(file = "plotRads.jpg", width = 750, height = 1000, pointsize = 24, quality = 100)
layout(1L:3)
tmp <- lapply(seq_along(piOne$"piVec"), 
    FUN = function(i, counts, piOne, piTwo, mains) 
    {
      drawPiPlot(
              countsData = counts[[i]], piOneObj = piOne$"piVec"[[i]], 
              piTwoObj = piTwo$"piTwo"[[i]], 
              main = mains[[i]], ylab = "Abundance Proportions", 
              theta = round(piOne$"theta"[[i]], 3L))
    }, 
    counts = generatedCounts$"piDirList", 
    piOne = piOne, piTwo = piTwo, mains = mains)
#dev.off()



