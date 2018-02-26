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
#setwd("C:/Users/Federico/Documents/SharedDocs/Dropbox/Microbiome/Architect/shinyMB")
#setwd("path/to/local/shinyApp/directory")
source("./www/auxCode.R")

input <- list(entero = "YES", MC = 1000,
    kindPiOne = "Geometric", nOTUs = 50, 
    mostLeastAb1 = "most abundant", mostLeastAb2 = "most abundant", 
    diffOTUs1 = 1, diffOTUs2 = 2,
    relAbund1 = 20, relAbund2 = 25, 
    n1 = 30, n2 = 30, theta = 0.01, totCounts = 1e4, 
    alpha = .1,
    sampleSizes = c(20, 40))


#selMostAb <- 1L:50
#piOne <- list("piVec" = list(piOneGen(K=max(selMostAb), kind = "geom")), 
#    "theta" = list(.01))

selMostAb <- 1L:50
ent1 <- readRDS("./data/stoolMostAbEnt1_3.rds")[, selMostAb]
ent2 <- readRDS("./data/stoolMostAbEnt2_3.rds")[, selMostAb]
ent3 <- readRDS("./data/stoolMostAbEnt3_3.rds")[, selMostAb]
piVec <- list(
    "ent1" = msWaldHMP:::piMoM4Wald(ent1), 
    "ent2" = msWaldHMP:::piMoM4Wald(ent2),
    "ent3" = msWaldHMP:::piMoM4Wald(ent3))
theta <- list(
    "ent1" = msWaldHMP:::weirMoM4Wald(ent1), 
    "ent2" = msWaldHMP:::weirMoM4Wald(ent2),
    "ent3" = msWaldHMP:::weirMoM4Wald(ent3))
piOne <- list("piVec" = piVec, "theta" = theta)


changedOTUs <- 1L:3
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
load("./data/librarySizes.RData")
tmpLibSize <- libSizesOrigRaref[["Stool"]]
inds <- tmpLibSize > 2000
tmpLibSize <- tmpLibSize[inds]
libSizes <- list(
    sample(tmpLibSize, size = input$n1, replace = TRUE), 
    sample(tmpLibSize, size = input$n2, replace = TRUE))


aux <- lapply(seq_along(piOne$"piVec"), FUN = function(i, obj1, obj2)
    {
      countsGen(
          sampleSizes = c(input$n1, input$n2),
          alphas = cbind(obj1[[i]], obj2[[i]]),
          theta = piOne$"theta"[[i]],
          K = length(obj1[[i]]),
          N = input$totCounts, libSizes = libSizes,
          seed = 12345)
    }, obj1 = piOne$"piVec", obj2 = piTwo$"piTwo")
generatedCounts <- list(
    "dmDataList" = lapply(aux, elNamed, name = "dmDataList"),
    "piDirList" = lapply(aux, elNamed, name = "piDirList"))

#mains <- "R.A.D."
mains <- c(
    "Enterotype: 'Bacteroides'", 
    "Enterotype: 'Prevotella'", 
    "Enterotype: 'Ruminococcus'")



### draw graph of RADs and estimated average counts
drawPiPlot2 <- function(piOneObj, piTwoObj, theta, mars = c(4, 4, 3, .1) + .1, ...)
{
  par(mar = mars)
  plot(piOneObj, type = "n", 
      ylim = c(0, max(piOneObj, piTwoObj)),
      xlab = "OTU Index", 
#                main = expression(paste("Ranked Abundances for ", pi[1])),
      cex.lab = 1.2, lty = 1, ...)
  
  points(piOneObj, lwd = 1.5, cex = 2, col = "black", 
      pch = 20)
  lines(piOneObj, lty = 1, col = "black")
  points(piTwoObj, lwd = 2, cex = 2, col = "red", 
      pch = 3)
  lines(piTwoObj, lty = 2, col = "red")
  
  abline(h = 0, lty = 4, col = "gray60", lwd = 2)
  
  usr <- par("usr")
  xCenter <- (usr[1L] + usr[2L])/2
  shift <- strwidth("x", cex = 2)
  text(x = xCenter - shift, y = usr[4L],
      labels = expression(theta), pos = 1, cex = 1.5)
  ## add *offset = 1.1* if hat(theta) 
  text(x = (usr[1L] + usr[2L])/2 + shift, y = usr[4L], 
      labels = paste("  =", theta), pos = 1, cex = 1.1)
  
  legend(x = "topright", pch = c(19, 3, 1, 4), 
      col = c("black", "red", "gray60", "#0000FF60"), pt.lwd = 3, pt.cex = 2,
      lwd = 2, lty = 1:2, cex = 1.5,
      legend = c(
          expression(pi[1]), 
          expression(pi[2]))) 
}





#pdf(file = "plotRads2.pdf", width = 9, height = 6)
#png(file = "plotRads.png", width = 750, height = 1200, pointsize = 24)
#jpeg(file = "plotRads.jpg", width = 750, height = 1000, pointsize = 24, quality = 100)
#jpeg(file = "Prevotella.jpg", width = 1000, height = 600, pointsize = 24, quality = 100)
#layout(t(1L:3))
layout(1L:3)
tmp <- lapply(seq_along(piOne$"piVec"), 
    FUN = function(i, counts, piOne, piTwo, mains) 
    {
#      mars <- if(i > 1)
#            c(4, 2, 3, 0) + .1 else
#            c(4, 4, 3, .1) + .1
      drawPiPlot2(
          piOneObj = piOne$"piVec"[[i]], 
          piTwoObj = piTwo$"piTwo"[[i]], 
          main = mains[[i]], ylab = "Abundance Proportions", 
          theta = round(piOne$"theta"[[i]], 3L), mars = mars)
    }, 
    counts = generatedCounts$"piDirList", 
    piOne = piOne, piTwo = piTwo, mains = mains)
#dev.off()




### power vs. sample size graph
set.seed(12345)
nSubsets <- 5L
sampleSizes <- round(
    seq(from = input$sampleSizes[1L], 
        to = input$sampleSizes[2L], 
        length = nSubsets))
nGroups <- length(generatedCounts$"dmDataList"[[1L]])
nOtus <- length(piOne$"piVec"[[1L]])
nStrata <- length(piOne$"piVec")

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
alphaDM <- lapply(seq_along(piOne$"piVec"), 
    FUN = function(i, one, two)
    {
      cbind(one[[i]], two[[i]])
    }, one = piOne$piVec, two = piTwo$piTwo)
names(alphaDM) <- names(piOne$piVec)

## create *thetaDM* for the function
thetaDM <- matrix(rep.int(unlist(piOne$"theta"), nGroups), 
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
### MonteCarlo simulations
tmp <- lapply(seq_len(input$MC), 
    FUN = function(x)
    {
      msWaldHMP:::msWald(
          nReads = nReads, alphaDM = alphaDM, thetaDM = thetaDM)
    })

res <- array(unlist(tmp), dim = c(nStrata, nSubsets, input$MC), 
    dimnames = list(
        rownames(el(tmp)), colnames(el(tmp)), 
        paste0("MC", seq_len(input$MC))
    ))


rejRes <- rowSums(res > qAlpha, na.rm = TRUE, dims = 2) / input$MC

## if multilple strata, sum up together
if (nStrata > 1L)
{
  rejResGlob <- rowSums(colSums(res, na.rm = TRUE) > qAlphaGlob) / input$MC
  rejRes <- rbind(rejRes, "Global" = rejResGlob)
} else {}

mcHmpWaldResults <- list("pow" = tail(rejRes, n = 1), 
    "seqSizes" = sampleSizes)

findSpecificPowFun <- approxfun(
    x = mcHmpWaldResults$"seqSizes", 
    y = mcHmpWaldResults$"pow")

xCoord <- 30# min(sampleSizes) + .75 * diff(range(sampleSizes))

resCoords <- list(
    "x" = xCoord, 
    "y" = findSpecificPowFun(xCoord)
)


#pdf(file = "plotPowSampSize.pdf", width = 9, height = 6)
layout(1)
par(mar = c(4, 4, 3, 1))
plot(mcHmpWaldResults$"seqSizes", mcHmpWaldResults$"pow", 
    xlim = input$sampleSizes + c(-2L, 2L), 
    ylim = c(0, 1.1), pch = 19, lwd = 1,
    type = "o", main = "Power vs. Sample Size",
    xlab = "sample size", ylab = "power")
abline(h = c(0, input$alpha, 1), lty = 4, col = "gray70", lwd = 2)
text(x =  min(input$sampleSizes), y = input$alpha, pos = 3, 
    labels = paste0("alpha = ", input$alpha), cex = 1)

lines(x = c(resCoords[["x"]], resCoords[["x"]]), 
    y = c(0, resCoords[["y"]]), 
    lty = 4, col = "red2", lwd = 2)
lines(x = c(0, resCoords[["x"]]), 
    y = c(resCoords[["y"]], resCoords[["y"]]), 
    lty = 4, col = "red2", lwd = 2)
points(resCoords, pch = 3, lwd = 3, cex = 2.5, col = "red2")

text(resCoords, #y = resCoords[["y"]], 
#                text(x = resCoords["x"], y = resCoords["y"], 
    labels = paste0(
        "Size=", round(resCoords[["x"]], 2L),
#                        "Size=", round(resCoords["x"]),
        "\n Power=", 
        round(resCoords[["y"]], 3L)), 
    pos = 3, offset = 1, cex = 1.3)
#dev.off()



