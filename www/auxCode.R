################################################################################
### Title: utilities.R                                              ###
###                                                                          ###
### Project: MicrobiomeTry                                                ###
###                                                                          ###
### Version: 0.1 - 19/mar/2015                                    ###
###                                                                          ###
### Description: short description                                           ###
###                                                                          ###
### Authors: Federico Mattiello <Federico.Mattiello@UGent.be>                ###
###                                                                          ###
### Maintainer: Federico Mattiello <Federico.Mattiello@UGent.be>             ###
###                                                                          ###
### Versions:                                                                ###
###     > 0.1 - 19/mar/2015 - 08:13:08:                                      ###
###         creation                                                         ###
###                                                                          ###
################################################################################
###

### generation of abundance probability distribution for sample one
piOneGen <- function (K = 100, pBreak = .15, rates = c(1, .5), kind = NULL,  diss = 0.05)
{
  if(missing(kind))
  {
    kind <- "geom"
  } else
  {
    kind <- match.arg(kind, c("exp", "geom", "bstick", "brGeom"))
  }
  
  if(K < 1)
  {
    stop("K (num. of taxa) should be at least 1.")
  }
  
  if(diss > 1 || diss < 0){
    stop("\"diss\" should be a number between 0 and 1.")
  }
  
  switch(kind,
         "exp" = {
           maxQuant <- qexp(.99, rate = min(rates))
           hiProbTaxa <- dexp(
             seq.int(from = 0, to = maxQuant * pBreak, 
                     length = round(K * pBreak)),
             rate = max(rates))
           
           loProbTaxa <- dexp(
             seq.int(from = maxQuant * pBreak, to = maxQuant * .5, 
                     length = round(K * (1 - pBreak))),
             rate = min(rates))
           rr <- c(hiProbTaxa, loProbTaxa)
           sq <- round(K * pBreak):min(K, round(K * pBreak) + 1L)
           rr[sq] <- mean(rr[sq])
         }, 
         "geom" = {
           const <- 1/(1 - ((1 - diss)^K))
           rr <- 1000 * const * diss * (1 - diss)^seq_len(K)
         }, 
         "bstick" = {
           rr <- bstick(K, 1)
         }, 
         "brGeom" = {    # broken geometric progression with change-point 
           # at one fifth of the number of taxa
           changePoint <- K %/% 5
           diss <- diss * 1.25
           const1 <- 1/(1 - ((1 - diss)^changePoint))
           rr1 <- 1000 * const1 * diss * (1 - diss)^seq_len(changePoint)
           
           diss2 <- max(diss / 2, .005)
           const2 <- 1/(1 - ((1 - diss2)^(K - changePoint)))
           rr2 <- 1000 * const2 * diss2 * 
             (1 - diss2)^(seq_len(K - changePoint) + changePoint)
           mult <- .95 * min(rr1) / max(rr2)
           rr2 <- rr2 * mult
           rr <- c(rr1, rr2)
         }
  )
  rr <- rr[order(rr, decreasing = TRUE)]
  rr / sum(rr)
}# END - function: piOneGen

#### tests
#ex <- piOneGen(kind = "exp")
#gs <- piOneGen(kind = "geom")
#bs <- piOneGen(kind = "bstick")
#bgs <- piOneGen(kind = "brGeom")
#matplot(cbind(ex, gs, bs, bgs), cex= .8, type = "o", 
#        xlab = "OTU index", ylab = "Relative Abundance")
#legend(x = "topright", pch = as.character(1:4), lty = 1:4, col = 1:4, pt.cex = 1, 
#        pt.lwd = 2, lwd = 2, 
#        legend = c("exponential decay \n(with change-point)", 
#                "geometric progr.", "broken stick", 
#                "broken geometric progr.\n(with change-point)"))


### rho_j generation, pi2_j = rho_j * pi1_j
#' Generate \eqn{\\rho}{"rho"} values: abundance proportions of first 
#' sample w.r.t. the second one
#' 
#' @param K 
#'     total number of OTUs, default is \code{100}
#' @param m1 
#'     number of OTUs with \emph{high} differential abundance between the two samples
#' @param m2 
#'     number of OTUs with \emph{low} differential abundance between the two samples
#' @param relAbund1 
#'     higher value for differential abundance, default is \code{1.33} equal to 
#'     one third difference
#' @param relAbund2 
#'     lower value for differential abundance, default is \code{1.25} equal to
#'     one fourth difference
#' @param otus2Change 
#'     indices of OTUs that need to be changed, it has to be equal to \code{m1 + m2},
#'     as default the most abundant OTUs are considered 
#' @return 
#' @author Federico Mattiello <Federico.Mattiello@@UGent.be>
#' @export
#' 
rhoGenSpecOtus <- function(K, m1 = 2, m2 = 3, relAbund1 = 1.33, relAbund2 = 1.25, 
    otus2Change = NA)
{
  
  rho <- rep.int(1, K)
  
  if (is.na(otus2Change[1L]))
  {
    #        otus2Change <- (K - m1 - m2 + 1):K
    otus2Change <- seq_len(m1 + m2)
  } else {}
  
  rho[otus2Change[seq_len(m1)]] <- relAbund1
  rho[otus2Change[seq_len(m2) + m1]] <- relAbund2
  
  list("rho" = rho, "changedOTUs" = otus2Change, 
       "relAbund1" = relAbund1, "relAbund2" = relAbund2, 
       "m1" = m1, "m2" = m2)
}# END - function: rhoGen





piTwoGenSimple <- function(piOne, rho, leastAbund = FALSE)
{
  ### classify OTUs in: (i) non DA, (ii) DA, and (iii) DA for compensation
  otuType <- rep.int(0L, length(rho$rho))
  indEqual <- abs(rho$rho - 1) <= .Machine$double.eps
  ind1 <- abs(rho$rho - rho$relAbund1) <= .Machine$double.eps & !indEqual
  ind2 <- abs(rho$rho - rho$relAbund2) <= .Machine$double.eps & !indEqual
  otuType[ind1] <- 2L
  otuType[ind2] <- 1L
  
  ### actually modify piOne
  piTwo <- piOne * c(rho$rho, rep.int(1, length(piOne) - length(rho$rho)))
  toCompensate <- sum(piTwo - piOne)
  
  ### default is _leastAbund=FALSE_ and compensation is done by dissipation among 
  ### all unrequested OTUs, otherwise least abundant OTUs are set to 0 until
  ### amount has been fully compensated.
  if (leastAbund)
  {
    startInd <- max(which(otuType == 0))
    indRun <- startInd
    numIter <- 1L
    
    while(abs(sum(piTwo) - 1) > 1e-12 && (numIter <= sum((otuType == 0))))
    {
      aux <- piTwo[indRun]
      piTwo[indRun] <- max(aux - toCompensate, 0)
      toCompensate <- sign(toCompensate) * abs(aux - toCompensate)
      indRun <- indRun - 1L
      numIter <- numIter + 1L
    }# END - while: compensate sequentially
    
    if (startInd != indRun)
    {
      otuType[startInd:indRun] <- -1
    } else {}
  } else
  {
    ### compute how much probability mass we have to compensate and
    ### check if compensation is too big, in that case rescale requested OTUs
    ## n0 = num. of 0-type OTUs, i.e. non-requested OTUs, procedure:
    ## Let _d_ = _toCompensate_, _m0_ = sum of rel.abund. of unrequested OTUs,
    ## _m1_ = sum of rel.abund. of requested OTUs,
    ## if _d > m0_ then rescale requested OTUs multiplying each element by
    ## [m1 - (d-m0)]/m1, as a consequence all unrequested OTUs are set to 0.
    ## Otherwise, each unrequested OTU of _piTwo_ is multiplied by (m0 - d)/m0. 
    m0 <- sum(piOne[otuType == 0L])
    m1 <- sum(piTwo[otuType != 0L])
    compensDiff <- toCompensate - m0 
    if (compensDiff > 1e-12)
    {
      piTwo[otuType != 0L] <- piTwo[otuType != 0L] * (m1 - compensDiff) / m1
    } else {}
    
    piTwo[otuType == 0L] <- piTwo[otuType == 0L] * (m0 - toCompensate)/m0
  }# END - previous compensation strategy or _dissipate_
  
  list("piTwo" = piTwo / sum(piTwo), "otuType" = as.factor(otuType))
}# END - function: new piTwoGen, simple version



piTwoGen <- function(piOne, rho, 
    compensation = c("dissipate", "simple", "mostAbund", "relDiff", "zeros", "piAux"), 
    otus4Balance = NA)
{
  if(missing(compensation))
  {
#    compensation <- "relDiff"
#     compensation <- "mostAbund"
    compensation <- "dissipate"
  } else {}
  
  ### classify OTUs in: (i) non DA, (ii) DA, and (iii) DA for compensation
  otuType <- rep.int(0L, length(rho$rho))
  indEqual <- abs(rho$rho - 1) <= .Machine$double.eps
  ind1 <- abs(rho$rho - rho$relAbund1) <= .Machine$double.eps & !indEqual
  ind2 <- abs(rho$rho - rho$relAbund2) <= .Machine$double.eps & !indEqual
  otuType[ind1] <- 2L
  otuType[ind2] <- 1L
  
  ### actually modify piOne
  piTwo <- piOne * rho$rho
#  changedOTUs <- rho$changedOTUs
  changedOTUs <- which(otuType != 0)
  
  ### compute how much probability mass we have to compensate and
  ### check if compensation is too big
  toCompensate <- sum(piTwo - piOne)
  compensDiff <- toCompensate - sum(piOne[otuType == 0])
  increase <- sum(abs(piTwo - piOne))
  if (compensDiff > 1e-12)
  {
    equalOTUs <- otuType != 0L
    piTwo[equalOTUs] <- piTwo[equalOTUs] - compensDiff / sum(equalOTUs)
    piTwo[piTwo < 0] <- 0
    rho$rho <- piTwo / piOne
    rho$relAbund1 <- el(rho$rho[otuType == 2L])
    rho$relAbund2 <- el(rho$rho[otuType == 1L])
  } else {}
  
  if (!is.na(otus4Balance[1L]))
  {
    compensation <- "fixedOTUs"
  } else {}
  
  ## method of the fixed relative difference
  switch (compensation, 
      "dissipate" = {
        ## n0 = num. of 0-type OTUs, i.e. non-requested OTUs, procedure:
        ## Let _d_ = _toCompensate_, _m0_ = sum of rel.abund. of unrequested OTUs,
        ## and _r_ = (m0 - d)/m0. Each unrequested OTU of _piTwo_ is multiplied by _r_.
        ## NOTE: this is the same strategy applied in the previous _if_ statement where
        ## the amount to compensate is compared to the total mass available for 
        ## compensation.
        m0 <- sum(piTwo[otuType == 0L])
        piTwo[otuType == 0L] <- piTwo[otuType == 0L] * (m0 - toCompensate)/m0
        piTwo[piTwo < 0] <- 0
      },
      "simple" = {
        ## n0 = num. of 0-type OTUs, i.e. non-requested OTUs, procedure:
        ## (i) subtract _toCompensate/n0_ from all 0-type OTUs, OTUs with abundances
        ##   lower than _toCompensate/n0_ are set to 0
        ## (ii) if still there is a difference, all necessary least abundant OTUs are 
        ##   set to 0 as well. In case _toCompensate < 0_, only (i) is needed
        auxPiTwo <- piTwo[otuType == 0L]
        auxAmount <- toCompensate / sum(otuType == 0L)
        auxPiTwo <- ifelse(auxPiTwo < auxAmount, yes = 0, no = auxPiTwo - auxAmount)
        piTwo[otuType == 0L] <- auxPiTwo
        
        newAmount <- sum(piTwo - piOne)
        if (newAmount > 1e-10)
        {
          revAuxPiTwoCS <- cumsum(rev(auxPiTwo))
          indToCompensate <- max(which(revAuxPiTwoCS <= newAmount))
          auxPiTwo[(length(auxPiTwo) - indToCompensate):length(auxPiTwo)] <- 0
          piTwo[otuType == 0L] <- auxPiTwo
        } else {}
        
#         otuType[piTwo < 10 * .Machine$double.eps] <- -1L
      },        # END - simple compensation: spread changes on other OTUs
      "mostAbund" = {
        
        startInd <- min(which(otuType == 0))
        indRun <- startInd
        numIter <- 1L
        
        while(abs(sum(piTwo) - 1) > 1e-12 && numIter <= (length(piTwo) - startInd))
        {
          aux <- piTwo[indRun]
          piTwo[indRun] <- max(aux - toCompensate, 0)
          toCompensate <- sign(toCompensate) * abs(aux - toCompensate)
          indRun <- indRun + 1L
          numIter <- numIter + 1L
        }# END - while: compensate sequentially
        
        if (startInd != indRun)
        {
          otuType[startInd:indRun] <- -1
        } else {}
      },
      
      "zeros" = {
        ### set to zero all least abundant OTUs or the ones next to the changed ones
        otus4Balance <- setdiff(seq_along(piTwo), changedOTUs)
        piTwoCumsum <- cumsum(rev(piTwo[-changedOTUs]))
        indBal <- max(which(piTwoCumsum <= increase))
        otus4Balance <- tail(otus4Balance, n = indBal)
        piTwo[otus4Balance] <- 0
        otuType[otus4Balance] <- -1
      },
      
      "relDiff" = {
        piOneAux <- piOne[-changedOTUs]
        
        ## if OTUs that are need for compensation are not specified, 
        ## take the ones just before the changed ones, but with relative 
        ## abundance difference that is less of *relAbundance2* when present
        otus4Balance <- setdiff(seq_along(piOne), changedOTUs)
        halfRelAbund2 <- max(1.25, .75 * (rho$relAbund2 - 1) + 1)
        piAuxCumsum <- cumsum(rev(piOneAux - piOneAux / halfRelAbund2))
        
        if(any(piAuxCumsum > 2 * increase))
        {
          indBal2 <- max(which(piAuxCumsum <= increase))
          otus4Balance <- tail(otus4Balance, n = indBal2)
          piTwo[otus4Balance] <- piOne[otus4Balance] / halfRelAbund2
        } else
        {
          ## shrink to zero the least abundant OTUs
          piOneCumsum <- cumsum(rev(piOneAux))
          indBal1 <- max(which(piOneCumsum <= increase/2) + 1L)
          indBal2 <- max(which(piAuxCumsum[-seq_len(indBal1)] <= increase/2))
          otus4Balance <- tail(otus4Balance, n = indBal1 + indBal2 + 1)
          piTwo[otus4Balance[seq_len(indBal2)]] <- 
              piOne[otus4Balance[seq_len(indBal2)]] / halfRelAbund2
          piTwo[otus4Balance[seq_len(indBal1) + indBal2 + 1]] <- 0
        }# END - ifelse: not enough OTUs for balancing
        
        otuType[otus4Balance] <- -1
      },        # END of *relDiff* compensation method
      
      "piAux" = {
        piAux <- piOneGen(length(piOne), pBreak = .05, rates = c(1.2, .2))[-changedOTUs]
        piTwoAux <- piTwo[-changedOTUs]
        cumSumAux <- cumsum(rev(abs(piTwoAux - piAux)))
        
        nOtus4Balance <- max(which(cumSumAux <= increase))
        
        revPiAux <- rev(piAux)
        revPiTwoAux <- rev(piTwoAux)
        
        revPiTwoAux[seq_len(nOtus4Balance)] <- revPiAux[seq_len(nOtus4Balance)]
        
        if (piTwoAux[1L] < piTwo[changedOTUs[1L]])
        {
          piTwo <- c(piTwo[changedOTUs], rev(revPiTwoAux))
        } else
        {
          piTwo <- c(rev(revPiTwoAux), piTwo[changedOTUs])
        }
        
        indOtus4Balance <- tail(which(otuType < 1e-9), n = nOtus4Balance)
        otuType[indOtus4Balance] <- -1
      },
      
      "fixedOTUs" = {
        piOneAux <- piOne[-changedOTUs]
        aux <- piTwo[otus4Balance] - increase / length(otus4Balance)
        if (any(aux < 0))
        {
          aux[aux < 0] <- 0
          warning ("Not enough OTUs for re-normalisation.")
        } else {}
        
        piTwo[otus4Balance] <- aux
        otuType[otus4Balance] <- -1
      }
  )# END *compensation* methods
  
  otuType <- as.factor(otuType)
  
  piTwo <- piTwo / sum (piTwo)
  list("piTwo" = piTwo, "otuType" = otuType)
}



### generate Dirichlet realisations, taken from gtools (identical in MCMCpack)
rDirichlet <- function (n, alpha) 
{
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x/as.vector(sm)
}

### generate stick-breaking distribution
bstick <- function (n, tot.var = 1, ...)
{
  res <- rev(cumsum(tot.var/n:1)/n)
  names(res) <- paste("Stick", seq(len = n), sep = "")
  res
}


### generate simulated counts matrix, so far works only with two samples
##totCounts <- 1e4 + rnorm(1e3, sd = 1e3)
countsGen <- function(sampleSizes, alphas, theta, K, N, seed = 1234, libSizes = NULL)
{
  
  if (!is.na(seed))
  {
    set.seed(seed)
  } else {}
  
  #    require(gtools)
  dmDataList <- piDirList <- vector("list", length(sampleSizes))
  
  for(nRun in seq_along(sampleSizes))
  {
    alphaTmp <- alphas[, nRun] * (1-theta)/theta
    alphaTmp[alphaTmp < 0] <- 0
    tmpPiDir <- rDirichlet(n = sampleSizes[nRun], alpha = alphaTmp)
    inds <- (tmpPiDir <= .Machine$double.eps) | is.na(tmpPiDir)
    tmpPiDir[inds] <- 0
    
    tmp <- matrix(NA, nrow = sampleSizes[nRun], ncol = K)
    
    for(iRun in seq_len(sampleSizes[nRun]))
      tmp[iRun, ] <- rmultinom(n = 1L, size = libSizes[[nRun]][iRun], 
          prob = tmpPiDir[iRun, ])
    
    dmDataList[[nRun]] <- tmp
    piDirList[[nRun]] <- tmpPiDir
  }
  
  list("dmDataList" = dmDataList, "piDirList" = piDirList)
}# END - function: countsGen



### draw graph of RADs and estimated average counts
drawPiPlot <- function(countsData, piOneObj, piTwoObj, theta, ...)
{
  #  avgEstPiData <- sapply(countsData, FUN = colMeans)
  #  sdEstPiData <- sapply(countsData, FUN = function(x) apply(x, 2, sd))
  #  
  #  yAuxCoords1 <- matrix(c(
  #              avgEstPiData[, 1] - sdEstPiData[, 1],
  #              avgEstPiData[, 1], 
  #              avgEstPiData[, 1] + sdEstPiData[, 1]), ncol = 3L)
  #  yAuxCoords2 <- matrix(c(
  #              avgEstPiData[, 2] - sdEstPiData[, 2],
  #              avgEstPiData[, 2], 
  #              avgEstPiData[, 2] + sdEstPiData[, 2]), ncol = 3L)
  #  pchs <- c(45, 4, 45)
  
  #  sdCex <- 2.5 * (sdEstPiData / max(sdEstPiData)) + .5
  
  par(mar = c(4, 4, 3, .1) + .1)
  plot(piOneObj, type = "n", 
       ylim = range(0, piOneObj, piTwoObj), #yAuxCoords1, yAuxCoords2), 
       xlab = "OTU Index of Controls", cex.lab = 1.2, 
       #                main = expression(paste("Ranked Abundances for ", pi[1])),
       ...)
  
  #  points(piOneObj, lwd = 1.5, col = "black", pch = 3, cex = sdCex[, 1])
  points(piOneObj, lwd = 2, col = "black", pch = 3, cex = 2)
  lines(piOneObj, lty = 1, col = "black")
  #  points(piTwoObj, lwd = 2, col = "red", pch = 4, cex = sdCex[, 2])
  points(piTwoObj, lwd = 2, col = "red", pch = 4, cex = 2)
  lines(piTwoObj, lty = 2, col = "red")
  
  #  for (i in seq_along(piOneObj))
  #  {
  #    lines(c(i, i, i), yAuxCoords1[i, ], lty = 1, col = "black")
  #    points(c(i, i, i), yAuxCoords1[i, ], col = "black", pch = pchs, cex = c(2, 1, 2), lwd = 3)
  #    lines(c(i, i, i), yAuxCoords2[i, ], lty = 1, col = "red")# "#0000FF60"
  #    points(c(i, i, i), yAuxCoords2[i, ], col = "red", pch = pchs, cex = c(2, 1, 2), lwd = 3)
  #  }
  
  
  abline(h = 0, lty = 4, col = "gray60", lwd = 2)
  
  usr <- par("usr")
  xCenter <- (usr[1L] + usr[2L])/2
#   shift <- strwidth(expression(theta), cex = 2)
  shift <- strwidth("M", cex = 1.5)
  text(x = xCenter - shift, y = usr[4L],
       labels = expression(theta), pos = 1, cex = 1.5)
  ## add *offset = 1.1* if hat(theta) 
  text(x = (usr[1L] + usr[2L])/2 + shift, y = usr[4L], 
       labels = paste("  =", theta), pos = 1, cex = 1.1)
  
  legend(x = "topright", pch = c(3, 4), 
         col = c("black", "red", "gray60", "#0000FF60"), pt.lwd = 3, pt.cex = 2,
         lwd = 2, lty = 1:2, cex = 1.5,
         legend = c(
           expression(pi[1]), 
           expression(pi[2])))#, 
  #          expression(bar(pi[1])), 
  #          expression(bar(pi[2]))))
}




#
#
#
#data(saliva)
#data(throat)
#data(tonsils)
#
#### Get a list of dirichlet-multinomial parameters for the data
#fit.saliva <- DM.MoM(saliva) 
#fit.throat <- DM.MoM(throat)
#fit.tonsils <- DM.MoM(tonsils)
#
#### Set up the number of Monte-Carlo experiments
#### We set MC=1 due to CRAN restrictions, Please set MC to be at least 1,000
#MC <- 1000
#
#### Generate a random vector of number of reads per sample
#Nrs1 <- rep(12000, 9)
#Nrs2 <- rep(12000, 11)
#Nrs3 <- rep(12000, 12)
#group.Nrs <- list(Nrs1, Nrs2, Nrs3)
#
#### Computing size of the test statistics (Type I error)
#mc.xdc_check1 <- MC.Xdc.statistics(group.Nrs, MC, fit.saliva$gamma, 3, "hnull", .05, "mom")
#mc.xdc_check1
#
#### Computing Power of the test statistics (1 - Type II error)
#group.alphap <- rbind(fit.saliva$gamma, fit.throat$gamma, fit.tonsils$gamma)
#mc.xdc_check2 <- MC.Xdc.statistics(group.Nrs, MC, group.alphap, 3, "ha", 0.05, "mom")
#mc.xdc_check2

#library(DESeq2)
#
#data1 <- read.csv("shinyMB-0.2/data/stoolMostAbundEnt1_3.csv")[, -1]
#data3 <- read.csv("shinyMB-0.2/data/stoolMostAbundEnt3_3.csv")[, -1]
#nams <- intersect(colnames(data1), colnames(data3))
#minN <- min(nrow(data1), nrow(data3))
#
#countData <- t(rbind(data3[seq_len(minN), nams], data1[seq_len(minN), nams]))
#colData <- data.frame(enterotype = 
#        factor(rep.int(c("ent1","ent3"), c(minN, minN))))
#dds <- DESeqDataSetFromMatrix(countData = countData, colData, formula(~ enterotype))
#dsObj <- DESeq(dds, test = "Wald")
#resDs <- results(dsObj)
#
#ord <- order(resDs@listData$padj)[seq_len(sum(resDs@listData$padj < .1))]
#cbind(resDs@rownames[ord], resDs@listData$padj[ord])




#rm(list = ls()); gc(reset = TRUE)



