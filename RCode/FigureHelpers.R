library(R.matlab)
library(tidyverse)
library(ggplot2)
library(ggpubr)


plotFluxDiv = function(RX, title) {
  mets = as.character(unlist(RX[4,1,1]))
  #print(mets)
  fluxDivs = RX[2,1,1]$fluxDivUbs
  a = as.numeric(unlist(RX[5,1,1]))
  
  #metsToPlot = c("glucose", "glutamine", "glycine", "D-lactate", "oxygen")
  glucose = fluxDivs[,mets == "glucose"]
  glutamine = fluxDivs[,mets == "glutamine"]
  NEFA = fluxDivs[,mets == "NEFA blood pool in"]
  chol = fluxDivs[,mets == "cholesterol"]
  
  #glycine = fluxDivs[,mets == "glycine"]
  lactate = fluxDivs[,mets == "L-lactate"]
  oxygen = fluxDivs[,mets == "O2"]
  #print(glucose)
  
  numMet = 6
  group = factor(rep(1:numMet, 1, each=length(a)), 1:numMet, c("glucose", "glutamine", "lipid pool", "cholesterol", "oxygen", "lactate"))
  ds = tibble(x=rep(a,numMet), y=c(glucose, glutamine, NEFA, chol, oxygen, lactate), Metabolite = group)
  col = c(rgb(0,0,0),gg_color_hue(6),rgb(0.5,0.5,0.5))
  labels = c("glucose", "glutamine", "lipid pool", "cholesterol       ", "oxygen", "lactate")
  
  pA = ggplot(ds,ggplot2::aes(x=x,y=y, group=Metabolite, size=Metabolite, colour=Metabolite)) +
    geom_line(ggplot2::aes(x=x,y=y, group=Metabolite, colour=Metabolite)) +
    scale_linetype_manual(values = rep(1,6), labels = labels) +
    scale_size_manual(values = rep(1,6), labels = labels) +
    scale_color_manual(values = col[c(2,3,4,8,6,7)], labels = labels) +
    labs(y="Min required fraction", x="a", title=title) +
    theme_bw() +
    theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank())
  
  
  return(pA)  
}

plotRed = function(RX, title) {
  mets = as.character(unlist(RX[4,1,1]))
  resultSolutions = RX[1,1,1]$resultSolutionsBasic
  
  a = as.numeric(unlist(RX[5,1,1]))
  
  nPoints = length(a)
  
  vals = matrix(NA,nrow=nPoints, ncol=length(mets))
  lowestAInd = NA

  redGrowth = RX[3,1,1]$metsRedGrowth
  
  for (i in 1:nPoints) {
    if (!is.null(resultSolutions[[i]])) {
      
      normGrowth = -as.numeric(resultSolutions[i][[1]][[1]][2,1,1])
      status = as.numeric(resultSolutions[i][[1]][[1]][3,1,1])
      if (status == 1) {
        if (is.na(lowestAInd)) {
          lowestAInd = i;
        }
        for (f in 1:length(mets)) {
          rg = redGrowth[i,f]
          vals[i, f] = rg/normGrowth
        }
      }
    }
  }
  
  
  
  glucose = vals[,mets == "glucose"]
  glutamine = vals[,mets == "glutamine"]
  NEFA = vals[,mets == "NEFA blood pool in"]
  chol = vals[,mets == "cholesterol"]
  
  lactate = vals[,mets == "L-lactate"]
  oxygen = vals[,mets == "O2"]

  numMet = 6
  group = factor(rep(1:numMet, 1, each=length(a)), 1:numMet, c("glucose", "glutamine", "lipid pool", "cholesterol", "oxygen", "lactate"))
  ds = tibble(x=rep(a,numMet), y=c(glucose, glutamine, NEFA, chol, oxygen, lactate), Metabolite = group)
  col = c(rgb(0,0,0),gg_color_hue(6),rgb(0.5,0.5,0.5))
  labels = c("glucose", "glutamine", "lipid pool", "cholesterol       ", "oxygen", "lactate")
  
  pA = ggplot2::ggplot(ds,ggplot2::aes(x=x,y=y, group=Metabolite, size=Metabolite, colour=Metabolite)) +
    ggplot2::geom_line(ggplot2::aes(x=x,y=y, group=Metabolite, colour=Metabolite)) +
    scale_linetype_manual(values = rep(1,6), labels = labels) +
    scale_size_manual(values = rep(1,6), labels = labels) +
    scale_color_manual(values = col[c(2,3,4,8,6,7)], labels = labels) +
    ggplot2::labs(y="Growth ratio reduced vs normal", x="a", title=title) +
    ggplot2::theme_bw() +
    theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank())
  
  return(pA)  
}

plotFluxDiv2 = function(RX, title) {
  mets = as.character(unlist(RX[4,1,1]))
  fluxDivs = RX[2,1,1]$fluxDivUbs
  a = as.numeric(unlist(RX[5,1,1]))
  
  glucose = fluxDivs[,mets == "glucose"]
  glutamine = fluxDivs[,mets == "glutamine"]
  NEFA = fluxDivs[,mets == "NEFA blood pool in"]

  lactate = fluxDivs[,mets == "L-lactate"]
  oxygen = fluxDivs[,mets == "O2"]

  numMet = 5
  group = factor(rep(1:numMet, 1, each=length(a)), 1:numMet, c("glucose", "glutamine", "lipid pool", "oxygen", "lactate"))
  ds = tibble(x=rep(a,numMet), y=c(glucose, glutamine, NEFA, oxygen, lactate), Metabolite = group)
  col = c(rgb(0,0,0),gg_color_hue(6),rgb(0.5,0.5,0.5))
  labels = c("glucose", "glutamine", "lipid pool", "oxygen", "lactate")
  
  pA = ggplot(ds,ggplot2::aes(x=x,y=y, group=Metabolite, size=Metabolite, colour=Metabolite)) +
    geom_line(ggplot2::aes(x=x,y=y, group=Metabolite, colour=Metabolite)) +
    scale_linetype_manual(values = rep(1,numMet), labels = labels) +
    scale_size_manual(values = rep(1,numMet), labels = labels) +
    scale_color_manual(values = col[c(2,3,4,6,7)], labels = labels) +
    labs(y="Min required fraction", x="a", title=title) +
    theme_bw()  + 
    theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank())
  
  
  return(pA)  
}



plotFluxes = function(RX, title) {
  a = as.numeric(unlist(RX[1,1,1]))
  rxns = as.character(unlist(RX[3,1,1]$rxns))
  resultSolutions = RX[2,1,1]$resultSolutions
  resultSolutions[[16]]
  nPoints = length(a)
  
  biomass = rep(NA, nPoints)
  glucoseUptake = rep(NA, nPoints)
  glutamineUptake = rep(NA, nPoints)
  glycineUptake = rep(NA, nPoints)
  oxygenUptake = rep(NA, nPoints)
  lactateExport = rep(NA, nPoints)
  complexI = rep(NA, nPoints)
  complexIII = rep(NA, nPoints)
  complexIV = rep(NA, nPoints)
  complexV = rep(NA, nPoints)
  protPool = rep(NA, nPoints)
  
  selBiomass = rxns == 'biomass_human'
  selGlucoseUptake = rxns == 'HMR_9034_REV'
  selGlutamineUptake = rxns == 'HMR_9063_REV'
  selGlycineUptake = rxns == 'HMR_9067_REV'
  selOxygenUptake = rxns == 'HMR_9048_REV'
  selLactateExport = rxns == 'HMR_9135'
  selComplexI = rxns == 'HMR_6921No1'
  selComplexIII = rxns == 'HMR_6918No1'
  selComplexIV = rxns == 'HMR_6914No1'
  selComplexV = rxns == 'HMR_6916No1'
  selProtPool = rxns == 'prot_pool_exchange'
  
  
  for (i in 1:nPoints) {
    if (!is.null(resultSolutions[[i]])) {
      fluxes = as.numeric(unlist(resultSolutions[[i]][[1]][1,1,1]))
      status = as.numeric(unlist(resultSolutions[[i]][[1]][3,1,1]))
      if (status == 1) {
        biomass[i] = fluxes[selBiomass];
        glucoseUptake[i] = fluxes[selGlucoseUptake]/100;
        glutamineUptake[i] = fluxes[selGlutamineUptake]/100;
        glycineUptake[i] = fluxes[selGlycineUptake]/100;
        oxygenUptake[i] = fluxes[selOxygenUptake]/100;
        lactateExport[i] = fluxes[selLactateExport]/100;
        complexI[i] = fluxes[selComplexI]/100;
        complexIII[i] = fluxes[selComplexIII]/100;
        complexIV[i] = fluxes[selComplexIV]/100;
        complexV[i] = fluxes[selComplexV]/100;
        protPool[i] = fluxes[selProtPool]/10;
        
      }
    }
  }
  
  #immuno-deficient rat tumor growth
  irat = rep(0.0080976, length(a))
  
  
  numFluxes = 12
  labels = c("biomass", "glucose upt.", "glutamine upt.", "glycine upt.", "oxygen upt.", "lactate exp.", "complex I", "complex III", "complex IV", "complex V", "Prot. pool","exp. growth")
  group = factor(rep(1:numFluxes, 1, each=length(a)), 1:numFluxes, labels)
  lst = c(1,1,1,1,1,1,2,2,2,2,2,3)
  cs = c(1,2,3,4,5,6,1,2,3,4,5,1)
  
  ds = tibble(x=rep(a,numFluxes), 
              y=c(biomass, glucoseUptake, glutamineUptake, glycineUptake, oxygenUptake, lactateExport, complexI, complexIII, complexIV, complexV, protPool, irat), 
              Fluxes = group
              )
  
  pA = ggplot(ds, aes(x = x, y = y, colour = Fluxes, linetype = Fluxes)) +
    geom_line() +
    scale_linetype_manual(values = lst, labels = labels) +
    scale_color_manual(values = cs, labels = labels) +
    #ggplot2::labs(y="Flux (mmol/gDW/h, 1/gDW/h)", x="a", title=title) +
    ggplot2::labs(y=expression("Flux (mmol gDW"^"-1"*"h"^"-1"*", gDW"^"-1"*"h"^"-1"*")"), x="a", title=title) +
    ggplot2::theme_bw() #+ ggplot2::theme(legend.position=legendPos, legend.title = element_blank())
  
  return(pA)

}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

plotFluxesGen = function(RX, title, rxnNames, fluxNames, includeIRat, linestyles, colors, fluxScaling, lineSizes, logA = FALSE, subtractLists = FALSE, hideBiomassUnit = FALSE) {
  col = c(rgb(0,0,0),gg_color_hue(6),rgb(0.5,0.5,0.5))
  if (logA) {
    a = log2(as.numeric(unlist(RX[1,1,1])))
  } else {
    a = as.numeric(unlist(RX[1,1,1]))
  }
  rxns = as.character(unlist(RX[3,1,1]$rxns))
  resultSolutions = RX[2,1,1]$resultSolutions
  nPoints = length(a)
  
  vals = rep(NA, length(rxnNames)*nPoints)
  lowestAInd = NA
  
  for (i in 1:nPoints) {
    if (!is.null(resultSolutions[[i]])) {
      fluxes = as.numeric(unlist(resultSolutions[[i]][[1]][1,1,1]))
      status = as.numeric(unlist(resultSolutions[[i]][[1]][3,1,1]))
      if (status == 1) {
        if (is.na(lowestAInd)) {
          lowestAInd = i;
        }
        for (f in 1:length(rxnNames)) {
          val = 0;
          rxnNamesRecord = rxnNames[[f]]
          for (k in 1:length(rxnNamesRecord)) {
            if (!subtractLists) {
              val = val + fluxes[rxns == rxnNamesRecord[k]]*fluxScaling[f]
            } else {
              if (k == 1) {
                val = val + fluxes[rxns == rxnNamesRecord[k]]*fluxScaling[f]
              } else {
                val = val - fluxes[rxns == rxnNamesRecord[k]]*fluxScaling[f]
              }
            }
          }
          vals[(f-1)*nPoints + i] = val
        }
      }
    }
  }

  #immuno-deficient rat tumor growth
  

  if (includeIRat){
    numFluxes = length(rxnNames) + 1
    labels = c(fluxNames, "exper. growth")
    lst = c(linestyles, 3)
    cs = c(colors, 1)
    irat = c(rep(NA, lowestAInd), rep(0.0080976, length(a)-lowestAInd))
    vals = c(vals, irat)
    lszs = c(lineSizes,1)
  } else {
    numFluxes = length(rxnNames)
    labels = fluxNames
    lst = linestyles
    cs = colors
    lszs = lineSizes
  }  
  group = factor(rep(1:numFluxes, 1, each=length(a)), 1:numFluxes, labels)
  

  ds = tibble(x=rep(a,numFluxes), 
              y=vals, 
              Flux = group)
  
  labls = rep(NA,length(labels))
  for(i in 1:length(labls)) {
    if (labels[i] == "biomass") {
      labls[i] = expression("biomass (h"^"-1"*")")
    } else {
      labls[i] = labels[i]#expression("test")#expression(labels[i])
    }
  }
  
  
  pA = ggplot() +
    geom_line(data=ds, mapping=aes(x = x, y = y, colour = Flux, linetype = Flux, size=Flux)) +
    scale_linetype_manual(values = lst, labels = labls) +
    scale_size_manual(values = lszs, labels = labls) +
    scale_color_manual(values = col[cs], labels = labls)
    if (logA) {
      pA = pA + ggplot2::labs(y=expression("Flux (mmol gDW"^"-1"*"h"^"-1"*", h"^"-1"*")"), x=expression("log"[2] * "(a)"), title=title)
        
    } else {
      if (hideBiomassUnit) {
        pA = pA + ggplot2::labs(y=expression("Flux (mmol gDW"^"-1"*"h"^"-1"*")"), x="a", title=title)
      } else {
        pA = pA + ggplot2::labs(y=expression("Flux (mmol gDW"^"-1"*"h"^"-1"*", h"^"-1"*")"), x="a", title=title)
      }
    }
    #ggplot2::labs(y="Flux (mmol/gDW/h, 1/gDW/h)", x="a", title=title) +
    pA = pA + ggplot2::theme_bw() + theme(legend.text.align = 0)#+ ggplot2::theme(legend.position=legendPos, legend.title = element_blank())
    pA = pA + theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank())
    
#    if (addRects) {
#      d = tibble(x1=0,x2=1,y1=0,y2=0.5);
#      pA = pA + geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2,fill=rgb(1,0,0)), color=rgb(0,0,0), alpha = 0.2)
#      
#    }
    
    
  return(pA)
  
}

plotFluxesGenSigned = function(RX, title, rxnNames, fluxNames, includeIRat, linestyles, colors, fluxScaling, lineSizes, logA = FALSE, subtractLists = FALSE, hideBiomassUnit = FALSE) {
  col = c(rgb(0,0,0),gg_color_hue(6),rgb(0.5,0.5,0.5))
  if (logA) {
    a = log2(as.numeric(unlist(RX[1,1,1])))
  } else {
    a = as.numeric(unlist(RX[1,1,1]))
  }
  rxns = as.character(unlist(RX[3,1,1]$rxns))
  resultSolutions = RX[2,1,1]$resultSolutions
  resultSolutions[[16]]
  nPoints = length(a)
  
  vals = rep(NA, length(rxnNames)*nPoints)
  lowestAInd = NA
  
  for (i in 1:nPoints) {
    if (!is.null(resultSolutions[[i]])) {
      fluxes = as.numeric(unlist(resultSolutions[[i]][[1]][1,1,1]))
      status = as.numeric(unlist(resultSolutions[[i]][[1]][3,1,1]))
      if (status == 1) {
        if (is.na(lowestAInd)) {
          lowestAInd = i;
        }
        for (f in 1:length(rxnNames)) {
          val = 0;
          rxnNamesRecord = rxnNames[[f]]
          for (k in 1:length(rxnNamesRecord)) {
            #sum together output and input, so output becomes negative
            fluxNeg = fluxes[rxns == rxnNamesRecord[k]]*fluxScaling[f]
            fluxPos = fluxes[rxns == paste0(rxnNamesRecord[k],"_REV")]*fluxScaling[f]
            val = val + fluxPos - fluxNeg
          }
          vals[(f-1)*nPoints + i] = val
        }
      }
    }
  }
  
  #immuno-deficient rat tumor growth
  
  
  if (includeIRat){
    numFluxes = length(rxnNames) + 1
    labels = c(fluxNames, "exper. growth")
    lst = c(linestyles, 3)
    cs = c(colors, 1)
    irat = c(rep(NA, lowestAInd), rep(0.0080976, length(a)-lowestAInd))
    vals = c(vals, irat)
    lszs = c(lineSizes,1)
  } else {
    numFluxes = length(rxnNames)
    labels = fluxNames
    lst = linestyles
    cs = colors
    lszs = lineSizes
  }  
  group = factor(rep(1:numFluxes, 1, each=length(a)), 1:numFluxes, labels)
  
  
  ds = tibble(x=rep(a,numFluxes), 
              y=vals, 
              Metabolite = group)
  
  
  pA = ggplot() +
    geom_line(data=ds, mapping=aes(x = x, y = y, colour = Metabolite, linetype = Metabolite, size=Metabolite)) +
    scale_linetype_manual(values = lst, labels = labels) +
    scale_size_manual(values = lszs, labels = labels) +
    scale_color_manual(values = col[cs], labels = labels)
  if (logA) {
    pA = pA + ggplot2::labs(y=expression("Flux (mmol gDW"^"-1"*"h"^"-1"*", h"^"-1"*")"), x=expression("log"[2] * "(a)"), title=title)
    
  } else {
    if (hideBiomassUnit) {
      pA = pA + ggplot2::labs(y=expression("Flux (mmol gDW"^"-1"*"h"^"-1"*")"), x="a", title=title)
    } else {
      pA = pA + ggplot2::labs(y=expression("Flux (mmol gDW"^"-1"*"h"^"-1"*", h"^"-1"*")"), x="a", title=title)
    }
  }
  pA = pA + theme_bw()
  pA = pA + theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank())

  
  return(pA)
  
}

extractRxnFluxes = function(RX, rxn) {
  a = as.numeric(unlist(RX[1,1,1]))
  rxns = as.character(unlist(RX[3,1,1]$rxns))
  resultSolutions = RX[2,1,1]$resultSolutions
  selRxn = rxns == rxn
  nPoints = length(a)
  res = rep(NA, nPoints)
  
  for (i in 1:nPoints) {
    if (!is.null(resultSolutions[[i]])) {
      fluxes = as.numeric(unlist(resultSolutions[[i]][[1]][1,1,1]))
      status = as.numeric(unlist(resultSolutions[[i]][[1]][3,1,1]))
      if (status == 1) {
        res[i] = fluxes[selRxn];
      }
    }
  }
  
  return(res)  
}



