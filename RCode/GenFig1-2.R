library(R.matlab)
library(tidyverse)
library(ggplot2)
library(ggpubr)

fig___path = "Z:/projects/TMEModeling/figures/"

setwd("C:/Work/MatlabCode/projects/TMEModeling/TMEModeling")

source("RCode/FigureHelpers.R")

D1_1 = readMat("data/D1_1.mat")

pA = plotFluxesGen(D1_1$D1.1, "Metabolite uptake and export", 
                   list('MAR13082', 'MAR09034_REV', 'MAR09063_REV', 'MAR09033_REV', 'MAR09285_REV', 'MAR09048_REV', 'MAR09135'),  #'EX_propropro[e]'
                   c("biomass", "glucose upt.", "glutamine upt.", "lipid pool upt.", "cholesterol upt.", "oxygen upt.", "lactate exp."), 
                   FALSE, 
                   c(1,1,1,1,1,1,1,1), 
                   c(1,2,3,4,8,6,7,5),
                   fluxScaling = c(1,0.01, 0.05, 0.2, 10, 0.1, 0.01, 0.05),
                   lineSizes = c(1.3,1,1,1,1,1,1,1),
                   hideBiomassUnit = TRUE)
pA = pA + geom_vline(xintercept = 0.000020, linetype="dashed")
pA = pA + geom_vline(xintercept = 0.000068, linetype="dashed")
pA = pA + annotate(geom="text", x=0.000007, y=0.045, label="Necrotic", color="black", angle = 90, size=5)
pA = pA + annotate(geom="text", x=0.000045, y=0.075, label="Hypoxic", color="black", size=5)
pA = pA + annotate(geom="text", x=0.000085, y=0.050, label="Enzyme-\nlimited", color="black", size=5)
pA

#pB = plotFluxesGen(D1_1$D1.1, "Internal fluxes", 
#                   list('MAR13082', 'MAR04358', 'MAR06921', 'MAR06918', 'MAR06914', 'MAR06916'), 
#                   c("biomass", "glycolysis        ", "complex I", "complex III", "complex IV", "complex V"), 
#                   FALSE, 
#                   c(1,1,1,1,1,1), 
#                   c(1,2,3,4,5,6),
#                   fluxScaling = c(1,0.01, rep(0.1,4)),
#                   lineSizes = c(1.3,1,1,1,1,1),
#                   hideBiomassUnit = TRUE)
#pB

pB = plotFluxesGen(D1_1$D1.1, "Internal fluxes", 
                   list('MAR13082', 'MAR04358', 'MAR06921', c('MAR04652_REV','MAR08743_REV'),'MAR06918', 'MAR06914', 'MAR06916'), 
                   c("biomass", "glycolysis        ", "complex I", "complex II", "complex III", "complex IV", "complex V"), 
                   FALSE, 
                   c(1,1,1,1,1,1,1), 
                   c(1,2,3,4,5,6,7),
                   fluxScaling = c(1,0.01, rep(0.1,5)),
                   lineSizes = c(1.3,1,1,1,1,1,1),
                   hideBiomassUnit = TRUE)
pB

D1_2 = readMat("data/D1_2.mat");
pC = plotFluxDiv(D1_2$D1.2, "Growth limitation from metabolites");
pC


D1_3 = readMat("data/D1_3.mat");
pD = plotFluxesGen(D1_3$D1.3, "Oxygen unconstrained", 
                   #list('biomass_human', 'HMR_9034_REV', 'HMR_9063_REV', 'HMR_9033_REV', 'HMR_9285_REV', 'HMR_9048_REV', 'HMR_9135', 'HMR_9730_REV'), 
                   list('MAR13082', 'MAR09034_REV', 'MAR09063_REV', 'MAR09033_REV', 'MAR09285_REV', 'MAR09048_REV', 'MAR09135', 'MAR09730_REV'), 
                   c("biomass", "glucose upt.", "glutamine upt.", "lipid pool upt.", "cholesterol upt.", "oxygen upt.", "lactate exp.", "albumin upt."), 
                   FALSE, 
                   c(1,1,1,1,1,1,1,2), 
                   c(1,2,3,4,8,6,7,3),
                   fluxScaling = c(1,0.01, 0.05, 0.2, 10, 0.1, 0.01, 200),
                   lineSizes = c(1.3,1,1,1,1,1,1,1),
                   hideBiomassUnit = TRUE)
pD

ggsave(
  paste0(fig___path, "Fig1b.eps"),
  plot = pA,
  width = 5, height = 3.5, dpi = 300)

ggsave(
  paste0(fig___path, "Fig1c.eps"),
  plot = pB,
  width = 5, height = 3.5, dpi = 300)

#pA2 = pA + ggplot2::theme(legend.position="bottom", legend.title = element_blank())
#pA2 #This might be an alternative

ggsave(
  paste0(fig___path, "FigSupFVA1.png"),
  plot = pC + ggtitle(""),
  width = 5, height = 3.5, dpi = 300)

ggsave(
  paste0(fig___path, "FigSupFullO2.png"),
  plot = pD + ggtitle(""),
  width = 5, height = 3.5, dpi = 300)


#################################################################
#supplementary fig for growth reduction when reducing a metabolite
#################################################################

D1_2 = readMat("data/D1_2.mat");
pSRed = plotRed(D1_2$D1.2, "");
pSRed

ggsave(
  paste0(fig___path, "FigSupRedMet.png"),
  plot = pSRed,
  width = 5, height = 3.5, dpi = 300)

#################################################################
#supplementary fig for biomass equation with reduced ATP cost
#################################################################

D1_2B = readMat("data/D1_2B.mat");
pSFVA2 = plotFluxDiv(D1_2B$D1.2B, "Growth limitation, 0.5 ATP");
pSFVA2
pSRed2 = plotRed(D1_2B$D1.2B, "Met. reduction, 0.5 ATP");
pSRed2
D1_2C = readMat("data/D1_2C.mat");
pSGrowth2 = plotFluxesGen(D1_2C$D1.2C, "Fluxes, 0.5 ATP", 
                          list('MAR13082', 'MAR09034_REV', 'MAR09063_REV', 'MAR09033_REV', 'MAR09285_REV', 'MAR09048_REV', 'MAR09135'),  #'EX_propropro[e]'
                          c("biomass", "glucose upt.", "glutamine upt.", "lipid pool upt.", "cholesterol upt.", "oxygen upt.", "lactate exp."), 
                          FALSE, 
                          c(1,1,1,1,1,1,1,1), 
                          c(1,2,3,4,8,6,7,5),
                          fluxScaling = c(1,0.01, 0.05, 0.2, 10, 0.1, 0.01, 0.05),
                          lineSizes = c(1.3,1,1,1,1,1,1,1),
                          hideBiomassUnit = TRUE)

pSGrowth2

D1_2D = readMat("data/D1_2D.mat");
pSFVA3 = plotFluxDiv(D1_2D$D1.2D, "Growth limitation, 0.25 ATP");
pSFVA3
pSRed3 = plotRed(D1_2D$D1.2D, "Met. reduction, 0.25 ATP");
pSRed3
D1_2E = readMat("data/D1_2E.mat");
pSGrowth3 = plotFluxesGen(D1_2E$D1.2E, "Fluxes, 0.25 ATP", 
                          list('MAR13082', 'MAR09034_REV', 'MAR09063_REV', 'MAR09033_REV', 'MAR09285_REV', 'MAR09048_REV', 'MAR09135'),  #'EX_propropro[e]'
                          c("biomass", "glucose upt.", "glutamine upt.", "lipid pool upt.", "cholesterol upt.", "oxygen upt.", "lactate exp."), 
                          FALSE, 
                          c(1,1,1,1,1,1,1,1), 
                          c(1,2,3,4,8,6,7,5),
                          fluxScaling = c(1,0.01, 0.05, 0.2, 10, 0.1, 0.01, 0.05),
                          lineSizes = c(1.3,1,1,1,1,1,1,1),
                          hideBiomassUnit = TRUE)

pSGrowth3


figSupRedATP = ggarrange(pSGrowth2,pSGrowth3,pSFVA2,pSFVA3,pSRed2,pSRed3, nrow=3, ncol=2, labels=c("A","B","C","D","E","F"))
figSupRedATP




ggsave(
  paste0(fig___path, "FigSupRedATP.png"),
  plot = figSupRedATP,
  width = 10, height = 10.5, dpi = 300)

#################################################################
#supplementary fig for which processes are most limiting for growth
#################################################################

#A figure showing which parts of the biomass function that are limiting for growth

D1_1 = readMat("data/D1_1.mat")
D1_9 = readMat("data/D1_9.mat")
D1_10 = readMat("data/D1_10.mat")
D1_11 = readMat("data/D1_11.mat")
D1_12 = readMat("data/D1_12.mat")
D1_13 = readMat("data/D1_13.mat")
D1_14 = readMat("data/D1_14.mat")
D1_15 = readMat("data/D1_15.mat")
D1_16 = readMat("data/D1_16.mat")
D1_17 = readMat("data/D1_17.mat")
D1_18 = readMat("data/D1_18.mat")
D1_22 = readMat("data/D1_22.mat")


fB = extractRxnFluxes(D1_1$D1.1, 'MAR13082') #biomass_human
noATP = extractRxnFluxes(D1_9$D1.9, 'incomplete_biomass')
noNucl = extractRxnFluxes(D1_10$D1.10, 'incomplete_biomass')
noGlycogen = extractRxnFluxes(D1_11$D1.11, 'incomplete_biomass')
noCofact = extractRxnFluxes(D1_12$D1.12, 'incomplete_biomass')
noProt = extractRxnFluxes(D1_13$D1.13, 'incomplete_biomass')
noLip = extractRxnFluxes(D1_14$D1.14, 'incomplete_biomass')
noMet = extractRxnFluxes(D1_15$D1.15, 'incomplete_biomass')
noATPProt = extractRxnFluxes(D1_16$D1.16, 'incomplete_biomass')
noATP2x = extractRxnFluxes(D1_17$D1.17, 'incomplete_biomass')
noATPNoProt = extractRxnFluxes(D1_18$D1.18, 'incomplete_biomass')
no22 = extractRxnFluxes(D1_22$D1.22, 'incomplete_biomass')

aS5 = as.numeric(unlist(D1_1$D1.1[1,1,1]))
numFluxes = 10
labels = c("No ATP", "No nucleotides", "No glycogen", "No cofactors", "No protein", "No lipids", "No metabolites", "No ATP prot.", "No 2xATP", "No 2xATP, lipids")

group = factor(rep(1:numFluxes, each = length(aS5)), 1:numFluxes, labels)
lst = c(1,1,1,1,1,1,1,2,2,2)
cs = c(1,2,3,4,5,6,8,6,1,2)

ds = tibble(x=rep(aS5,numFluxes), 
            y=log2(c(noATP/fB, noNucl/fB, noGlycogen/fB, noCofact/fB, noProt/fB, noLip/fB, noMet/fB, noATPProt/fB, noATP2x/fB, no22/fB)), 
            Removed = group
)

pX = ggplot(ds, aes(x = x, y = y, colour = Removed, linetype = Removed)) +
  geom_line() +
  scale_linetype_manual(values = lst, labels = labels) +
  scale_color_manual(values = cs, labels = labels) +
  ggplot2::labs(y=expression(Log[2]*"(Growth ratio vs full biomass)"), x="a") +
  ggplot2::theme_bw() + ggplot2::theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pX

ggsave(
  paste0(fig___path, "SupReducedBiomass.png"),
  plot = pX,
  width = 5, height = 4, dpi = 300)


#############################
# Amino acid figures
#############################

D1_1 = readMat("data/D1_1.mat")

#the amino acids with high fluxes
p2A = plotFluxesGenSigned(D1_1$D1.1, "", 
                          list("MAR09070","MAR09063","MAR09067","MAR09068","MAR09069","MAR09044"),
                          c("aspartate","glutamine","glycine","proline",   "serine","threonine"), 
                          FALSE, 
                          c(1,1,1,1,1,1), 
                          c(1,2,3,4,5,6),
                          fluxScaling = c(1,1,1,1,1,1),
                          lineSizes = c(1,1,1,1,1,1),
                          hideBiomassUnit = TRUE)
p2A



#the amino acids with low fluxes and that are not proportional to growth
p2Sup1 = plotFluxesGenSigned(D1_1$D1.1, "", 
                          list("MAR09066","MAR09062",  "MAR09065", "MAR09071",  "MAR09038", "MAR09046"),
                          c(   "arginine","asparagine","cysteine", "glutamate", "histidine","valine"), 
                          FALSE, 
                          c(1,1,1,1,1,1), 
                          c(1,2,3,4,5,6),
                          fluxScaling = c(1,1,1,1,1,1),
                          lineSizes = c(1,1,1,1,1,1),
                          hideBiomassUnit = TRUE)
p2Sup1



#the rest of the amino acids (proportional to growth)
p2Sup2 = plotFluxesGenSigned(D1_1$D1.1, "", 
                          list("MAR09061","MAR09039",  "MAR09040","MAR09041","MAR09042",  "MAR09043",     "MAR09045",  "MAR09064"),
                          c(   "alanine", "isoleucine","leucine", "lysine",  "methionine","phenylalanine","tryptophan","tyrosine"), 
                          FALSE, 
                          c(1,1,1,1,2,2,2,2), 
                          c(1,2,3,4,1,2,3,4),
                          fluxScaling = c(1,1,1,1,1,1,1,1),
                          lineSizes = c(1,1,1,1,1,1,1,1),
                          hideBiomassUnit = TRUE)
p2Sup2

p2Sup = ggarrange(p2Sup1,p2Sup2, nrow=1, ncol=2, labels=c("A","B"))


ggsave(
  paste0(fig___path, "Fig2A.eps"),
  plot = p2A,
  width = 5, height = 4, dpi = 300)

ggsave(
  paste0(fig___path, "Fig2Sup.png"),
  plot = p2Sup,
  width = 10, height = 4, dpi = 300)

##############################################
#Supplementary for glutamine addiction
##############################################
D3_1 = readMat("data/D3_1.mat")$s[,,1]
mets = as.character(unlist(D3_1$mets))
atpHyp = as.numeric(D3_1$ATPLowO2NoProdh)
tcaHyp = as.numeric(D3_1$TCACycleFluxesLowO2NoProdh)
c1Hyp = as.numeric(D3_1$complexILowO2NoProdh)
c3Hyp = as.numeric(D3_1$complexIIILowO2NoProdh)
c5Hyp = as.numeric(D3_1$complexVLowO2NoProdh)
atpEnz = as.numeric(D3_1$ATPEnzLim)
tcaEnz = as.numeric(D3_1$TCACycleFluxesEnzLim)
c1Enz = as.numeric(D3_1$complexIEnzLim)
c3Enz = as.numeric(D3_1$complexIIIEnzLim)
c5Enz = as.numeric(D3_1$complexVEnzLim)

#compare lactate and glutamine
#

y = c(atpHyp[1],atpHyp[3],
      tcaHyp[3],tcaHyp[3],
      c1Hyp[1],c1Hyp[3],
      c3Hyp[1],c3Hyp[3],
      c5Hyp[1],c5Hyp[3])
Condition = factor(rep(c(1,2),5), c(1,2), c("Lactate","Glutamine"))
x = factor(rep(c(1,2,3,4,5), 1, each=2), c(1,2,3,4,5), c("Total ATP", "TCA Cycle","Complex I", "Complex III", "Complex V"))

dfPlot = data.frame(x, y, Condition)


pSBPA = ggplot(data=dfPlot, aes(x=x, y=y, fill=Condition)) +
  geom_bar(stat="identity",position=position_dodge()) +
  labs( y=expression("Flux (mmol gDW"^"-1"*"h"^"-1"*")"), x="", fill="Substrate: ", title = "Hypoxic") +
  scale_fill_manual(values=c("#82E182", "#229A22", "#E47060", "#AC210E")) + 
  theme_bw() +
  theme(axis.text = element_text(size = 8), legend.position="bottom")

pSBPA

y = c(atpEnz[1],atpEnz[3],
      tcaEnz[1],tcaEnz[3],
      c1Enz[1],c1Enz[3],
      c3Enz[1],c3Enz[3],
      c5Enz[1],c5Enz[3])
Condition = factor(rep(c(1,2),5), c(1,2), c("Lactate","Glutamine"))
x = factor(rep(c(1,2,3,4,5), 1, each=2), c(1,2,3,4,5), c("Total ATP", "TCA Cycle","Complex I", "Complex III", "Complex V"))

dfPlot = data.frame(x, y, Condition)


pSBPB = ggplot(data=dfPlot, aes(x=x, y=y, fill=Condition)) +
  geom_bar(stat="identity",position=position_dodge()) +
  labs( y=expression("Flux (mmol gDW"^"-1"*"h"^"-1"*")"), x="", fill="Substrate: ", title = "Enzyme-limited") +
  scale_fill_manual(values=c("#82E182", "#229A22", "#E47060", "#AC210E")) +
  theme_bw() +
  theme(axis.text = element_text(size = 8), legend.position="bottom")

pSBPB


pSX = ggarrange(pSBPA,pSBPB, nrow=1, ncol=2, labels=c("A","B"))
pSX
ggsave(
  paste0(fig___path, "FigSupAA.png"),
  plot = pSX,
  width = 10, height = 6, dpi = 300)




##############################################
#Supplementary for PRODH/PYCR
##############################################


pX = plotFluxesGen(D1_1$D1.1, "", 
                   list('MAR13082', 'MAR03838', 'MAR03835'), 
                   c("biomass", "PRODH", "PYCR3"), 
                   FALSE, 
                   c(1,1,1), 
                   c(1,2,3),
                   fluxScaling = c(1,1, 1),
                   lineSizes = c(1.3,1,1),
                   hideBiomassUnit = TRUE)
pX

ggsave(
  paste0(fig___path, "SupProline.png"),
  plot = pX,
  width = 5, height = 3.5, dpi = 300)


#####################################
# Supplementary figure about ATP production 
# when some reactions are allowed to go in reverse
#####################################

#First increase in ATP production

D5_6_1 = readMat("data/D5_6_1.mat")
D5_6_2 = readMat("data/D5_6_2.mat")
D5_6_3 = readMat("data/D5_6_3.mat")
D5_6_4 = readMat("data/D5_6_4.mat")
D5_6_5 = readMat("data/D5_6_5.mat")

a5 = as.numeric(unlist(D5_6_1$D5.6.1[1,1,1]))

b1 = extractRxnFluxes(D5_6_1$D5.6.1, 'MAR03964')
b2 = extractRxnFluxes(D5_6_2$D5.6.2, 'MAR03964')
b3 = extractRxnFluxes(D5_6_3$D5.6.3, 'MAR03964')
b4 = extractRxnFluxes(D5_6_4$D5.6.4, 'MAR03964')
b5 = extractRxnFluxes(D5_6_5$D5.6.5, 'MAR03964')

c2 = (b2/b1 - 1)*100
c3 = (b3/b1 - 1)*100
c4 = (b4/b1 - 1)*100
c5 = (b5/b1 - 1)*100


#We skip the Complex I in reverse with blocked citrate synthase - it has a very limited effect and makes the plot less readable

numFluxes = 3
labels = c("PRODH", "C II",  "Both")

group = factor(rep(1:numFluxes, 1, each=length(a5)), 1:numFluxes, labels)
lst = c(1,1,2)
cs = c(1,2,4)

ds = tibble(x=rep(a5,numFluxes), 
            y=c(c2,c4,c5), 
            Reaction = group
)

pSA = ggplot(ds, aes(x = x, y = y, colour = Reaction, linetype = Reaction)) +
  geom_line(size=1.3) +
  scale_linetype_manual(values = lst, labels = labels) +
  scale_color_manual(values = cs, labels = labels) +
  theme_bw() + 
  theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  labs(y=expression("ATP incr.(%)"), x="a", title="")
pSA

ggsave(
  paste0(fig___path, "Fig3E.eps"),
  plot = pSA,
  width = 5, height = 2, dpi = 300)


#######################################
# Now plot succinate export and where that comes from.
# #####################################

dSuccExp = extractRxnFluxes(D5_6_4$D5.6.4, 'succExp')
dCIIRev1 = extractRxnFluxes(D5_6_4$D5.6.4, 'MAR04652') #the reaction is defined in reverse direction, version 1
dCIIRev2 = extractRxnFluxes(D5_6_4$D5.6.4, 'MAR08743') #the reaction is defined in reverse direction, version 2
dCIIRev = dCIIRev1 + dCIIRev2
dATPGen = extractRxnFluxes(D5_6_4$D5.6.4, 'MAR04152') #there are two reactions, one that generates GTP, catch them both
dGTPGen = extractRxnFluxes(D5_6_4$D5.6.4, 'MAR04147_REV') #Defined in the opposite order for some reason
dForwardDir = dATPGen + dGTPGen

numFluxes = 3
labels = c("C II rev",  "Fwd TCA", "Succ. Exp.")

group = factor(rep(1:numFluxes, 1, each=length(a5)), 1:numFluxes, labels)
lst = c(1,1,2)
cs = c(2,3,1)

ds = tibble(x=rep(a5,numFluxes), 
            y=c(dCIIRev,dForwardDir,dSuccExp), 
            Reaction = group
)

pSB = ggplot(ds, aes(x = x, y = y, colour = Reaction, linetype = Reaction)) +
  geom_line(size=1.3) +
  scale_linetype_manual(values = lst, labels = labels) +
  scale_color_manual(values = cs, labels = labels) +
  theme_bw() + 
  theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank()) +
  labs(y=expression("Flux (mmol gDW"^"-1"*"h"^"-1"*")"), x="a", title="")

pSB

ggsave(
  paste0(fig___path, "FigSupCII.png"),
  plot = pSB,
  width = 5, height = 3.75, dpi = 300)




