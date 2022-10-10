library(R.matlab)
library(tidyverse)
library(ggplot2)
library(ggpubr)

fig___path = "Z:/projects/TMEModeling/figures/"

setwd("C:/Work/MatlabCode/projects/TMEModeling/TMEModeling")

source("RCode/FigureHelpers.R")

##############################################
#Supplementary fig for the blood flow model - single model
##############################################

##############
# Fig A
##############

D1_1_bf = readMat("data/D1_1_bloodflow.mat")

pA = plotFluxesGen(D1_1_bf$D1.1.bloodflow, "Metabolite uptake and export", 
                   list('MAR13082', 'MAR09034_REV', 'MAR09063_REV', 'MAR09033_REV', 'MAR09285_REV', 'MAR09048_REV', 'MAR09135'),  #'EX_propropro[e]'
                   c("biomass", "glucose upt.", "glutamine upt.", "lipid pool upt.", "cholesterol upt.", "oxygen upt.", "lactate exp."), 
                   FALSE, 
                   c(1,1,1,1,1,1,1,1), 
                   c(1,2,3,4,8,6,7,5),
                   fluxScaling = c(1,0.01, 0.05, 0.2, 10, 0.1, 0.01, 0.05),
                   lineSizes = c(1.3,1,1,1,1,1,1,1),
                   hideBiomassUnit = TRUE)
pA


##############
# Fig B
##############

D1_2_bf = readMat("data/D1_2_bloodflow.mat");
pB = plotFluxDiv(D1_2_bf$D1.2.bloodflow, "Growth limitation from metabolites");
pB

##############
# Fig C
##############
pC = plotRed(D1_2_bf$D1.2.bloodflow, "");
pC


##############
# Fig D
##############


D1_1_bf = readMat("data/D1_1_bloodflow.mat")
D1_9_bf = readMat("data/D1_9_bloodflow.mat")
D1_10_bf = readMat("data/D1_10_bloodflow.mat")
D1_11_bf = readMat("data/D1_11_bloodflow.mat")
D1_12_bf = readMat("data/D1_12_bloodflow.mat")
D1_13_bf = readMat("data/D1_13_bloodflow.mat")
D1_14_bf = readMat("data/D1_14_bloodflow.mat")
D1_15_bf = readMat("data/D1_15_bloodflow.mat")
D1_16_bf = readMat("data/D1_16_bloodflow.mat")
D1_17_bf = readMat("data/D1_17_bloodflow.mat")
D1_18_bf = readMat("data/D1_18_bloodflow.mat")
D1_22_bf = readMat("data/D1_22_bloodflow.mat")


fB = extractRxnFluxes(D1_1_bf$D1.1.bloodflow, 'MAR13082') #biomass_human
noATP = extractRxnFluxes(D1_9_bf$D1.9.bloodflow, 'incomplete_biomass')
noNucl = extractRxnFluxes(D1_10_bf$D1.10.bloodflow, 'incomplete_biomass')
noGlycogen = extractRxnFluxes(D1_11_bf$D1.11.bloodflow, 'incomplete_biomass')
noCofact = extractRxnFluxes(D1_12_bf$D1.12.bloodflow, 'incomplete_biomass')
noProt = extractRxnFluxes(D1_13_bf$D1.13.bloodflow, 'incomplete_biomass')
noLip = extractRxnFluxes(D1_14_bf$D1.14.bloodflow, 'incomplete_biomass')
noMet = extractRxnFluxes(D1_15_bf$D1.15.bloodflow, 'incomplete_biomass')
noATPProt = extractRxnFluxes(D1_16_bf$D1.16.bloodflow, 'incomplete_biomass')
noATP2x = extractRxnFluxes(D1_17_bf$D1.17.bloodflow, 'incomplete_biomass')
noATPNoProt = extractRxnFluxes(D1_18_bf$D1.18.bloodflow, 'incomplete_biomass')
no22 = extractRxnFluxes(D1_22_bf$D1.22.bloodflow, 'incomplete_biomass')

aS5 = as.numeric(unlist(D1_1_bf$D1.1.bloodflow[1,1,1]))
numFluxes = 10
labels = c("No ATP", "No nucleotides", "No glycogen", "No cofactors", "No protein", "No lipids", "No metabolites", "No ATP prot.", "No 2xATP", "No 2xATP, lipids")

group = factor(rep(1:numFluxes, each = length(aS5)), 1:numFluxes, labels)
lst = c(1,1,1,1,1,1,1,2,2,2)
cs = c(1,2,3,4,5,6,8,6,1,2)

ds = tibble(x=rep(aS5,numFluxes), 
            y=log2(c(noATP/fB, noNucl/fB, noGlycogen/fB, noCofact/fB, noProt/fB, noLip/fB, noMet/fB, noATPProt/fB, noATP2x/fB, no22/fB)), 
            Removed = group
)

pD = ggplot(ds, aes(x = x, y = y, colour = Removed, linetype = Removed)) +
  geom_line() +
  scale_linetype_manual(values = lst, labels = labels) +
  scale_color_manual(values = cs, labels = labels) +
  labs(y=expression(Log[2]*"(Growth ratio vs full biomass)"), x="a") +
  theme_bw() + 
  theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pD

pSupBf1 = ggarrange(pA,pB,pC,pD, nrow=2, ncol=2, labels=c("A","B","C","D"))

ggsave(
  paste0(fig___path, "SupBf1.png"),
  plot = pSupBf1,
  width = 10, height = 7, dpi = 300)



##############################################
#Supplementary fig for the blood flow model - combined model
##############################################



##################################################################
# Fig A - Growth advantage of literature collaboration metabolites
##################################################################

D2_7_bf = readMat("data/D2_7_bloodflow.mat") #m1 with blocked collaboration
D2_11_bf = readMat("data/D2_11_bloodflow.mat") #m1 with literature collaboration mets only

normGrowth2 = extractRxnFluxes(D2_7_bf$D2.7.bloodflow, 'MAR13082')
literatureGrowth = extractRxnFluxes(D2_11_bf$D2.11.bloodflow, 'MAR13082')
as3 = as.numeric(unlist(D2_7_bf$D2.7.bloodflow[1,1,1]))

max(literatureGrowth/normGrowth2, na.rm = TRUE) #1.008611, so 0.86% - becomes that large because the growthrate is close to 0. not that relevant

#for some reason there is an NA in the middle somewhere, don't display that point
diff = literatureGrowth - normGrowth2;
sel = rep(TRUE, length(diff))
sel[(1:length(sel)) > 30 & is.na(diff)] = FALSE

ds = tibble(x=as3[sel], 
            y=diff[sel]
)

formatC(ds$y[90], digits = 30, format = "f") #1.000000000122490684262288596074, maybe just roundoff uncertainties, really nothing

pA2 = ggplot(ds, aes(x = x, y = y)) +
  geom_line() +
  labs(y=expression("Growth increase (h"^"-1"*")"), x="a", title="Effect of literature collab. metabolites") +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pA2
#as can be seen, there is really no benefit at all, just roundoff noise - no point showing a plot


#####################
# Fig B - macrophages
#####################

D2_9_bf = readMat("data/D2_9_bloodflow.mat")
D2_10_bf = readMat("data/D2_10_bloodflow.mat")


normGrowth = extractRxnFluxes(D2_10_bf$D2.10.bloodflow, 'MAR13082')
macroGrowth = extractRxnFluxes(D2_9_bf$D2.9.bloodflow, 'MAR13082')
plot(normGrowth, macroGrowth)
as2 = as.numeric(unlist(D2_9_bf$D2.9.bloodflow[1,1,1]))

ds = tibble(x=as2, 
            y=macroGrowth/normGrowth
)

pB2 = ggplot(ds, aes(x = x, y = y)) +
  geom_line() +
  labs(y=expression("Growth ratio: collab. vs no collab."), x="a", title="Macrophage collaboration") +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pB2

pSupBf2 = ggarrange(pA2,pB2, nrow=1, ncol=2, labels=c("A","B"))

ggsave(
  paste0(fig___path, "SupBf2.png"),
  plot = pSupBf2,
  width = 8, height = 4, dpi = 300)
