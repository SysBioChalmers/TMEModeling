#code for plotting Fig. 3 and associated suppliementary plots.
##############################################################
library(R.matlab)
library(tidyverse)
library(ggplot2)
library(ggpubr)

fig___path = "Z:/projects/TMEModeling/figures/"

setwd("C:/Work/MatlabCode/projects/TMEModeling/TMEModeling")

source("RCode/FigureHelpers.R")

#plot growth of the small model vs full model, with low ECM

D1_1 = readMat("data/D1_1.mat")
D2_0 = readMat("data/D2_0.mat")
D2_1A = readMat("data/D2_1A.mat")
D2_2A = readMat("data/D2_2A.mat")
D2_3A = readMat("data/D2_3A.mat")
D2_4A = readMat("data/D2_4A.mat")
D2_5A = readMat("data/D2_5A.mat")
D2_8 = readMat("data/D2_8.mat")

aFull = as.numeric(unlist(D2_1A$D2.1A[1,1,1]))
aS = as.numeric(unlist(D1_1$D1.1[1,1,1]))

smGrowth = extractRxnFluxes(D2_0$D2.0, 'MAR13082')
m1Growth = extractRxnFluxes(D2_1A$D2.1A, 'MAR13082')
m2Growth = extractRxnFluxes(D2_2A$D2.2A, 'MAR13082')
m3Growth = extractRxnFluxes(D2_3A$D2.3A, 'MAR13082')
m4Growth = extractRxnFluxes(D2_4A$D2.4A, 'MAR13082')
m5Growth = extractRxnFluxes(D2_5A$D2.5A, 'MAR13082')


numFluxes = 6
labels = c("m0", "m1", "m2", "m3", "m4", "m5")

group = factor(rep(1:numFluxes, 1, each=length(aFull)), 1:numFluxes, labels)
lst = c(2,1,1,1,1,1)
cs = c(1,2,3,4,5,6)

ds = tibble(x=c(aS,rep(aFull,numFluxes-1)), 
            y=c(smGrowth, m1Growth, m2Growth, m3Growth, m4Growth, m5Growth), 
            Model = group
)

pC = ggplot(ds, aes(x = x, y = y, colour = Model, linetype = Model)) +
  geom_line() +
  scale_linetype_manual(values = lst, labels = labels) +
  scale_color_manual(values = cs, labels = labels) +
  labs(y=expression("Growth rate (gDW"^"-1"*"h"^"-1"*")"), x="a", title="Growth rates") +
  theme_bw() + 
  theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank())

pC

ggsave(
  paste0(fig___path, "Fig3C.eps"),
  plot = pC,
  width = 5, height = 3.75, dpi = 300)


D2_6 = readMat("data/D2_6.mat")

succeeded = D2_6$D2.6[1,1,1]
a = as.numeric(unlist(D2_6$D2.6[2,1,1]))


collaborationMets = D2_6$D2.6[3,1,1]$collaborationMets
numCollabMets = colSums(collaborationMets)
#rxns              List,46776     
#metNames          List,21694     
metNames = unlist(D2_6$D2.6[5,1,1]$metNames)

ds = tibble(x=a, 
            y=numCollabMets
)


pD = ggplot(ds, aes(x = x, y = y)) +
  geom_line() +
  ggplot2::labs(y=expression("No. collaboration mets"), x="a", title="Collaboration metabolites") +
  ggplot2::theme_bw() + 
  theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank())
pD

collabMetsLog = rowSums(collaborationMets) > 0
totNumCollabMets = sum(collabMetsLog)
totNumCollabMets#229 - this number can vary slightly due to numerical instabilities

#Export collaboration mets in a table, which is then imported to an excel sheet.
#Metabolite a0 a1 a2 ... a60
#ATP 0 0 0 0 0 0 0 0 1 1 1 1 1 1 ....

#filter the non-collaboration mets

metNamesFilt = metNames[collabMetsLog]
collaborationMetsFilt = collaborationMets[collabMetsLog,]

df = cbind(metNamesFilt, as.data.frame(collaborationMetsFilt))
cn = c("Metabolite", paste0("a", 1:(ncol(df)-1)))
colnames(df) = cn

write.table(df, file = paste0(fig___path, "SupTabCollab.txt"), quote = FALSE, sep = "\t",row.names = FALSE)

#Show oxygen consumption and lactate production in the different cells
#investigate how the "other cells" can contribute - it is possible for them to for example consume lactate, if they are few?
#check what m6 looks like (i.e. D2_8)

cScale = 1/0.9699
oScale = 1/0.03

#other cells: no lipids, little glutamine, no cholesterol, no albumin
pSD = plotFluxesGen(D2_8$D2.8, "Collab. with other cells", 
                   list('MAR09034_REV', 'MAR09048_REV', 'MAR09135', 'o_MAR09034_REV', 'o_MAR09048_REV', 'o_MAR09135'), 
                   c("C glucose upt.", "C oxygen upt.", "C lactate exp.", "O glucose upt.", "O oxygen upt.", "O lactate exp."), 
                   FALSE, 
                   c(1,1,1,2,2,2), 
                   c(1,6,7,1,6,7),
                   fluxScaling = c(cScale,cScale*10,cScale,oScale,oScale*10,oScale),
                   lineSizes = c(1,1,1,1,1,1),
                   hideBiomassUnit = TRUE)
pSD


################
#macrophages
D2_9 = readMat("data/D2_9.mat")
D2_10 = readMat("data/D2_10.mat")


normGrowth = extractRxnFluxes(D2_10$D2.10, 'MAR13082')
macroGrowth = extractRxnFluxes(D2_9$D2.9, 'MAR13082')
plot(normGrowth, macroGrowth)
as2 = as.numeric(unlist(D2_9$D2.9[1,1,1]))

ds = tibble(x=as2, 
            y=macroGrowth/normGrowth
)

pSE = ggplot(ds, aes(x = x, y = y)) +
  geom_line() +
  ggplot2::labs(y=expression("Growth ratio: collab. vs no collab."), x="a", title="Macrophage collaboration") +
  ggplot2::theme_bw() + 
  theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank())

pSE



#Limiting metabolites

D2_3B = readMat("data/D2_3B.mat")
pSA = plotFluxDiv2(D2_3B$D2.3B, "Growth limitation from metabolites");
pSA


#Growth advantage of literature collaboration metabolites
D2_7 = readMat("data/D2_7.mat") #m1 with blocked collaboration
D2_11 = readMat("data/D2_11.mat") #m1 with literature collaboration mets only

normGrowth2 = extractRxnFluxes(D2_7$D2.7, 'MAR13082') #compare to the s model
literatureGrowth = extractRxnFluxes(D2_11$D2.11, 'MAR13082')
as3 = as.numeric(unlist(D2_7$D2.7[1,1,1]))

#for some reason there is an NA in the middle somewhere, don't display that point
diff = literatureGrowth - normGrowth2;
sel = rep(TRUE, length(diff))
sel[(1:length(sel)) > 30 & is.na(diff)] = FALSE

ds = tibble(x=as3[sel], 
            y=diff[sel]
)

formatC(ds$y[90], digits = 30, format = "f") #1.000000000122490684262288596074, maybe just roundoff uncertainties, really nothing

pSB = ggplot(ds, aes(x = x, y = y)) +
  geom_line() +
  #ylim(1,1.0015) +
  labs(y=expression("Growth increase (h"^"-1"*")"), x="a", title="Effect of literature collab. metabolites") +
  theme_bw() + 
  theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank())

pSB
#as can be seen, there is really no benefit at all, just roundoff noise - no point showing a plot
D2_11b = readMat("data/D2_11b.mat") #m1 with literature collaboration mets only, H2O2 added.
literatureGrowth2 = extractRxnFluxes(D2_11b$D2.11b, 'MAR13082')
diff2 = literatureGrowth2 - normGrowth2;
sel = rep(TRUE, length(diff2))
sel[(1:length(sel)) > 30 & is.na(diff2)] = FALSE

ds = tibble(x=as3[sel], 
            y=diff2[sel]
)

formatC(ds$y[90], digits = 30, format = "f") #0.000000000000000097144514654701, maybe just roundoff uncertainties, really nothing

pSC = ggplot(ds, aes(x = x, y = y)) +
  geom_line() +
  #ylim(1,1.0015) +
  labs(y=expression("Growth increase (h"^"-1"*")"), x="a", title="Effect of lit. collab. metabolites + H2O2") +
  theme_bw() + 
  theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank())

pSC



D1_8 = readMat("data/D1_8.mat")

#plot to test only
plotFluxesGen(D1_8$D1.8, "Metabolite uptake and export", 
              #list('biomass_human', 'HMR_9034_REV', 'HMR_9063_REV', 'HMR_9033_REV', 'HMR_9285_REV', 'HMR_9048_REV', 'HMR_9135', 'HMR_9068'),  #'EX_propropro[e]'
              list('MAR13082', 'MAR09034_REV', 'MAR09063_REV', 'MAR09033_REV', 'MAR09285_REV', 'MAR09048_REV', 'MAR09135', 'MAR09068'),  #'EX_propropro[e]'
              c("biomass", "glucose upt.", "glutamine upt.", "lipid pool upt.", "cholesterol upt.", "oxygen upt.", "lactate exp.", 'proline exp.'), 
              TRUE, 
              c(1,1,1,1,1,1,1,1), 
              c(1,2,3,4,8,6,7,5),
              fluxScaling = c(1,0.01, 0.05, 0.2, 10, 0.1, 0.01, 0.05),
              lineSizes = c(1.3,1,1,1,1,1,1,1))
#lactate export is half of glucose uptake, ok

normGrowth = extractRxnFluxes(D1_1$D1.1, 'MAR13082')
lacRedGrowth = extractRxnFluxes(D1_8$D1.8, 'MAR13082')


as2 = as.numeric(unlist(D1_1$D1.1[1,1,1]))

ds = tibble(x=rep(as2,2), y=c(normGrowth, lacRedGrowth), Setup = factor(c(rep(0,length(as2)), rep(1,length(as2))), c(0,1), c("normal", "constr. lactate")))

pSF = ggplot(ds, aes(x = x, y = y, colour=Setup)) +
  geom_line(size=1) +
  scale_color_manual(values = c(rgb(0,0,0),gg_color_hue(6)[1]), labels = c("normal", "constr. lactate")) +
  #ylim(1,1.0015) +
  labs(y=expression("Growth rate (h"^"-1"*")"), x="a", title="") +
  theme_bw() + 
  theme(panel.background = element_rect("white", "white", 0, 0, "white"), panel.grid.major= element_blank(),panel.grid.minor= element_blank())

pSF



ggsave(
  paste0(fig___path, "Fig3D.eps"),
  plot = pD,
  width = 5, height = 3.75, dpi = 300)


figSupCollab = ggarrange(pSA,pSB,pSC,pSD,pSE, nrow=3, ncol=2, labels=c("A","B","C","D","E"))

ggsave(
  paste0(fig___path, "FigSupCollab.png"),
  plot = figSupCollab,
  width = 10, height = 10, dpi = 300)

ggsave(
  paste0(fig___path, "FigSupLact.eps"),
  plot = pSF,
  width = 5, height = 2.75, dpi = 300)





