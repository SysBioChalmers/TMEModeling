# Code for generating Fig. S2
##########################################
library(R.matlab)
library(tidyverse)
library(ggplot2)
library(ggpubr)

fig___path = "Z:/projects/TMEModeling/figures/"

setwd("C:/Work/MatlabCode/projects/TMEModeling/TMEModeling")

DCData = readMat("data/diffCoeff.mat")

x = as.vector(DCData$plotExp[1,1,1]$x)
y = as.vector(DCData$plotExp[2,1,1]$y)
mets = unlist(DCData$plotExp[3,1,1]$mets)
vjust = rep(-1,16)
vjust[mets == "threonine"] = 1.5


ds = tibble(x=x, y=y)

pA = ggplot(ds,aes(x=x,y=y)) +
  stat_smooth(method = "lm", col = "red") +
  geom_point() +
  geom_text(mapping = aes(label=mets), hjust=0.5, vjust=vjust) +
  labs(y=expression("Diff. Coeff * 10"^10*" ( m"^2*"/s)"), x="Molecular weight (g/mol)") +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pA


ggsave(
  paste0(fig___path, "FigSupDiffConst.png"),
  plot = pA,
  width = 6, height = 6, dpi = 300)




