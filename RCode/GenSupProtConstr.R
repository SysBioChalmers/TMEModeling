#This figure is not used in the paper, but can be useful when experimenting with fitting the protein constraint.

library(R.matlab)
library(tidyverse)
library(ggplot2)
library(ggpubr)

fig___path = "Z:/projects/TMEModeling/figures/"

setwd("C:/Work/MatlabCode/projects/TMEModeling/TMEModeling")

RPCData = readMat("data/FitProteinConstraintResults.mat")

sigmas = RPCData$sigmas
ptots = RPCData$ptotAll
proBounds = RPCData$proBounds
growthRates = RPCData$growthAll[1,]



ds = tibble(x=growthRates, y=proBounds)

pS2A = ggplot2::ggplot(ds,ggplot2::aes(x=x,y=y)) +
  ggplot2::geom_point() +
  ggplot2::labs(y="Min. prot. constr.(g/gDW)", x=expression("Growth rate (h"^"-1"*")")) +
  ggplot2::theme_bw() #+ ggplot2::theme(legend.position=legendPos, legend.title = element_blank())

pS2A

mean(proBounds) #0.07859261


ggsave(
  paste0(fig___path, "FigSupProtConstr.png"),
  plot = pS2A,
  width = 6, height = 6, dpi = 300)




