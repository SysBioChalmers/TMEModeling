#code for finding the amino acid composition in collagen
############################################################

library(tidyverse)

setwd("C:/Work/MatlabCode/projects/TMEModeling/TMEModeling")

getFracs = function(filename, AAConvTable) {
  fl = file(paste0("data/ECMBiomass/", filename))
  fasta = readLines(fl)
  close(fl)
  #skip the first line, it is some kind of header
  sequence = strsplit(paste0(fasta[-1], collapse=''), "")[[1]]
  #letters = unique(strsplit(sequence, "")[[1]])

  results = cbind(AAConvTable, tibble(fraction=rep(0,20)))
  
  for (i in 1:20) {
    results$fraction[i] = sum(sequence == results$letters[i])
  }
  results$fraction = results$fraction / sum(results$fraction)
  
  return (results)
}



AAFullNames = c("Alanine", "Arginine", "Asparagine", "Aspartic acid", "Cysteine", "Glutamine", "Glutamic acid", "Glycine", "Histidine", "Isoleucine",
              "Leucine", "Lysine", "Methionine", "Phenylalanine", "Proline", "Serine", "Threonine", "Tryptophan", "Tyrosine", "Valine")
AAShortNames = c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val")
AALetters = c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

AAConvTable = tibble(fullNames = AAFullNames, shortNames = AAShortNames, letters = AALetters)

col1ChainAlpha1 = getFracs("P02452.fasta.txt", AAConvTable)
col1ChainAlpha2 = getFracs("P08123.fasta.txt", AAConvTable)

#test
sum(col1ChainAlpha1$fraction) #1, ok
sum(col1ChainAlpha2$fraction) #1, ok

totCol1 = col1ChainAlpha1
totCol1$fraction = (col1ChainAlpha1$fraction*2 + col1ChainAlpha2$fraction)/3
write_tsv(totCol1, "data/ECMBiomass/Col1AAFractions.txt")


























