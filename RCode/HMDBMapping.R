################################
# This file contains a lot of manual 
# mapping work between HMDB and our 
# other source for blood concentrations
################################

library(qdapTools)
library(tidyverse)
library(R.matlab)
library(scales)
library(ggplot2)
library(ggpubr)

data____path = "C:/Work/MatlabCode/projects/TMEModeling/TMEModeling/data/"
fig___path = "Z:/projects/TMEModeling/figures/"

metaboliteTable = readRDS(paste0(data____path, "metTable.RDS"))
metaboliteSynonymes = readRDS(paste0(data____path, "metSynonyms.RDS"))
metaboliteVals = readRDS(paste0(data____path, "metVals.RDS"))
metaboliteValRefs = readRDS(paste0(data____path, "metRefs.RDS"))


#make a more advanced calculation
concCalcTable = tibble(met=character(), conc = numeric(), ref=character())
for (i in 1:nrow(metaboliteTable)) {
  if (!is.null(metaboliteVals[[i]])) {
    for (j in 1:length(metaboliteVals[[i]])) {
      #remove spaces in ref, otherwise comparison typically doesn't work
      concCalcTable = concCalcTable %>% add_row(met=metaboliteTable$name[i], conc = metaboliteVals[[i]][j], ref = gsub("[[:space:]]", "", metaboliteValRefs[[i]][j]))
    }
  }
}
dim(concCalcTable)#9023    3
concCalcTableFilt = concCalcTable[!is.na(concCalcTable$ref),]
concCalcTableFilt2 = concCalcTableFilt[concCalcTableFilt$ref != "",]
dim(concCalcTableFilt2)#doesn't filter anything

#first calculate medians within each study
concCalcTableColl = concCalcTableFilt2 %>% group_by(met,ref) %>% summarize(conc=median(conc))
dim(concCalcTableColl)#6423    3
concCalcTableColl
#then medians for each metabolite
concCalcTableColl2 = concCalcTableColl %>% group_by(met) %>% summarize(conc=median(conc), numRefs=n())
dim(concCalcTableColl2)#2753    2
concCalcTableColl2

#identify outlier refs
lookupTab = concCalcTableColl2[, 1:2]
colnames(lookupTab) = c("x", "y")
lookupTab2 = concCalcTableColl2[, c(1,3)]
colnames(lookupTab2) = c("x", "y")
outlTib = concCalcTableColl
outlTib$medianConc = lookup(outlTib$met, lookupTab)
outlTib$numRefs = lookup(outlTib$met, lookupTab2)
outlTibFilt = outlTib[outlTib$numRefs > 1,]
outlTibFilt2 = outlTibFilt[outlTibFilt$medianConc > 50,] #the effect of variation below 50 uM of a single metabolite is not that bad
outlTibFilt2$outlierMetric = abs(log2(outlTibFilt2$conc/outlTibFilt2$medianConc))
outlTibble = outlTibFilt2 %>% group_by(ref) %>% summarize(outlierMetric=mean(outlierMetric), medianConc=median(medianConc)) #medianConc is the same for all in the group

outlTibbleFilt = outlTibble[outlTibble$outlierMetric > 2, ]
#ref                                                                                                                           outlierMetric medianConc
#<chr>                                                                                                                                 <dbl>      <dbl>
#1 AmianoP,DorronsoroM,deRenobalesM,RuizdeGordoaJC,IrigoienI:Very-long-chainomega-3fattyacidsasmarkersforhabitualfishintakeinap~         14.4        89.0
#2 Baselt,RC.Dispositionoftoxicdrugsandchemicalsinman.1982.2ndedition.BiomedicalPublications.Davis,CA.                                    2.60       77.4
#3 ClaytonPT,EatonS,Aynsley-GreenA,EdgintonM,HussainK,KrywawychS,DattaV,MalingreHE,BergerR,vandenBergIE:Hyperinsulinisminshort-~          4.14      150  
#4 DurantonF,CohenG,DeSmetR,RodriguezM,JankowskiJ,VanholderR,ArgilesA:Normalandpathologicconcentrationsofuremictoxins.JAmSocNep~          3.78      171. 
#5 GeigyScientificTables,8thRevedition,pp.123.EditedbyC.Lentner,WestCadwell,N.J.:MedicaleducationDiv.,Ciba-GeigyCorp.,Basel,Swi~          2.37      116. 
#6 HagenfeldtL,vonDobelnU,HolmeE,AlmJ,BrandbergG,EnockssonE,LindebergL:3-Hydroxydicarboxylicaciduria--afattyacidoxidationdefect~          9.84      110. 
#7 HansenJL,FreierEF:Directassaysoflactate,pyruvate,beta-hydroxybutyrate,andacetoacetatewithacentrifugalanalyzer.ClinChem.1978M~          2.06      150  
#8 HjartakerA,LundE,BjerveKS:Serumphospholipidfattyacidcompositionandhabitualintakeofmarinefoodsregisteredbyasemi-quantitativef~          2.54      110  
#9 HuckJH,VerhoevenNM,StruysEA,SalomonsGS,JakobsC,vanderKnaapMS:Ribose-5-phosphateisomerasedeficiency:newinbornerrorinthepentos~          5.72       52.6
#10 KageS,IkedaH,IkedaN,TsujitaA,KudoK:Fatalhydrogensulfidepoisoningatadyeworks.LegMed(Tokyo).2004Jul;6(3):182-6.                          5.09       51.2
#11 KargasG,RudyT,SpennettaT,TakayamaK,QuerishiN,ShragoE:Separationandquantitationoflong-chainfreefattyacidsinhumanserumbyhigh-p~          5.80       89.0
#12 KurikiK,NagayaT,TokudomeY,ImaedaN,FujiwaraN,SatoJ,GotoC,IkedaM,MakiS,TajimaK,TokudomeS:Plasmaconcentrationsof(n-3)highlyunsa~          3.29       89.0
#13 MannN,PirottaY,O'ConnellS,LiD,KellyF,SinclairA:Fattyacidcompositionofhabitualomnivoreandvegetariandiets.Lipids.2006Jul;41(7)~          2.04      110  
#14 MansoorMA,SvardalAM,UelandPM:Determinationoftheinvivoredoxstatusofcysteine,cysteinylglycine,homocysteine,andglutathioneinhum~          2.16       57  
#15 McNaughtonSA,HughesMC,MarksGC:ValidationofaFFQtoestimatetheintakeofPUFAusingplasmaphospholipidfattyacidsandweighedfoodsrecor~          2.53      110  
#16 MinY,GhebremeskelK,LowyC,ThomasB,CrawfordMA:Adverseeffectofobesityonredcellmembranearachidonicanddocosahexaenoicacidsingesta~          2.55      122. 
#17 NakajimaM,YamatoS,WakabayashiH,ShimadaK:High-performanceliquidchromatographicdeterminationofcholesterolandcholestanolinhuman~          3.11     4320  
#18 SakuraN,MizoguchiN,EguchiT,OnoH,MawatariH,NaitouK,ItoK:Elevatedplasmabileacidsinhypergalactosaemicneonates:adiagnosticclueto~          2.38       52.6
#19 ThompsonGN,HsuBY,PittJJ,TreacyE,StanleyCA:Fastinghypoketoticcomainachildwithdeficiencyofmitochondrial3-hydroxy-3-methylgluta~          3.54      150  

#1 should definitely be removed, it is unclear which unit the data is reported in + it diverges a lot

#So, test to just remove them all

concCalcTable #9,023 x 3

concCalcTableFilt3 = concCalcTableFilt2[!(concCalcTableFilt2$ref %in% outlTibbleFilt$ref),]
#also remove this study, the measurements are done during a xylose test:
toRem = gsub("[[:space:]]", "", "Ehrenpreis ED, Salvino M, Craig RM: Improving the serum D-xylose test for the identification of patients with small intestinal malabsorption. J Clin Gastroenterol. 2001 Jul;33(1):36-40.")
concCalcTableFilt3 = concCalcTableFilt3[concCalcTableFilt3$ref != toRem,]
#and these ones, they are pretty old, and the results seem extreme
#It's about melanin
toRem = gsub("[[:space:]]", "", "Hegedus ZL, Frank HA, Steinman TI, Altschule MD, Nayak U: Elevated levels of plasma lipofuscins in patients with chronic renal failure. Arch Int Physiol Biochim. 1988 Dec;96(5):211-21.")
concCalcTableFilt3 = concCalcTableFilt3[concCalcTableFilt3$ref != toRem,]
toRem = gsub("[[:space:]]", "", "Hegedus ZL: The probable involvement of soluble and deposited melanins, their intermediates and the reactive oxygen side-products in human diseases and aging. Toxicology. 2000 Apr 14;145(2-3):85-101.")
concCalcTableFilt3 = concCalcTableFilt3[concCalcTableFilt3$ref != toRem,]
concCalcTableFilt3 #  8,816 x 3, so we remove roughly 200 unreliable values, that is not a large problem, so just remove all outliers

#redo the concentration calculations

#first calculate medians within each study
concCalcTableCollFilt = concCalcTableFilt3 %>% group_by(met,ref) %>% summarize(conc=median(conc))
dim(concCalcTableCollFilt)#6232    3
concCalcTableCollFilt
#then medians for each metabolite
concCalcTableCollFilt2 = concCalcTableCollFilt %>% group_by(met) %>% summarize(conc=median(conc), numRefs=n())
dim(concCalcTableCollFilt2)#2715    3 #we lost a few metabolites
concCalcTableCollFilt2
#look at the new table
srtNew = sort(concCalcTableCollFilt2$conc, decreasing=TRUE, index.return=TRUE)
print(concCalcTableCollFilt2[srtNew$ix,], n=100)
concCalcTableCollFilt2[srtNew$ix,]$conc


#The levels of Fructosamine, ATP and Glyceraldehyde seems very high


#
#tst = concCalcTable$ref[grepl( "10.1016/j.jaci.2015.01.022", concCalcTable$ref , fixed = TRUE)]
#remove white space
#tst2 = gsub("[[:space:]]", "", tst)
#tst
#unique(tst)
#unique(tst2)#works

#concCalcTableColl$ref[grepl( "10.1016/j.jaci.2015.01.022", concCalcTableColl$ref , fixed = TRUE)]

"10.1016/j.jaci.2015.01.022"



#first filter out NA concentrations
#in total 25,411 metabolites
#ntFilt1 = metaboliteTable[!is.na(metaboliteTable$conc),] #2,743 left
dsTmp = concCalcTableCollFilt2
colnames(dsTmp)[1] = "name"
dsTmp2 = metaboliteTable[,c(1,2,4)]
dsTmp2$index = 1:nrow(metaboliteTable)

ntFilt1 = inner_join(dsTmp, dsTmp2, by="name")

#get rid of metabolites with negligible concentrations.
ntFilt1 = ntFilt1[ntFilt1$conc >= 1,]
#get rid of all lipids and lipid-like molecules, we are not interested in them
ntFilt1 = ntFilt1[ntFilt1$superclass != "Lipids and lipid-like molecules",]


ntFilt1#307 metabolites

synonyms = metaboliteSynonymes[ntFilt1$index]


#read the matlab file with all metabolite names
setwd("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Johan/OptimalTMEGrowthStrategy")

modelMetsX = readMat("data/metsImportedByModel.mat")
modelMetsY = unlist(modelMetsX$metsImportedByModel)
modelMets = as.character(modelMetsY)


metNames = ntFilt1$name
metNamesLower = tolower(metNames)

matchingMetNames = rep("", length(modelMets))
sum(tolower(modelMets) %in% metNamesLower) #64
lowerModelMets = tolower(modelMets);
#first fill in directly matching names for which measurements exist
#sel = tolower(modelMets) %in% metNamesLower
for (i in 1:length(modelMets)) {
  sel = metNamesLower == lowerModelMets[i];
  if (sum(sel) == 1) {
    matchingMetNames[i] = metNames[sel]
  } else if (sum(sel) > 1) {
    print(paste("Multiple matches:", modelMets[i]))
  }
}
#test
print(tibble(modelMets, matchingMetNames),n=300) #, looks ok


#then fill in those that match on synonyms matching names for which measurements exist
synonymsFilt = synonyms
for (i in 1:length(synonymsFilt)) {
  synonymsFilt[[i]] = tolower(synonymsFilt[[i]])
}

multipleMatches = list();

for (i in 1:length(modelMets)) {
  if (matchingMetNames[i] == "") { #don't override direct matches
    sel = rep(FALSE, length(synonymsFilt))
    for(j in 1:length(synonymsFilt)) {
      sel[j] = sum(lowerModelMets[i] == synonymsFilt[[j]]) > 0
    }
    if (sum(sel) == 1) {
      matchingMetNames[i] = metNames[sel]
    } else if (sum(sel) > 1) {
      multipleMatches = append(multipleMatches, list(c(modelMets[i], metNames[sel])))
      print(paste("Multiple matches:", modelMets[i], ":", paste(metNames[sel])))
    }
  }
}

#setMMN = function(mmn, modelMets, modelName, matchedName) {
#  mmn[modelMets == modelName] = matchedName
#  return(mmn)
#}

#manually fix the multiple matches:
sum(matchingMetNames != "") #178
length(multipleMatches) #9
#ideal length after fix
178 + 9#187

matchingMetNames[modelMets == "acetate"] = "Acetic acid"
matchingMetNames[modelMets == "alanine"] = "L-Alanine" #this is the type incorporated into proteins
#matchingMetNames[modelMets == "CO"] = "Carbon monoxide"
#matchingMetNames[modelMets == "cystine"] = "L-Cystine"
#matchingMetNames[modelMets == "DHA"] = "Docosahexaenoic acid"
matchingMetNames[modelMets == "D-lactate"] = "D-Lactic acid"
#metaboliteTable[metaboliteTable$name == "Docosapentaenoic acid (22n-6)",]#10
#metaboliteTable[metaboliteTable$name == "Docosapentaenoic acid (22n-3)",]#30.7
#DPA can mean both of these, for simplicity we use the most highly abundant
#matchingMetNames[modelMets == "DPA"] = "Docosapentaenoic acid (22n-3)"
#matchingMetNames[modelMets == "eicosanoate"] = "Arachidic acid"
#metaboliteTable[metaboliteTable$name == "L-Glutamic acid",]#65.4
#metaboliteTable[metaboliteTable$name == "D-Glutamic acid",]#0.107
#we use the most abundant form for simplicity
#matchingMetNames[modelMets == "glutamate"] = "L-Glutamic acid"
matchingMetNames[modelMets == "glycerate"] = "Glyceric acid"
matchingMetNames[modelMets == "glycocholate"] = "" # a correction, this one is just mapped wrongly
#matchingMetNames[modelMets == "hexanoic acid"] = "Caproic acid"
matchingMetNames[modelMets == "isoleucine"] = "L-Isoleucine"
#metaboliteTable[metaboliteTable$name == "alpha-Linolenic acid",]#28.7
#metaboliteTable[metaboliteTable$name == "gamma-Linolenic acid",]#12.1
#we use the most abundant form for simplicity
#matchingMetNames[modelMets == "linolenate"] = "gamma-Linolenic acid"
matchingMetNames[modelMets == "L-lactate"] = "L-Lactic acid"
matchingMetNames[modelMets == "nicotinate"] = "Nicotinic acid"
#metaboliteTable[metaboliteTable$name == "5-Tetradecenoic acid",]#1.5
#metaboliteTable[metaboliteTable$name == "5Z-Tetradecenoic acid",]#1.75
#they are basically the same, use the one without Z
#matchingMetNames[modelMets == "physeteric acid"] = "5-Tetradecenoic acid"
#metaboliteTable[metaboliteTable$name == "Prostaglandin F2a",]#0.000251
#metaboliteTable[metaboliteTable$name == "11b-PGF2a",]#0.000181
#use the first
#matchingMetNames[modelMets == "prostaglandin F2alpha"] = "Prostaglandin F2a"
#metaboliteTable[metaboliteTable$name == "all-trans-Retinoic acid",]#0.107
#metaboliteTable[metaboliteTable$name == "9-cis-Retinoic acid",]#0.0049
#metaboliteTable[metaboliteTable$name == "13-cis-Retinoic acid",]#0.003
#all low, use the fist
#matchingMetNames[modelMets == "retinoate"] = "all-trans-Retinoic acid"
#matchingMetNames[modelMets == "benzoate"] = "Benzoic acid"
matchingMetNames[modelMets == "glucuronate"] = "D-Glucuronic acid"
matchingMetNames[modelMets == "fumarate"] = "Fumaric acid"
#metaboliteTable[metaboliteTable$name == "14R,15S-EpETrE",]#0.000134
#metaboliteTable[metaboliteTable$name == "14,15-Epoxy-5,8,11-eicosatrienoic acid",]#0.000942
#all low, use the second
#matchingMetNames[modelMets == "14,15-EET"] = "14,15-Epoxy-5,8,11-eicosatrienoic acid"
matchingMetNames[modelMets == "homogentisate"] = "" # "Homogentisic acid" is correct, but filtered, glycerol is just wrong, so clear it
#metaboliteTable[metaboliteTable$name == "17-Hydroxyprogesterone",]#0.00808
#metaboliteTable[metaboliteTable$name == "3-(3-Hydroxyphenyl)propanoic acid",]#0.144
#both low, use the first
#matchingMetNames[modelMets == "17alpha-hydroxyprogesterone"] = "17-Hydroxyprogesterone"
#matchingMetNames[modelMets == "13-cis-retinoate"] = "13-cis-Retinoic acid"
#metaboliteTable[metaboliteTable$name == "Ganglioside GM3 (d18:1/16:0)",]# 4.5
#metaboliteTable[metaboliteTable$name == "Ganglioside GM3 (d18:1/26:1(17Z)))",]#0.567
#both low, use the first
#matchingMetNames[modelMets == "GM3"] = "Ganglioside GM3 (d18:1/16:0)"

sum(matchingMetNames != "") #185, ok, we removed two that had the wrong mapping

sum(is.na(matchingMetNames))
length(unique(matchingMetNames[matchingMetNames != ""])) #181, so there are some model metabolites that are mapped to the same - fix that
matchedNames = matchingMetNames[matchingMetNames != ""]
dupl = matchedNames[duplicated(matchedNames)]
dupl
#[1] "L-Lysine"        "L-Phenylalanine" "L-Kynurenine"    "L-Arabitol"  

#[1] "gamma-Linolenic acid"     "L-Lysine"                 "L-Phenylalanine"          "Phytanic acid"            "Vitamin A"               
#[6] "Chenodeoxycholic acid"    "L-Kynurenine"             "L-Arabitol"               "3-Hydroxyisovaleric acid"

#[1] "Adenosine"                 "gamma-Linolenic acid"      "L-Lysine"                  "Testosterone"              "L-Phenylalanine"          
#[6] "Phytanic acid"             "Vitamin A"                 "Chenodeoxycholic acid"     "L-Kynurenine"              "Norepinephrine"           
#[11] "20-Hydroxy-leukotriene B4" "Thromboxane B2"            "3-Hydroxyisovaleric acid"  "Dihydrobiopterin"          "15-HEPE"  

#investigate them 1 by 1:
for (i in 1:length(dupl)) {
  print(paste(dupl[i], paste(modelMets[matchingMetNames == dupl[i]])))
} 

# [1] "L-Lysine L-lysine" "L-Lysine lysine"  
# [1] "L-Phenylalanine phenylalanine" "L-Phenylalanine THF"          
# [1] "L-Kynurenine kynurenine"             "L-Kynurenine 3-hydroxy-L-kynurenine"
# [1] "L-Arabitol L-arabitol" "L-Arabitol D-Arabitol"
#matchingMetNames[modelMets == "deoxyadenosine"] = "" #wrongly matched to adenosine
#matchingMetNames[modelMets == "linolenate"] = "" 
matchingMetNames[modelMets == "L-lysine"] = ""#L-lysine doesn't really exist, this is a bug in Gecko, it should be [protein]-L-lysine
#matchingMetNames[modelMets == "testosterone sulfate"] = ""
matchingMetNames[modelMets == "THF"] = "" 
#matchingMetNames[modelMets == "Phytanate"] = "" #So, these two are the same thing more or less, I'm keeping the one with most references in the model
#matchingMetNames[modelMets == "11-cis-retinol"] = ""
#matchingMetNames[modelMets == "chenodiol"] = "" #So, these two are the same thing more or less, I'm keeping the one with most references in the model
matchingMetNames[modelMets == "3-hydroxy-L-kynurenine"] = ""
#matchingMetNames[modelMets == "2-oxoadipate"] = ""
#matchingMetNames[modelMets == "5,12,20-TriHETE"] = "" #So, these two are the same thing more or less, I'm keeping the one with most references in the model
#matchingMetNames[modelMets == "Thromboxane B2"] = "" #So, these two are obviously the same thing, I'm keeping the one with most references in the model
matchingMetNames[modelMets == "D-Arabitol"] = "" 
#matchingMetNames[modelMets == "beta-hydroxy-beta-methylbutyrate"] = "" 
#matchingMetNames[modelMets == "quinonoid dihydrobiopterin"] = ""
#matchingMetNames[modelMets == "15(S)-HEPE"] = "" #15(R)-HEPE is more commonly used in the model

sum(matchingMetNames != "") #181, ok

#We have now done the following
# * Removed the ones with concentration below 1uM to save some work, they really don't matter
# * Removed all lipids, we already have a way to handle those. 
# We now want to filter out the rest of the metabolites that we have already handled, which are in the BloodData.txt file

metaboliteCandidates = matchingMetNames
modelMetaboliteCandidates = modelMets
metaboliteCandidatesFilt = metaboliteCandidates[metaboliteCandidates != ""]
modelMetaboliteCandidatesFilt = modelMetaboliteCandidates[metaboliteCandidates != ""]
print(tibble(metaboliteCandidatesFilt,modelMetaboliteCandidatesFilt), n=300)

#get the metabolites currently included in the simulation file
currMetTable = read_tsv(paste0(data____path, "BloodData.txt"))
currMetTable = currMetTable[,1:4]
sel1 = modelMetaboliteCandidatesFilt %in% currMetTable$Metabolite
sum(sel1) #56, so we have found most of them (not the free fatty acids and a few more). Remove them from the list
metaboliteCandidatesFilt2 = metaboliteCandidatesFilt[!sel1]
modelMetaboliteCandidatesFilt2 = modelMetaboliteCandidatesFilt[!sel1]
#also check which ones used that are not found in HMDB
currMetTable[!(currMetTable$Metabolite %in% modelMetaboliteCandidatesFilt),]

print(tibble(metaboliteCandidatesFilt2,modelMetaboliteCandidatesFilt2), n=300)
existingInBoth = inner_join(ntFilt1, tibble(name=metaboliteCandidatesFilt[sel1], modelName = modelMetaboliteCandidatesFilt[sel1]), by="name")
existingInBothMerged = inner_join(existingInBoth, tibble(modelName=currMetTable[[1]], concFromBloodData = currMetTable[[2]]), by="modelName")
plot(log2(existingInBothMerged$concFromBloodData), log2(existingInBothMerged$conc))
print(existingInBothMerged, n=100)

#create supplementary figure
labels = rep("", length(existingInBothMerged$concFromBloodData))
labels[existingInBothMerged$name == "D-Glucose"] = "Glucose"
labels[existingInBothMerged$name == "L-Lactic acid"] = "Lactate"
labels[existingInBothMerged$name == "D-Glucuronic acid"] = "Glucuronic acid"
labels[existingInBothMerged$name == "Oxygen"] = "Oxygen"

ds = tibble(x=existingInBothMerged$concFromBloodData, y=existingInBothMerged$conc, label=labels)
#replace the glucose value with that from Akinci et al:
ds$x[existingInBothMerged$name == "D-Glucose"] = 4860
dsLine = tibble(x=c(1,5500), y=c(1,5500))

pA = ggplot(ds,aes(x=x,y=y)) +
  geom_point() +
  geom_text(mapping = aes(label=labels), hjust=0.5, vjust=-1) +
  ggplot2::geom_line(data=dsLine) +
  labs(y=expression("HMDB blood conc ( "*mu*"M )"), x=expression("Used plasma conc ( "*mu*"M )")) +
  theme_bw() + 
  theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_cartesian(xlim = c(1,11000), ylim = c(1,11000)) +
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2')

pA

ggsave(
  paste0(fig___path, "FigSupBloodConc.png"),
  plot = pA,
  width = 6, height = 6, dpi = 300)



#investigate the candidates to add to the BloodData.txt file
cands = inner_join(ntFilt1, tibble(name=metaboliteCandidatesFilt2, modelName = modelMetaboliteCandidatesFilt2), by="name")
print(cands, n=200)
#Remove some metabolites that seem way off:
#1. #ATP, ADP and AMP concentrations are way off. Investigation of one source, Physiological concentrations of purines and pyrimidines (Traut)
#gave that it seems that this is measured mostly within cells, or possibly in whole blood. However, Human plasma ATP concentration (Gorman)
#says ATP levels is roughly 1 uM, compared to 2271 in HMDB. So, this is likely negligible, just remove them.
cands  = cands[!(cands$name %in% c("Adenosine triphosphate", "ADP", "Adenosine monophosphate", "Guanosine diphosphate", "Guanosine triphosphate", "Uridine 5'-diphosphate", "Uridine 5'-monophosphate")),]
#D-alanine seems to be there, in about the same concentrations as L-alanine. The problem is that D-alanine is a rare thing that is not really
#there, so this must be a mistake. Remove it.
cands  = cands[cands$name != "D-Alanine",]
#Glyceraldehyde levels are very high, and it is based on a single value. There are no other reports of this compound, so it is likely not very abundant. 
#Lactate levels in the same study are also very high 3k+ uM. We remove it because we don't have good data, although it may be a compound that should be included.
cands  = cands[cands$name != "Glyceraldehyde",]
srt = sort(cands$conc, decreasing=TRUE, index.return=TRUE )
sortedCands = cands[srt$ix, ]
print(sortedCands, n=600)
nrow(sortedCands)
#Now, we finally want to check that we have not missed any important metabolites (i.e. any highly abundant 
#metabolites where we failed to match the metabolite names to the model)
write_tsv(sortedCands, paste0(data____path, "extraMetCandidatesFromHMDB.txt"))

#matchedMets = matchingMetNames[matchingMetNames != ""]

#Now check that we didn't miss any metabolites with high concentration in blood
#metTableMatched = ntFilt1[ntFilt1$name %in% matchedMets, ]
#srtMatched = sort(metTableMatched$conc, decreasing=TRUE, index.return=TRUE)
#print(metTableMatched[srtMatched$ix,], n=600)
#metTableUnmatched = ntFilt1[!(ntFilt1$name %in% matchedMets), ]
#print(metTableUnmatched[srtMatched$ix,], n=600)

#first filter out the ones that we have managed to map to the model metabolites
unmCands = ntFilt1[!(ntFilt1$name %in% matchingMetNames[matchingMetNames != ""]), ]
srt = sort(unmCands$conc, decreasing=TRUE, index.return=TRUE )
sortedUnmCands = unmCands[srt$ix, ]
print(sortedUnmCands, n=600)
write_tsv(sortedUnmCands, paste0(data____path, "unmappedMetsFromHMDB.txt"))
