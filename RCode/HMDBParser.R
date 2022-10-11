#xml SAX parser for HMDB
####################

#install.packages("XML")
library(XML)
library(tidyverse)

#fileName <- system.file("exampleData", "mtcars.xml", package="XML")

metaboliteTable = tibble(name=character(), MW = numeric(), conc = numeric(), superclass = character())
metaboliteSynonymes = list()
metaboliteVals = list()
metaboliteValRefs = list()


currName = ""
currMW = NA
currConcVec = NULL
currConcRefVec = NULL
currConc = NA
currConcUnit = NA
currBioSpec = ""
currSynonymList = NULL
currConcRef = NA


onMet = FALSE
onName = FALSE
onMW = FALSE
inMet = FALSE
inNC = FALSE
inConc = FALSE
onSuperClass = FALSE;
onBioSpec = FALSE
onConcVal = FALSE
onConcUnit = FALSE
onSynonymList = FALSE
onSynonym = FALSE
onConcRef = FALSE
depth = 0
startElem = function(name, attrs) {
  depth <<- depth + 1
  if (name == "metabolite") {
    inMet <<- TRUE;
    currName <<- ""
    currMW <<- NA
    currConcVec <<- NULL
    currConcRefVec <<- NULL
    currAltNames <<- list();
    currSynonymList <<- NULL;
    currSuperClass <<- "";
  } else if (name == "name") {
    onName <<- TRUE; 
  } else if (name == "synonyms") {
    onSynonymList <<- TRUE; 
  } else if (name == "synonym" && onSynonymList) {
    onSynonym <<- TRUE; 
  } else if (name == "average_molecular_weight") {
    onMW <<- TRUE
  } else if (name == "normal_concentrations") {
    inNC <<- TRUE
  } else if (name == "concentration" && inNC) {
    currBioSpec <<- ""
    currConc <<- NA
    currConcUnit <<- NA;
    currConcRef <<- NA;
    inConc <<- TRUE
    onBioSpec <<- FALSE
    onConcVal <<- FALSE
  } else if (name == "biospecimen" && inConc) {
    onBioSpec <<- TRUE
  } else if (name == "concentration_value" && inConc) {
    onConcVal <<- TRUE
  } else if (name == "concentration_units" && inConc) {
    onConcUnit <<- TRUE
  } else if (name == "reference_text" && inConc) {
    onConcRef <<- TRUE
  } else if (name == "super_class") {
    onSuperClass <<- TRUE
  } else {
    inMet <<- FALSE
  }
}

txtElem = function(txt) {
  if (onName && depth == 3) { 
    currName <<- txt
  }
  if (onSynonym && depth == 4) { 
    currSynonymList <<- c(currSynonymList, txt)
  }
  else if (onMW) { 
    currMW <<- as.numeric(txt)
  }
  else if (onSuperClass) { 
    currSuperClass <<- txt
  }
  else if (onBioSpec) { 
    currBioSpec <<- txt
  }
  else if (onConcUnit) {
    currConcUnit <<- txt;
  }
  else if (onConcRef) {
    if(is.na(currConcRef)) {
      currConcRef <<- txt;
    } else {
      currConcRef <<- paste0(currConcRef, txt)
    }
  }
  else if (onConcVal) { 
    tstval = suppressWarnings(as.numeric(txt))
    success = FALSE
    if (!is.na(tstval)) {
      currConc <<- tstval
      success = TRUE
    } 
    #handle values of the type 9.1 +/- 7.2
    if (!success) {
      ind1 = str_locate(txt, "\\+/-")[1]
      if (!is.na(ind1)) {
        txtToProc = substr(txt, 1, ind1-1)
        val = as.numeric(txtToProc)
        if (!is.na(val)) {
          currConc <<- val
          success = TRUE
        }
      } 
    }
    if (!success) {
      #handle the format "137-148"
      ind1 = str_locate(txt, "-")[1]
      if (!is.na(ind1)) {
        n1 = suppressWarnings(as.numeric(substr(txt, 1, ind1-1)))
        n2 = suppressWarnings(as.numeric(substr(txt, ind1+1, str_length(txt))))
        mn = mean(c(n1,n2))
        if (!is.na(mn)) {
          currConc <<- mn
          success = TRUE
        } 
      }         
    }
    if (!success) {
      #47.06(17.31-146.8)
      ind1 = str_locate(txt, "\\(")[1]
      if (!is.na(ind1)) {
        txtToProc = substr(txt, 1, ind1-1)
        val = as.numeric(txtToProc)
        if (!is.na(val)) {
          currConc <<- val
          success = TRUE
        }
      } 
    }
  }
}

endElem = function(name, attrs) {
  depth <<- depth - 1
  if (name == "metabolite") {
    if (is.null(currConcVec)) cc = NA else cc = mean(currConcVec)
    metaboliteTable <<- metaboliteTable %>% add_row(name=currName, MW = currMW, conc = cc, superclass = currSuperClass)
    metaboliteSynonymes <<- append(metaboliteSynonymes, list(currSynonymList))
    metaboliteVals <<- append(metaboliteVals, list(currConcVec))
    metaboliteValRefs <<- append(metaboliteValRefs, list(currConcRefVec))
  } 
  else if (name == "name") onName <<- FALSE
  else if (name == "average_molecular_weight") onMW <<- FALSE
  else if (name == "normal_concentrations") inNC <<- FALSE
  else if (name == "concentration" && inConc) { 
    inConc <<- FALSE 
    if (currBioSpec == "Blood" && !is.na(currConc)) {
      currConcRefVec <<- c(currConcRefVec, currConcRef)
      if (currConcUnit == "uM") {
        currConcVec <<- c(currConcVec, currConc)
      } 
      else if (currConcUnit == "mM") {
        currConcVec <<- c(currConcVec, currConc*1000)
      }
    }
  }
  else if (name == "biospecimen") onBioSpec <<- FALSE
  else if (name == "concentration_value") onConcVal <<- FALSE
  else if (name == "concentration_units") onConcUnit <<- FALSE
  else if (name == "synonyms") onSynonymList <<- FALSE
  else if (name == "synonym" ) onSynonym <<- FALSE
  else if (name == "reference_text" && inConc) onConcRef <<- FALSE
  else if (name == "super_class") onSuperClass <<- FALSE
}


#test:
#xmlEventParse("F:/BloodData/test.xml",
#              list(startElement=startElem,
#                   endElement = endElem,
#                   text = txtElem),
#              useTagName=FALSE, addContext = FALSE)



#metaboliteTable
#currConcVec
#metaboliteSynonymes
#metaboliteVals
#metaboliteValRefs

#7.7, 14.4, 1.13-166.96, 47.06
#mean(c(7.7, 14.4, mean(c(1.13,166.96)), 47.06)) #38.30125

filename = "C:/Work/MatlabCode/projects/TMEModeling/TMEModeling/data/serum_metabolites.xml"


#The real run - takes a long time to run (15 minutes or so)
xmlEventParse(filename,
              list(startElement=startElem,
                   endElement = endElem,
                   text = txtElem),
              useTagName=FALSE, addContext = FALSE)
warnings() #16 warnings about NANs, that is ok


data____path = "C:/Work/MatlabCode/projects/TMEModeling/TMEModeling/data/"
saveRDS(metaboliteTable, paste0(data____path, "metTable.RDS"))
saveRDS(metaboliteSynonymes, paste0(data____path, "metSynonyms.RDS"))
saveRDS(metaboliteVals, paste0(data____path, "metVals.RDS"))
saveRDS(metaboliteValRefs, paste0(data____path, "metRefs.RDS"))

