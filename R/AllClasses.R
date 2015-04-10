setClass("MsaMetaData",
         slots=c(version="character",
                 params="list",
                 call="character"),
         contains="VIRTUAL")

setClass("MsaDNAMultipleAlignment",
         contains=c("DNAMultipleAlignment", "MsaMetaData"))

setClass("MsaRNAMultipleAlignment",
         contains=c("RNAMultipleAlignment", "MsaMetaData"))

setClass("MsaAAMultipleAlignment",
         contains=c("AAMultipleAlignment", "MsaMetaData"))
