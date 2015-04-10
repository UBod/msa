setClass("MsaMetaData", 
        representation(version="character", 
                       params="list", 
                       call="character"))

setClass("MsaDNAMultipleAlignment", 
        contains=c("DNAMultipleAlignment", "MsaMetaData"))

setClass("MsaRNAMultipleAlignment", 
        contains=c("RNAMultipleAlignment", "MsaMetaData"))

setClass("MsaAAMultipleAlignment", 
        contains=c("AAMultipleAlignment", "MsaMetaData"))
