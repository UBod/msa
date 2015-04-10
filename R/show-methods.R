setMethod("show", signature("MsaMetaData"),
          function(object) print(object, show=c("version", "standardParams",
                                                "algParams", "call")))

setMethod("show", signature("MsaDNAMultipleAlignment"),
          function(object) print(object))
setMethod("show", signature("MsaRNAMultipleAlignment"),
          function(object) print(object))
setMethod("show", signature("MsaAAMultipleAlignment"),
          function(object) print(object))
