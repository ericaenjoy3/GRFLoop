#' @include GRFLoopClass.R

#' @export ProteinCodingInfo
setGeneric(name = "ProteinCodingInfo",
  def = function(info.obj){
    standardGeneric("ProteinCodingInfo")
  }
)

#' @rdname ProteinCodingInfo-methods
setMethod(f = "ProteinCodingInfo",
  signature = c("info"),
  definition = function(info.obj) {
    # genes of protein coding
    gid <- copy(info.obj@gene[type == "protein_coding", gene])
    # update info.obj
    message(sum(!info.obj@gene[, gene] %in% gid), " non-protein-coding genes were filtered.")
    info.obj@gene <- copy(info.obj@gene[gene %in% gid])
    validObject(info.obj)
    return(info.obj = info.obj)
  }
)
