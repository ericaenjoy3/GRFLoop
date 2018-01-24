#' @include GRFClass.R
#' @export rmNonVarRNA
setGeneric(name = "rmNonVarRNA",
  def = function(loop.obj, fet.obj){
    standardGeneric("rmNonVarRNA")
  }
)

#' @export infoFilter
setGeneric(name = "infoFilter",
  def = function(loop.obj, fet.obj, info.obj){
    standardGeneric("infoFilter")
  }
)

#' @export ProteinCodingInfo
setGeneric(name = "ProteinCodingInfo",
  def = function(info.obj){
    standardGeneric("ProteinCodingInfo")
  }
)

#' @export TPMInfo
setGeneric(name = "TPMInfo",
  def = function(info.obj){
    standardGeneric("TPMInfo")
  }
)

#' @export quantRm
setGeneric(name = "quantRm",
  def = function(loop.obj, fet.obj, dedup){
    standardGeneric("quantRm")
  }
)

#' @export orderLoop
setGeneric(name = "orderLoop",
  def = function(loop.obj, fet.obj, sm_nm, order_method = c("last_column", "diff", "last_to_first", "quant")){
    standardGeneric("orderLoop")
  }
)

#' @export ceilFet
setGeneric(name = "ceilFet",
  def = function(fet.obj){
    standardGeneric("ceilFet")
  }
)

#' @export multiHeat
setGeneric(name = "multiHeat",
  def = function(loop.obj, fet.obj, rng_list, pdffout, scheme = c("wr", "br")){
    standardGeneric("multiHeat")
  }
)

#' @export medianViolin
setGeneric(name = "medianViolin",
  def = function(loop.obj, fet.obj, pdffout){
    standardGeneric("medianViolin")
  }
)

#' @export conVBoxPlot
setGeneric(name = "conPlot",
  def = function(loop.obj, fet.obj, dout){
    standardGeneric("conPlot")
  }
)

#' @export gene2direction
setGeneric(name = "gene2direction",
  def = function(gene_list, info.obj){
    standardGeneric("gene2direction")
  }
)

#' @export coordShulf
setGeneric(name = "coordShulf",
  def = function(coord, info.obj, dout, nshulf, nmin, nmax){
    standardGeneric("coordShulf")
  }
)

#' @export inTADShulf
setGeneric(name = "inTADShulf",
  def = function(gene_list, info.obj){
    standardGeneric("inTADShulf")
  }
)

#' @export TADshipPlot
setGeneric(name = "TADshipPlot",
  def = function(gene_list_list, info.obj, nms, pdffout){
    standardGeneric("TADshipPlot")
  }
)

#' @export shufPlot
setGeneric(name = "shufPlot",
  def = function(loop.obj, info.obj, nmin, nmax, dout, tadStatpdf, coregBoxpdf){
    standardGeneric("shufPlot")
  }
)
