#' @include GRFClass.R
#' @export rmNonVarRNA
setGeneric(name = "rmNonVarRNA",
  def = function(loop.obj, fet.obj){
    standardGeneric("rmNonVarRNA")
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

#' @export coordShulf
setGeneric(name = "coordShulf",
  def = function(coord, info.obj, dout, nshulf, ncon){
    standardGeneric("coordShulf")
  }
)

#' @export gene2direction
setGeneric(name = "gene2direction",
  def = function(gene_list, info.obj){
    standardGeneric("gene2direction")
  }
)

#' @export gene2direction
setGeneric(name = "shufPlot",
  def = function(loop.obj, info.obj, ncon, dout, pdffout){
    standardGeneric("shufPlot")
  }
)
