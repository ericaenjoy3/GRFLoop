#' @include GRFLoopClass.R GRFLoopGeneric.R

#' @export geneCor
setGeneric(name = "geneCor",
  def = function(info.obj){
    standardGeneric("geneCor")
  }
)

#' @rdname geneCor-methods
setMethod(f = "geneCor",
  signature = c("info"),
  definition = function(info.obj) {
  	col_nm <- colnames(info.obj@gene)[grep("tpm_", colnames(info.obj@gene))]
  	stopifnot(length(col_nm) > 1) 
  	stopifnot(all(info.obj@gene[, apply(.SD, 1, var), .SDcols = col_nm] > 0)) 
  	gcor_dat <- dcast(melt(info.obj@gene[, c("gene", col_nm), with = FALSE], id.vars = "gene"), variable ~ gene)[, !"variable"] %>% cor(method = "spearman")
  }
)