#' @include GRFClass.R

#' @export TPMInfo
setGeneric(name = "TPMInfo",
  def = function(info.obj){
    standardGeneric("TPMInfo")
  }
)

#' @rdname TPMInfo-methods
setMethod(f = "TPMInfo",
  signature = c("info"),
  definition = function(info.obj) {
    # genes with TPM >= 1 at any stage during reprogramming
    col_nm <- colnames(info.obj@gene)[grep("^tpm_(MEF|ESC)", colnames(info.obj@gene))]
    idx <- info.obj@gene[, apply(.SD, 1, function(vec)any(vec>1)), .SDcols = col_nm]
    stopifnot(all(!idx) == FALSE)
    # update info.obj
    message(sum(!idx), " genes that did not reach TPM > 1 at any stage of reprogramming were filtered.")
    info.obj@gene <- info.obj@gene[idx]
    validObject(info.obj)
    return(info.obj = info.obj)
  }
)
