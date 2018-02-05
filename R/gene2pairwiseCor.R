#' @include GRFLoopClass.R

#' @export gene2pairwiseCor
setGeneric(name = "gene2pairwiseCor",
  def = function(gene_list, info.obj){
    standardGeneric("gene2pairwiseCor")
  }
)

#' @rdname gene2pairwiseCor-methods
setMethod(f = "gene2pairwiseCor",
  signature = c("list", "info"),
  definition = function(gene_list, info.obj){
    stopifnot(!is.null(info.obj@gcor))
    cor_vec <- unlist(lapply(gene_list, function(x){
      mat <- combn(x, 2)
      sapply(1:ncol(mat), function(j){
        stopifnot(mat[1,j] %in% rownames(info.obj@gcor))
        return(info.obj@gcor[mat[1,j], mat[2,j]])
      })
    }))
    return(cor_vec)
  }
)