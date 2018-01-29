#' @include GRFLoopClass.R

#' @export inTADShulf
setGeneric(name = "inTADShulf",
  def = function(gene_list, info.obj){
    standardGeneric("inTADShulf")
  }
)

#' @rdname inTADShulf-methods
setMethod(f = "inTADShulf",
  signature = c("list", "info"),
  definition = function(gene_list, info.obj){
    genep_list <- lapply(seq_along(gene_list), function(j){
      gid <- gene_list[[j]]
      tads <- unique(info.obj@gene[gene %in% gid, tadid])
      set.seed(j)
      rand_gid <- info.obj@gene[tadid %in% tads & !gene %in% gid, sample(gene, size = length(gid), replace = FALSE)]
      return(rand_gid)
    })
    return(genep_list)
  }
)
