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
    message("beginning inTADShulf")
    genep_list <- lapply(seq_along(gene_list), function(j) {

      gid <- unique(unlist(gene_list[[j]]))
      tads <- copy(info.obj@gene[gene %in% gid, unique(tadid)])

      set.seed(j)
      rand_gid <- if (nrow(info.obj@gene[tadid %in% tads & !gene %in% gid]) < length(gid)) {
        if (nrow(info.obj@gene[tadid %in% tads]) < length(gid)) {
          copy(info.obj@gene[tadid %in% tads, gene])
        } else {
          copy(info.obj@gene[tadid %in% tads, sample(gene, size = length(gid), replace = FALSE)])
        }
      } else {
        copy(info.obj@gene[tadid %in% tads & !gene %in% gid, sample(gene, size = length(gid), replace = FALSE)])
      }  
      return(rand_gid)
    })
    message("completing inTADShulf")
    return(genep_list)
  }
)
