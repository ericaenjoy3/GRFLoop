#' @include GRFLoopClass.R

#' @export infoFilter
setGeneric(name = "infoFilter",
  def = function(loop.obj, fet.obj, info.obj){
    standardGeneric("infoFilter")
  }
)

#' @rdname infoFilter-methods
setMethod(f = "infoFilter",
  signature = c("loop", "fet", "info"),
  definition = function(loop.obj, fet.obj, info.obj) {
    idx <- loop.obj@loop[["PromGene"]] %in% info.obj@gene[["gene"]]
    if (all(idx)) {
      return(list(loop.obj = loop.obj, fet.obj = fet.obj))
    }
    message(sum(!idx), " loops (including duplicated) filtered out in infoFilter.")
    # update loop.obj
    stopifnot(length(unique(loop.obj@loop[["loop"]][!idx])) <= length(E(loop.obj@g)))    
    loop.obj@loop <- loop.obj@loop[idx]
    loop.obj@g <- delete.edges(loop.obj@g, which(!gsub("\\|", "_", as_ids(E(loop.obj@g))) %in% unique(loop.obj@loop[["loop"]])))
    loop.obj@g <- delete.vertices(loop.obj@g, which(igraph::degree(loop.obj@g)<1))
    if (is.null(loop.obj@split)) {
      loop.obj@split <- factor(loop.obj@split[idx], levels = unique(loop.obj@split[idx]))
    }
    validObject(loop.obj)    
    # update fet.obj
    stopifnot(length(idx) <= nrow(fet.obj@dat_list[[1]]))
    fet.obj@dat_list <- lapply(fet.obj@dat_list, function(dat)dat[idx])
    validObject(fet.obj)
    return(list(loop.obj = loop.obj, fet.obj = fet.obj))
  }
)
