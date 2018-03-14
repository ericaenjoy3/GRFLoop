#' @include GRFLoopClass.R

#' @export infoFilter
setGeneric(name = "infoFilter",
  def = function(loop.obj, fet.obj, info.obj){
    standardGeneric("infoFilter")
  }
)

#' @rdname infoFilter-methods
# updated
setMethod(f = "infoFilter",
  signature = c("loop", "fet", "info"),
  definition = function(loop.obj, fet.obj, info.obj) {
    kpt_idx <- copy(loop.obj@loop[, (gene1 %in% info.obj@gene[["gene"]] | gene1=="") & (gene2 %in% info.obj@gene[["gene"]] | gene2=="")])
    # filter loop slot of loop object by gene1 and gene2 either empty or in gene slot of info object
    loop.obj@loop <- copy(loop.obj@loop[kpt_idx])
    loop.obj@loop[, rowid := 1:nrow(loop.obj@loop)]
    # use filtered loops to subset edges of g slot
    message("remove ", sum(!E(loop.obj@g)$loop %in% unique(loop.obj@loop[["loop"]])), " edges")
    loop.obj@g <- delete.edges(loop.obj@g, which(!E(loop.obj@g)$loop %in% unique(loop.obj@loop[["loop"]])))
    # update vertex of g slot
    loop.obj@g <- delete.vertices(loop.obj@g, which(igraph::degree(loop.obj@g)<1))    
    if (!is.null(loop.obj@split)) {
      loop.obj@split <- factor(loop.obj@split[kpt_idx], levels = unique(loop.obj@split[kpt_idx]))
    }   
    # update fet.obj
    stopifnot(sum(kpt_idx) <= nrow(fet.obj@dat_list[[1]]))
    fet.obj@dat_list <- lapply(fet.obj@dat_list, function(dat)copy(dat[kpt_idx]))
    validObject(fet.obj)
    return(list(loop.obj = loop.obj, fet.obj = fet.obj))
  }
)

# updated
setMethod(f = "infoFilter",
  signature = c("loop", "missing", "info"),
  definition = function(loop.obj, fet.obj, info.obj) {
    kpt_idx <- loop.obj@loop[, (gene1 %in% info.obj@gene[["gene"]] | gene1=="" | is.na(gene1)) & (gene2 %in% info.obj@gene[["gene"]] | gene2=="" | is.na(gene2))]
    # filter loop slot of loop object by gene1 and gene2 either empty or in gene slot of info object
    loop.obj@loop <- loop.obj@loop[kpt_idx]
    loop.obj@loop[, rowid := 1:nrow(loop.obj@loop)]
    # use filtered loops to subset edges of g slot
    message("remove ", sum(!E(loop.obj@g)$loop %in% unique(loop.obj@loop[["loop"]])), " edges")
    loop.obj@g <- delete.edges(loop.obj@g, which(!E(loop.obj@g)$loop %in% unique(loop.obj@loop[["loop"]])))
    # update vertex of g slot
    loop.obj@g <- delete.vertices(loop.obj@g, which(igraph::degree(loop.obj@g)<1))    
    if (!is.null(loop.obj@split)) {
      loop.obj@split <- factor(loop.obj@split[kpt_idx], levels = unique(loop.obj@split[kpt_idx]))
    }
    validObject(loop.obj)
    return(list(loop.obj = loop.obj))
  }
)