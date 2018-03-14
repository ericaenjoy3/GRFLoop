#' @include GRFLoopClass.R GRFLoopGeneric.R

#' @export geneListProc
setGeneric(name = "geneListProc",
  def = function(loop.obj, info.obj, nmin, nmax, type, uniqueLoopGene = FALSE){
    standardGeneric("geneListProc")
  }
)

#' @rdname geneList-methods
setMethod(f = "geneListProc",
  signature = c("loop", "info"),
  definition = function(loop.obj, info.obj, nmin, nmax, type, uniqueLoopGene) {
    # dedup loop in the loop slot of loop.obj
    kpt.idx <- copy(!duplicated(loop.obj@loop[["loop"]]))
    loop_hash <- copy(loop.obj@loop[kpt.idx])
    setkeyv(loop_hash, "loop")
    # identify incident loop
    stopifnot(sum(V(loop.obj@g)$vtype == type) > 0)
    ve <- V(loop.obj@g)$name[V(loop.obj@g)$vtype == type]
    ed <- incident_edges(loop.obj@g, ve)
    idx <- sapply(ed, length) >= nmin & sapply(ed, length) <= nmax
    # extract PromGene from loop slot of loop.obj
    gene_list <- lapply(ed[idx], function(es, loop_hash){
      lp <- as_ids(es)
      gs <- copy(loop_hash[loop %in% lp, gene1])
      stopifnot(all(!is.na(gs)))
      return(gs)
    }, loop_hash = loop_hash)
    if (uniqueLoopGene) {
      gene_list <- gene_list[!duplicated(gene_list)]
    }
    if (type == "Enh") {
    	idx <- which(sapply(gene_list, function(vec)length(unique(vec))) >= nmin)
    	gene_list <- gene_list[idx]
    }
    if (type == "Prom") {
    	gene_list <- lapply(gene_list, function(vec)unique(vec))
    }
    return(gene_list)
  }
)