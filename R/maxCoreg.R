#' @include GRFLoopClass.R GRFLoopGeneric.R

#' @export maxCoreg
setGeneric(name = "maxCoreg",
  def = function(loop.obj, coregfout, uniqueLoopGene = FALSE){
    standardGeneric("maxCoreg")
  }
)

#' @rdname geneList-methods
setMethod(f = "maxCoreg",
  signature = c("loop"),
  definition = function(loop.obj, coregfout, uniqueLoopGene) {
    type <- "Enh"

    # dedup loop in the loop slot of loop.obj
    kpt.idx <- !duplicated(loop.obj@loop[["loop"]])
    loop_hash <- loop.obj@loop[kpt.idx]
    setkeyv(loop_hash, "loop")

    # identify incident loop
    stopifnot(sum(V(loop.obj@g)$vtype == type) > 0)
    ve <- V(loop.obj@g)$name[V(loop.obj@g)$vtype == type]
    ed <- incident_edges(loop.obj@g, ve)
    
    # output genes coregulated by the same enhancer hub
    dd <- rbindlist(lapply(seq_along(ed), function(i) {
      data.table(hub = names(ed)[i],
        partner = gsub(paste0("\\|{0,1}", names(ed)[i], "\\|{0,1}"), "", as_ids(ed[[i]])),
        connum = length(ed[[i]]), 
        genes = loop_hash[chmatch(as_ids(ed[[i]]), loop), gene1])
    }))
    write.table(dd, file = coregfout, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

    # max connection for which the number of hub > 1
    maxco <- data.table(connum = sapply(ed, length))
    maxco <- maxco[, .N, by = connum]

    return(maxco[N>1, max(connum)])
  }
)