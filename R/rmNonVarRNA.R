#' @include GRFLoopClass.R

#' @export rmNonVarRNA
setGeneric(name = "rmNonVarRNA",
  def = function(loop.obj, fet.obj){
    standardGeneric("rmNonVarRNA")
  }
)

#' @rdname rmNonVarRNA-methods
setMethod(f = "rmNonVarRNA",
  signature = c("loop", "fet"),
  definition = function(loop.obj, fet.obj) {
    idx <- grep("RNA", fet.obj@hash[["sms"]], ignore.case = TRUE)[1]
    kpt.idx <- which(apply(fet.obj@dat_list[[idx]], 1, var) > 0)
    stopifnot(length(kpt.idx) <= nrow(fet.obj@dat_list[[idx]]))
    # update dat_list slot of fet class
    # update loop slot of loop class
    # update g slot of loop class
    # update split slot of loop class
    if (length(kpt.idx) > 0 & length(kpt.idx) < nrow(fet.obj@dat_list[[idx]])) {
      fet.obj@dat_list <- lapply(fet.obj@dat_list, function(dat)dat[kpt.idx,])
      validObject(fet.obj)
      loop.obj@loop <- loop.obj@loop[kpt.idx,]
      loop.obj@loop[["rowid"]] <- seq_len(nrow(loop.obj@loop))
      loop.obj@g <- delete.edges(loop.obj@g, which(!E(loop.obj@g)$loop %in% loop.obj@loop[["loop"]]))
      loop.obj@g <- delete.vertices(loop.obj@g, which(igraph::degree(loop.obj@g)<1))
      if (!is.null(loop.obj@split)) {
        loop.obj@split <- factor(as.character(loop.obj@split)[kpt.idx], levels = unique(as.character(loop.obj@split)[kpt.idx]))
      }
      validObject(loop.obj)
    }
    # row scale RNA data
    invisible(sapply(grep("RNA", fet.obj@hash[["sms"]], ignore.case = TRUE), function(i){
      fet.obj@dat_list[[i]] <<- data.table(t(scale(t(fet.obj@dat_list[[i]]))))
    }))
    validObject(fet.obj)
    return(list(loop.obj = loop.obj, fet.obj = fet.obj))
  }
)
