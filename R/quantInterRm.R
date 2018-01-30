#' @include GRFLoopClass.R

#' @export quantInterRm
setGeneric(name = "quantInterRm",
  def = function(loop.obj, fet.obj, dedup){
    standardGeneric("quantInterRm")
  }
)

#' @rdname quantInterRm-methods
setMethod(f = "quantInterRm",
  signature = c("loop", "fet", "logical"),
  definition = function(loop.obj, fet.obj, dedup){
    score <- E(loop.obj@g)$score
    names(score) <- E(loop.obj@g)$loop
    if (!dedup) {
      score <- score[loop.obj@loop[["loop"]]]
    }
    thresh <- as.numeric(quantile(score, probs = c(0.25, 0.75)))
    bottom_loop <- E(loop.obj@g)$loop[E(loop.obj@g)$score <= thresh[1]]
    top_loop <- E(loop.obj@g)$loop[E(loop.obj@g)$score >= thresh[2]]
    # update g slot of loop object
    loop.obj@g <- delete.edges(loop.obj@g, which(!E(loop.obj@g)$loop %in% c(bottom_loop, top_loop)))
    loop.obj@g <- delete.vertices(loop.obj@g, which(igraph::degree(loop.obj@g)<1))
    # update loop slot of loop object
    kpt.idx <- which(loop.obj@loop[["loop"]] %in% c(bottom_loop, top_loop))
    stopifnot(length(kpt.idx) %between% c(0, nrow(loop.obj@loop)))
    loop.obj@loop <- loop.obj@loop[kpt.idx, ]
    loop.obj@loop[["rowid"]] <- seq_len(nrow(loop.obj@loop))
    split <- rep(0, length(loop.obj@loop[["loop"]]))
    split[loop.obj@loop[["loop"]] %in% top_loop] <- 1
    loop.obj@split <- factor(split, levels = unique(split))
    validObject(loop.obj)
    # update dat_list slot of fet object
    fet.obj@dat_list <- lapply(fet.obj@dat_list, function(dat)dat[kpt.idx,])
    validObject(fet.obj)
    return(list(loop.obj = loop.obj, fet.obj = fet.obj))
  }
)
