#' @include GRFLoopClass.R

#' @export quantConRm
setGeneric(name = "quantConRm",
  def = function(loop.obj, fet.obj, vType = c("Prom", "Enh")){
    standardGeneric("quantConRm")
  }
)

#' @rdname quantInterRm-methods
setMethod(f = "quantConRm",
  signature = c("loop", "fet"),
  definition = function(loop.obj, fet.obj, vType = c("Prom", "Enh")){
    vType <- match.arg(vType)
    stopifnot(vType %in% unique(V(loop.obj@g)$type))
  	# incident edges of a particular vertex type
    ve <- V(loop.obj@g)$name[V(loop.obj@g)$type == vType]
    ed <- incident_edges(loop.obj@g, ve)
    ed_num <- sapply(ed, length)
    thresh <- as.numeric(quantile(ed_num, probs = c(0.30, 0.90)))
    message("top thresh: ", thresh[2], "; bottom thresh:", thresh[1])
    stopifnot(!identical(thresh[1], thresh[2]))
    bottom_loop <- unlist(lapply(ed[ed_num <= thresh[1]], function(e){
    	vec<-gsub("\\|", "_", as_ids(e)); 
    	return(vec)
    }))
    top_loop <- unlist(lapply(ed[ed_num >= thresh[2]], function(e){
    	vec<-gsub("\\|", "_", as_ids(e)); 
    	return(vec)
    }))
    names(bottom_loop) <- names(top_loop) <- NULL
    # update g slot of loop object
    loop.obj@g <- delete.edges(loop.obj@g, which(!E(loop.obj@g)$loop %in% c(bottom_loop, top_loop)))
    loop.obj@g <- delete.vertices(loop.obj@g, which(igraph::degree(loop.obj@g) < 1))
    # update loop slot of loop object
    kpt.idx <- which(loop.obj@loop[["loop"]] %in% c(bottom_loop, top_loop))
    stopifnot(length(kpt.idx) %between% c(0, nrow(loop.obj@loop)))
    loop.obj@loop <- loop.obj@loop[kpt.idx, ]
    loop.obj@loop[["rowid"]] <- seq_len(nrow(loop.obj@loop))
    split <- rep(0, length(loop.obj@loop[["loop"]]))
    split[loop.obj@loop[["loop"]] %in% top_loop] <- 1
    message(sum(split ==0), " bottom loops")
     message(sum(split ==1), " top loops")   
    loop.obj@split <- factor(split, levels = unique(split))
    validObject(loop.obj)
    # update dat_list slot of fet object
    fet.obj@dat_list <- lapply(fet.obj@dat_list, function(dat)dat[kpt.idx,])
    validObject(fet.obj)
    return(list(loop.obj = loop.obj, fet.obj = fet.obj))
  }
)