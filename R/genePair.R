#' @include GRFLoopClass.R GRFLoopGeneric.R

#' @export genePair
setGeneric(name = "genePair",
  def = function(loop.obj, info.obj, type, nmin, nmax, random.it = FALSE){
    standardGeneric("genePair")
  }
)

#' @rdname genePair-methods
setMethod(f = "genePair",
  signature = c("loop", "info"),
  definition = function(loop.obj, info.obj, type, nmin, nmax, random.it = TRUE) {

  	glist2gpair <- function(glist, mess = NULL) {
  		# glist: split by promoter fragments
  		if (!is.null(mess)) {
  			message(mess, " from gene list to gene pair")
  		}

	    # cross connectom gene pairs
	    mat <- combn(1:length(glist), 2)
	    gpair_dat <- rbindlist(lapply(1:ncol(mat), function(j){
	    	return(expand.grid(glist[[mat[1, j]]], glist[[mat[2, j]]], stringsAsFactors = FALSE))
	    }))
	    
	    # remove identical gene at gene pair ends
	    rm_ends <- gpair_dat[, Var1 == Var2]

	    # filter for unique gene pairs
	    rm_pairs <- duplicated(gpair_dat)

	    if (sum(rm_ends | rm_pairs) == nrow(gpair_dat)) {
        return(NA)
      }

	    gpair_dat <- gpair_dat[!(rm_ends|rm_pairs)]

	    return(gpair_dat)
	  }

    # identify incident loop
    stopifnot(sum(V(loop.obj@g)$vtype == type) > 0)
    ve <- V(loop.obj@g)$name[V(loop.obj@g)$vtype == type]
    ed <- incident_edges(loop.obj@g, ve)
    kpt_num <- sapply(ed, length) >= nmin & sapply(ed, length) <= nmax

    # (genunine) extract PromGene from loop slot of loop.obj
    gene_list <- list()
    gene_list[[1]] <- lapply(ed[which(kpt_num)], function(es, loop_hash){
      # potentially multiple loops => mutiple fragments
      lp <- as_ids(es)
      dd <- loop_hash[loop %in% lp]
      stopifnot(all(dd[, all(!is.na(gene1)), by = loop][, V1]))
      dd <- dd[, unique(gene1), by = loop]
      gs <- split(dd[, V1], dd[, loop])
      names(gs) <- NULL
      return(gs)
    }, loop_hash = loop.obj@loop)
    names(gene_list[[1]]) <- NULL

    # (genunine) unique gene sets for connectomes
    kpt_gset <- !duplicated(gene_list[[1]])
    message(sum(kpt_gset), " hubs removed due to duplication in idnetical genes contacted overall.")

    # (genunine) final gene_list
    gene_list[[1]] <- gene_list[[1]][kpt_gset]
    ve <- ve[which(kpt_num)][which(kpt_gset)]
    ed <- ed[which(kpt_num)][which(kpt_gset)]

    # (genuine) unique 
    gpair_list <- list()
    g1_list <- lapply(gene_list[[1]], function(glist) {
      gpair <- glist2gpair(glist)
      return(gpair)
    })
    g1_list <- g1_list[!is.na(g1_list)]
    gpair_list[[1]] <- rbindlist(g1_list)


    if (!random.it) {
      return(gene_list)
    }

    # coordinates of genes
    coord <- rbindlist(lapply(gene_list[[1]], function(gs, info.obj){
        gs <- unique(unlist(gs))
        idx <- chmatch(gs, info.obj@gene[["gene"]])
        stopifnot(all(!is.na(idx)))
        unique(info.obj@gene[idx][, list(chr = chr, start = min(start), end = max(end))])
    }, info.obj = info.obj))

    # (global random) 
    gene_list[[2]] <- coordShulf(coord, info.obj, dout, nshulf = 100, nmin = nmin, nmax = nmax)
    
    # (TAD matched random)
    gene_list[[3]] <- inTADShulf(gene_list[[1]], info.obj)

    for (i in 2:3) {
    	message(ifelse(i==2, "Global random", "TAD matched random"))
    	gpair_list[[i]] <- rbindlist(lapply(gene_list[[i]], function(vec)glist2gpair(as.list(vec))))
    }

    return(gpair_list)
  }
 )