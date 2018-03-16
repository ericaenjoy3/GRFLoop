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

	    stopifnot(sum(rm_ends | rm_pairs) < nrow(gpair_dat))
	    gpair_dat <- gpair_dat[!(rm_ends|rm_pairs)]

	    return(gpair_dat)
	  }

    # identify incident loop
    stopifnot(sum(V(loop.obj@g)$vtype == type) > 0)
    ve <- V(loop.obj@g)$name[V(loop.obj@g)$vtype == type]
    ed <- incident_edges(loop.obj@g, ve)
    kpt_num <- sapply(ed, length) >= nmin & sapply(ed, length) <= nmax

    # (genunine) testing for more than 1 gene pair

    # (genunine) extract PromGene from loop slot of loop.obj
    gene_list <- list()
    gene_list[[1]] <- lapply(ed[which(kpt_num)], function(es, loop_hash){
      lp <- as_ids(es)
      gs <- copy(loop_hash[loop %in% lp, unique(gene1)])
      stopifnot(all(!is.na(gs)))
      return(gs)
    }, loop_hash = loop.obj@loop)

    # (genunine) unique gene sets for connectomes
    kpt_gset <- !duplicated(gene_list[[1]])

    # (genunine) final gene_list
    gene_list[[1]] <- gene_list[[1]][kpt_gset]
    ve <- ve[which(kpt_num)][which(kpt_gset)]
    ed <- ed[which(kpt_num)][which(kpt_gset)]

    gpair_list <- list()
    gpair_list[[1]] <- rbindlist(lapply(ed, function(es, loop_hash){
      lp <- as_ids(es)
      message(lp)
      dd <- loop_hash[loop %in% lp, .(loop, gene1)]
      glist <- split(dd[, gene1], dd[, loop])
      gpair <- glist2gpair(glist)
      return(gpair)
    }, loop_hash = loop.obj@loop))

    if (!random.it) {
      return(gene_list)
    }

    # coordinates of genes
    coord <- rbindlist(lapply(gene_list[[1]], function(gs, info.obj){
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