glist2gpair <- function(glist, mess = NULL) {
	# glist: split by promoter fragments
  	if (!is.null(mess)) {
		message(mess, " from gene list to gene pair")
	}

	# cross connectom gene pairs
	mat <- combn(1:length(glist), 2)

	gidx_dat <- rbindlist(lapply(1:ncol(mat), function(j){
		gvec1 <- glist[[mat[1, j]]]
		gvec2 <- glist[[mat[2, j]]]
		return(expand.grid(paste(mat[1,j], 1:length(gvec1), sep = "_"), paste(mat[2,j], 1:length(gvec2), sep = "_"), stringsAsFactors = FALSE))
	}))

	gpair_dat <- rbindlist(lapply(1:ncol(mat), function(j){
		return(expand.grid(glist[[mat[1, j]]], glist[[mat[2, j]]], stringsAsFactors = FALSE))
	}))
	    
	# remove identical gene at gene pair ends
	rm_ends <- copy(gpair_dat[, Var1 == Var2])

	# filter for duplication of Var1 and Var2 swap
	rm_swap <- duplicated(lapply(1:nrow(gpair_dat),function(i){
		vec <- as.vector(gpair_dat[i])
		names(vec) <- NULL
		mixedsort(vec)
	}))

	# filter for unique gene pairs
	rm_pairs <- copy(duplicated(gpair_dat))

	stopifnot(sum(rm_ends | rm_swap | rm_pairs) < nrow(gpair_dat))

	gpair_dat <- gpair_dat[!(rm_ends|rm_swap|rm_pairs)]
	gidx_dat <- gidx_dat[!(rm_ends|rm_swap|rm_pairs)]

	return(list(gidx_dat = gidx_dat, gpair_dat = gpair_dat))
}