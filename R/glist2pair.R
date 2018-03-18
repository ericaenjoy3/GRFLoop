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