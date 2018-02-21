#!/usr/bin/env Rscript

###
# Input data file into JSON for read into H3K27AC_LoopType.R
# input: 
# (1) H3K27ac HiChiP file (minimum 4 columns: locus1, locus2, gene1, gene2): hichip
# (2) ChIP-seq file(s) to overlap with at least one of hi-chip anchors for loop validation: vchip
# (3) ChIP-seq file(s) to overlap at non-gene anchr for anchor to be called enhancers: echip
# (4) Output bed-like file for downstream analyses.
# output:
# (1) JSON file
###

library(RJSONIO)

readInput <- function(prompt){ 
	str <- readline(prompt = prompt)
	return(str)
}

fs <- c("hichip", "vchip", "echip", "bedout")
fs_list <- list()
i <- 1
repeat {
	string <- readInput(prompt = paste0(fs[i], " file: "))
	if (nchar(string) == 0) {
		i <- i + 1
		if (i > length(fs)) break
	} else {
		fs_list[[i]] <- if (length(fs_list) < i) {
			structure(string, class = "character")
		} else {
			c(fs_list[[i]], structure(string, class = "character"))
		}
	}
}
names(fs_list) <- fs

repeat {
	string <- readInput(prompt = "output configuration: ")
	if (nchar(string) > 0) {
		fout <- string
		break
	}
}

exportJson <- toJSON(fs_list)
write(exportJson, fout)
