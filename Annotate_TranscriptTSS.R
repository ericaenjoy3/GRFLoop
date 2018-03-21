#!/usr/bin/env Rscript

###
# Re-annotate Transcript TSSs overlapped with anchors of H3K27AC Hi-ChIP loops
# Input: original Andreas's H3K27AC Hi-ChIP loop data (1-based coordinates)
# Output: updated gene annotation (1-based coordinates)
###

libfs <- c("data.table", "tidyverse", "gtools", "igraph", "argparse")
invisible(sapply(libfs, function(f)suppressPackageStartupMessages(require(f, character.only = T))))
options(scipen = 999)

parser <- ArgumentParser()
parser$add_argument("--fin", type = "character", required = FALSE,
  help = "An input K27ac loop file.")
parser$add_argument("--tssbed", type = "character", required = FALSE,
  help = "An input tss bed file.")
parser$add_argument("--fout", type = "character", required = FALSE,
  help = "An output K27ac loop file.")
args <- parser$parse_args()
attach(args)

dat <- fread(fin, header = TRUE, sep = "\t", na.strings=c("N/A", ""))
ori_colmn <- colnames(dat)

if (all(c("locus.1", "locus.2") %in% colnames(dat))) {
	loc_stand <- "two"
	dat <- dat %>% separate(
	locus.1, c("loc1Chr", "loc1Start", "loc1End"), remove = TRUE, convert = TRUE) %>% separate(
	locus.2, c("loc2Chr", "loc2Start", "loc2End"), remove = TRUE, convert = TRUE)
} else {
	loc_stand <- "one"
	dat <- dat %>% separate(
	locus, c("loc1Chr", "loc1Start", "loc1End", "loc2Chr", "loc2Start", "loc2End"), remove = TRUE, convert = TRUE)
}

# retrieve original loop for ordering reference
dat[, loop := paste0(loc1Chr, ":", loc1Start, "-", loc1End, "_", loc2Chr, ":", loc2Start, "-", loc2End)]
stopifnot(dat[, sum(duplicated(loop))==0])
ori_loop <- copy(dat[["loop"]])

# overlap of transcript TSSs
bed <- fread(tssbed, header = FALSE)[, 1:4, with = FALSE]
setnames(bed, c("chr", "start", "end", "gid"))
bed[, start:=start + 1]
setkeyv(bed, colnames(bed)[1:3])

overlap_proc <- function(dat, locSel = c("loc1", "loc2")) {
	locSel <- match.arg(locSel)
	setkeyv(dat, paste0(locSel, c("Chr", "Start", "End")))
	idx_tbl <- foverlaps(dat, bed, nomatch = 0, which = TRUE)
	gdat <- data.table(xid = idx_tbl[["xid"]], bed[idx_tbl[["yid"]], "gid"]) %>% unique()
	gdat <- gdat[, paste(gid, collapse = ","), by = xid]
	setnames(gdat, "V1", "gid")
	set(dat, i = gdat[["xid"]], j = paste0(locSel, ".genes"), value = gdat[["gid"]])
	return(dat)
}

dat <- overlap_proc(dat, locSel = "loc1")
dat <- overlap_proc(dat, locSel = "loc2")

# replace old columns with new
dat[, `:=`(locus1.genes = loc1.genes, locus2.genes = loc2.genes)][, `:=`(loc1.genes = NULL , loc2.genes = NULL)]

# matching row order as the original
idx <- chmatch(ori_loop, dat[["loop"]])
stopifnot(identical(dat[idx, loop], ori_loop))
dat <- dat[idx]

if (loc_stand == "two") {
	dat[, `:=`(locus.1 = paste0(loc1Chr, ":", loc1Start, "-", loc1End), 
		locus.2 = paste0(loc2Chr, ":", loc2Start, "-", loc2End))]
	ess_colmn <- c("locus.1", "locus.2", "locus1.genes", "locus2.genes")
} else if (loc_stand == "one") {
	dat[, `:=`(locus = loop)]
	ess_colmn <- c("locus", "locus1.genes", "locus2.genes") 

}

# select and order columns
ori_colmn <- c(ess_colmn, ori_colmn[!ori_colmn %in% ess_colmn])
dat <- dat[, ori_colmn, with = FALSE]
setcolorder(dat, ori_colmn)

write.table(dat, file = fout, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")



