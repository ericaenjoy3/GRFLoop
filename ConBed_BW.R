#!/usr/bin/env Rscript

###
# Bigwig signals under peaks that overlapped with large windows
###

libfs <- c("data.table", "tidyverse", "gtools", "argparse")
invisible(sapply(libfs, function(f)suppressPackageStartupMessages(require(f, character.only = TRUE))))
options(scipen = 999)

parser <- ArgumentParser()
parser$add_argument("--window", type = "character", required = FALSE,
  help = "A bed file of large windows (> 1 kb).")
parser$add_argument("--peak", type = "character", nargs='+', required = FALSE,
  help = "A peak file of peaks.")
parser$add_argument("--config", type = "character", required = FALSE,
  help = "Bigwig configuration files.")
parser$add_argument("--matf", type = "character", required = FALSE,
  help = "Output signal matrix file.")
args <- parser$parse_args()
attach(args)

# create temporary directory
dout <- file.path(dirname(config), "tmp")
dir.create(dout, showWarnings = FALSE, recursive = TRUE)

# intersect peak file by window
ow_dat <- fread(window, header = FALSE)[, 1:4, with = FALSE]
w_dat <- ow_dat[!duplicated(ow_dat)]
setkeyv(w_dat, paste0("V", 1:3))
ow2w_tbl <- data.table(xid = 1:nrow(ow_dat), yid = chmatch(ow_dat[, paste(V1, ":", V2, "-", V3)], w_dat[, paste(V1, ":", V2, "-", V3)]))

p_dat <- fread(peak, header = FALSE)[, 1:3, with = FALSE]
setkeyv(p_dat, paste0("V", 1:3))

# stopifnot(all(w_dat[['V1']] %in% p_dat[['V1']]) & all(p_dat[['V1']] %in% w_dat[['V1']]))

# overlap intersected peaks files by windows
idx_tbl <- foverlaps(p_dat, w_dat, which = TRUE, nomatch = 0)
pint_dat <- data.table(p_dat[idx_tbl$xid, 1:3, with = FALSE], w_dat[idx_tbl$yid, 1:3, with = FALSE])
setnames(pint_dat, 4:6, paste0("V", 4:6))

idx <- pint_dat[V2 < V5, which = TRUE]
pint_dat[idx, V2 := V5]

idx <- pint_dat[V3 > V6, which = TRUE]
pint_dat[idx, V3:= V6]

idx <- pint_dat[V2==V3, which = TRUE]
pint_dat <- pint_dat[!1:nrow(pint_dat) %in% idx]
idx_tbl <- idx_tbl[!1:nrow(idx_tbl) %in% idx]

tmp_pint <- paste0(tempfile(tmpdir = dout), ".bed")
write.table(pint_dat[, 1:3, with = FALSE], file = tmp_pint, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# average signals under intersected peaks
c_dat <- fread(config, header = FALSE)
cmd <- paste("export NSLOTS=1; ~/athena/ComScripts/CHIP/BW2matrix.sh -b", tmp_pint, "-f", config, "-s -1 -w 0")
system(cmd, intern = FALSE)

# heatmap matrix
mat_dat <- rbindlist(lapply(1:nrow(c_dat), function(i){
	dd <- fread(c_dat[i, V3], header = FALSE)
	stopifnot(nrow(dd) == nrow(idx_tbl))
	dd <- data.table(dd, yid = idx_tbl$yid, w_width = p_dat[idx_tbl$xid, V3-V2])
	nd <- dd[, sum(V1*w_width)/sum(w_width), by = yid]
	nd <- merge(ow2w_tbl, nd, by.x = "yid", by.y = "yid", all.x = TRUE, sort = FALSE)[, !"yid"]
	nd[is.na(V1), V1:=0]
	write.table(nd[,!"xid"], file = c_dat[i, V3], row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
	setnames(nd, 2, c_dat[i, V2])
	return(nd[,!"xid"])
}))

write.table(mat_dat, file = matf, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

unlink(dout, recursive = TRUE, force = TRUE)



