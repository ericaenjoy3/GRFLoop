#!/usr/bin/env Rscript

###
# Processing KLF4 Peak Classification into Bi-anchor/One-anchor/No-anchor
# Filters:
# (1) Remove duplicated loops
# (2) Filter hichip chr by chr in chip-seq peaks
# (3) igraph
###

libfs <- c("data.table", "tidyverse", "igraph", "argparse")
invisible(sapply(libfs, function(f)suppressPackageStartupMessages(require(f, character.only = T))))
options(scipen = 999)

# hichip: 1-based coordinate
# chip: 0-based coordinate
# fout: 0-based coordinate

parser <- ArgumentParser()
parser$add_argument("--hichip", type = "character", required = FALSE,
  help = "A query loop file with columns being locus 1 and 2, genes 1 and 2 and optional columns.")
parser$add_argument("--chip", type = "character", nargs='+', required = FALSE,
  help = "Chip-seq peaks of KLF4.")
parser$add_argument("--fout", type = "character", nargs='+', required = FALSE,
  help = "Output bed file with the 4th column indicating the anchor status.")
args <- parser$parse_args()
attach(args)


dat <- fread(hichip, header = TRUE, sep = "\t", na.strings=c("N/A", ""))[, 1:2, with = FALSE]
setnames(dat, colnames(dat), c("loc1", "loc2"))
dat <- dat %>% separate(
	loc1, c("loc1Chr", "loc1Start", "loc1End"), sep = "[:-]", remove = TRUE, convert = TRUE) %>% separate(
	loc2, c("loc2Chr", "loc2Start", "loc2End"), sep = "[:-]", remove = TRUE, convert = TRUE)
# dat[, `:=`(loc1Start = loc1Start - 1, loc2Start = loc2Start -1)]
dat[, loc1 := paste0(loc1Chr, ":", loc1Start, "-", loc1End)]
dat[, loc2 := paste0(loc2Chr, ":", loc2Start, "-", loc2End)]

# (1) Remove duplicated loops
stopifnot(all(dat[, loc1Chr == loc2Chr]))
idx <- dat[loc1Chr == loc2Chr & loc1Start > loc2Start, which = TRUE]
if (length(idx) > 0) {
	n_before <- nrow(dat)
	dat <- dat[idx, c("loc1Start", "loc1End", "loc2Start", "loc2End", "gene1", "gene2") := .(loc2Start, loc2End, loc1Start, loc1End, gene2, gene1)] %>% unique()
	n_after <- nrow(dat)
	message(n_after - n_before, " duplicated loops removed.")
}

# Chip-seq peaks
chip_dat <- fread(chip, header = FALSE)[, 1:3, with = FALSE]
setnames(chip_dat, c("chr", "start", "end"))
chip_dat[, start := start + 1]
setkeyv(chip_dat, c("chr", "start", "end"))
chip_dat[, rowid := 1:nrow(chip_dat)]

# (2) Filter hichip chr by chr in chip-seq peaks
message(dat[!loc1Chr %in% chip_dat[, unique(chr)], which = TRUE], " loops removed due to chr not present in the chip-seq file.")
dat <- dat[loc1Chr %in% chip_dat[, unique(chr)],]

# (3) igraph
e_dat <- dat[, .(loc1, loc2)]
v_dat <- rbind(dat[, .(loc1, loc1Chr, loc1Start, loc1End)], dat[, .(loc2, loc2Chr, loc2Start, loc2End)], use.names = FALSE) %>% unique()
setnames(v_dat, gsub("1", "", colnames(v_dat)))
setkeyv(v_dat, c("locChr", "locStart", "locEnd"))
idx_dat <- foverlaps(v_dat, chip_dat, nomatch = 0, which = TRUE)
idx_dat <- idx_dat[, paste(yid, collapse = ";"), by = xid]
v_dat[idx_dat[["xid"]], rowid := idx_dat[["V1"]]]
g <- graph_from_data_frame(e_dat, directed = FALSE, vertices = v_dat)
vp_rowid <- V(g)$rowid[which(!is.na(V(g)$rowid))]
v_list <- lapply(adjacent_vertices(g, which(!is.na(V(g)$rowid))), function(v)v$rowid)
lab <- sapply(v_list, function(vec){
	if(any(!is.na(vec))) {
		return("Bi-anchor")
	} else if (all(is.na(vec))) {
		return("One-anchor")
	}
})
ndat <- unique(data.table(rowid = vp_rowid, lab = lab))
ndat <- separate_rows(ndat, rowid, sep = ";", convert = TRUE)
ndat <- ndat[, ifelse(length(lab) >1, "Mixed", lab), by = rowid] # rowid =58 mixed lab
setnames(ndat, "V1", "lab")

res_dat <- merge(chip_dat, ndat, by.x = "rowid", by.y= "rowid", sort = FALSE, all.x = TRUE)
res_dat[is.na(lab), lab := "No-anchor"]
res_dat[, start :=  start -1]
write.table(res_dat[, !"rowid"], file = fout, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")