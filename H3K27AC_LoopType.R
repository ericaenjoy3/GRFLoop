#!/usr/bin/env Rscript

###
# Processing H3K27AC Hi-ChIP Interaction Data
# Filters:
# (1) Remove duplicated loops
# (2) Transform gene id into gencode IDs
# (3) Prom-Prom, Prom-Enh Looping, Enh-Enh Looping by overlaps with K27ac specified.
###

libfs <- c("data.table", "tidyverse", "igraph", "RJSONIO", "ggplot2", "argparse")
sapply(libfs, library, character.only = TRUE)
options(scipen = 999)

parser <- ArgumentParser()
parser$add_argument("--hichip", type = "character", required = FALSE,
  help = "A query loop file with columns being locus 1 and 2, genes 1 and 2 and optional columns.")
parser$add_argument("--vchip", type = "character", required = FALSE,
  help = "h3k27ac chip file to determine enhancers.")
parser$add_argument("--echip", type = "character", required = FALSE,
  help = "enhancers coordinate file to overlap anhcors for determining enhancers.")
parser$add_argument("--intsec", action="store_true",
  help = "anchors overlapped by multiple chipf to be called enhancers.")
parser$add_argument("--bedout", type = "character", required = FALSE,
  help = "Output bed-like file of interactions.")
parser$add_argument("--json", type = "character", required = FALSE,
  help = "json configuration file.")
args <- parser$parse_args()

if (!"json" %in% names(args)) {
	args <- c(args["intsec"], fromJSON(args[["json"]]))
}
stopifnot(all(c("hichip", "vchip", "bedout") %in% names(args)))
attach(args)

dat <- fread(hichip, header = TRUE, sep = "\t", na.strings="N/A")
setnames(dat, colnames(dat)[1:4], c("loc1", "loc2", "gene1", "gene2"))
dat <- dat %>% separate(
	loc1, c("loc1Chr", "loc1Start", "loc1End"), sep = "[:-]", remove = TRUE, convert = TRUE) %>% separate(
	loc2, c("loc2Chr", "loc2Start", "loc2End"), sep = "[:-]", remove = TRUE, convert = TRUE)

# (1) Remove duplicated loops
stopifnot(all(dat[, loc1Chr == loc2Chr]))
idx <- dat[loc1Chr == loc2Chr & loc1Start > loc2Start, which = TRUE]
if (length(idx) > 0) {
	dat <- dat[idx, c("loc1Start", "loc1End", "loc2Start", "loc2End", "gene1", "gene2") := .(loc2Start, loc2End, loc1Start, loc1End, gene2, gene1)] %>% unique()
}

# (2) Gene1 and 2 contains observations with multiple delimited values, this separates the values and places each one in its own row.
dat <- separate_rows(dat, gene1, sep = ",", convert = TRUE) %>% separate_rows(gene2, sep = ",", convert = TRUE)

# (3) Transform gene id into gencode IDs
gencodef <- "/home/liuyiyua/athena/Gencode/mm10/annotation/gencode.vM6.annotation.gene.bed"

gencode <- fread(gencodef, header = FALSE, sep = "\t")[, 4, with = FALSE] %>% 
	setnames(c("gid")) %>%
	separate(gid, c("id", "gname", "type"), sep = "\\|", remove = FALSE, convert = TRUE) %>%
	select(gid, gname)

dat <- dat %>% 
	merge(gencode, by.x = "gene1", by.y = "gname", all.x = TRUE, sort = FALSE) %>% 
	setnames("gid", "g1") %>%
	merge(gencode, by.x = "gene2", by.y = "gname", all.x = TRUE, sort = FALSE) %>%
	setnames("gid", "g2") %>%  
	unique()
dat <- dat[!(!is.na(gene1) & is.na(g1)) | (!is.na(gene2) & is.na(g2))]
set(dat, j = c("gene1", "gene2"), value = NULL)
setnames(dat, c("g1", "g2"), c("gene1", "gene2"))
dat[, loop := paste0(loc1Chr, ":", loc1Start, "-", loc1End, "_", loc2Chr, ":", loc2Start, "-", loc2End)]

# (4) Filter by at least one of the two anchor overlapped by H3K27AC ChIP-seq
anchorOlap <- function(dat, fs) {
	chip_list <- lapply(seq_along(fs), function(j){
		chip <- fread(fs[j], header = FALSE)[, 1:3] %>%
			setnames(c("chr", "start", "end"))
		chip[, start := start + 1]
		setkeyv(chip, colnames(chip))
		return(chip)
	})
	matchloc_list <- lapply(chip_list, function(chip){
		setkeyv(dat, c("loc1Chr", "loc1Start", "loc1End"))
		matchLoc1Loop <- dat[unique(foverlaps(dat, chip, nomatch = 0, which = TRUE)$xid), loop]
		setkeyv(dat, c("loc2Chr", "loc2Start", "loc2End"))
		matchLoc2Loop <- dat[unique(foverlaps(dat, chip, nomatch = 0, which = TRUE)$xid), loop]
		return(list(matchLoc1Loop = matchLoc1Loop, matchLoc2Loop = matchLoc2Loop)) 
	})
	return(matchloc_list)
}

vmatchloc_list <- anchorOlap(dat, vchip)
matchloc <- unlist(vmatchloc_list, use.names = FALSE)
dat <- dat[loop %in% matchloc]

# (4) Prom-Prom, Prom-Enh Looping, Enh-Enh Looping
if ("echip" %in% names(args) & !is.null(args[["echip"]])) {
	ematchloc_list <- anchorOlap(dat, echip)
} else {
	ematchloc_list <- vmatchloc_list
}
ematchloc_spec <- if (intsec & length(ematchloc_list) > 1) {
	func <-ifelse(intsec & length(ematchloc_list) >1, "intersect", "union")
	message(func, " on ematchloc_list")
	maxlen <- max(sapply(ematchloc_list, length))
	lapply(seq(maxlen),function(i) Reduce(func,lapply(ematchloc, "[[", i))) %>% 
		unlist()
} else {
	ematchloc_list[[1]]
}

dat <- dat[, c("loc1type", "loc2type") := list(ifelse(!is.na(gene1), "Prom", ifelse(loop %in% ematchloc_spec[[1]], "Enh", NA)),
	ifelse(!is.na(gene2), "Prom", ifelse(loop %in% ematchloc_spec[[2]], "Enh", NA)))] %>%
	subset(subset = !(is.na(loc1type) | is.na(loc2type)))

dat[loc1type == "Enh" & loc2type == "Prom", 
	c("loc1Chr", "loc1Start", "loc1End", "loc2Chr", "loc2Start", "loc2End", "gene1", "gene2", "loc1type", "loc2type") := 
	.(loc2Chr, loc2Start, loc2End, loc1Chr, loc1Start, loc2End, gene2, gene1, loc2type, loc1type)]
fwrite(dat[, c("loc1Chr", "loc1Start", "loc1End", "loc2Chr", "loc2Start", "loc2End", "gene1", "gene2", "loc1type", "loc2type")], 
	file = bedout, 
	sep = "\t", 
	row.names = FALSE,
	col.names = TRUE, 
	quote = FALSE)



