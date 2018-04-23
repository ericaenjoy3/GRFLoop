#!/usr/bin/env Rscript

###
# LOLA directory
###

libfs <- c("LOLA", "data.table", "tidyverse", "circlize", "ComplexHeatmap", "argparse")
invisible(sapply(libfs, function(f)suppressPackageStartupMessages(require(f, character.only = T))))

parser <- ArgumentParser()
parser$add_argument("--din", type = "character", required = TRUE,
  help = "An input LOLA directory.")
parser$add_argument("--pdffout", type = "character", required = TRUE,
  help = "Output heatmap.")
parser$add_argument("--self", type = "character", required = FALSE,
  help = "Select a subset of factors.")
args <- parser$parse_args()
attach(args)

fs <- dir(path = din, pattern = "userSet_*", full.names = TRUE, recursive = TRUE)
sms <- gsub("userSet_(.+)\\.tsv", "\\1", basename(fs))

dat_list <- lapply(fs, function(f) {
	dd <- fread(f, header = TRUE)
	dd <- copy(dd[!(is.na(cellType) | is.na(antibody))])
	return(dd)
})

# signficant features
tmp <- rbindlist(lapply(dat_list, function(dd){
	dd[, `:=`(cellType = toupper(cellType), antibody = toupper(antibody))]
	dd[grepl("Embryonic Stem Cell", cellType, ignore.case = TRUE), cellType := toupper("Embryonic Stem Cell")]
	dd <- copy(dd[pValueLog >= -log10(0.05)])
	return(dd)
}))
com_dat <- tmp %>% group_by(cellType, antibody) %>% top_n(1, pValueLog) %>% filter(grepl("Embryonic Stem Cell", cellType, ignore.case = TRUE))  %>%data.table()

sel_vec <- fread(self, header = FALSE)[[1]]
ridx <- com_dat[toupper(antibody) %in% sel_vec, which = TRUE]
com_dat <- com_dat[ridx]

# heatmap of qvalue and logOddsRatio
stopifnot(identical(dat_list[[1]][chmatch(com_dat[, filename], filename), filename], com_dat[, filename]))

proc <- function(dat_list, col_nm, com_dat) {
	or_dat <- do.call("cbind", lapply(dat_list, function(dd){
		stopifnot(col_nm %in% colnames(dd))
		sm <- dd[, unique(userSet)]
		dd <- copy(dd[chmatch(com_dat[, filename], filename)])
		nd <- dd[, c(col_nm), with = FALSE]
		setnames(nd, col_nm, sm)
		return(nd)
	}))
}

heatPlot <- function(or_dat, pv_dat, com_dat, pdffout, pcol) {
	gen <- function(mat, com_dat, name) {
		mat <- mat %>% t() %>% scale() %>% t()
		rownames(mat) <- com_dat[, paste(cellType, antibody, sep = "_")]
		max.val <- ceiling(max(max(mat, na.rm = TRUE), abs(min(mat, na.rm = TRUE))))
		ht_list <- Heatmap(mat,
			col = colorRamp2(c(-max.val, 0, max.val), c("#4575B4", "white", "#D73027")),
			heatmap_legend_param = list(color_bar = "continuous"),
			show_row_dend = FALSE,
			row_names_side = "left",
			row_names_max_width = unit(10, "cm"),
			name = name,
			cluster_columns = FALSE,
			clustering_distance_rows = "spearman",
			clustering_method_rows = "average")
		return(ht_list)
	}
	ht_list <- gen(pv_dat, com_dat, name = paste("Scaled", pcol))
	ht_list <- ht_list + gen(or_dat, com_dat, name = paste("Scaled", "logOddsRatio"))
	pdf(pdffout, width = 7*2, height = 8*2)
	draw(ht_list)
	dev.off()
}



or_dat <- proc(dat_list, col_nm = "logOddsRatio", com_dat = com_dat)
write.table(data.table(com_dat[, .(cellType, antibody)], or_dat), 
	file = file.path(din, "logOddsRatio.txt"), 
	row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

pv_dat <- proc(dat_list, col_nm = "pValueLog", com_dat = com_dat)
write.table(data.table(com_dat[, .(cellType, antibody)], pv_dat), 
	file = file.path(din, "pValueLog.txt"), 
	row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

qv_dat <- proc(dat_list, col_nm = "qValue", com_dat = com_dat)
write.table(data.table(com_dat[, .(cellType, antibody)], qv_dat), 
	file = file.path(din, "qValue.txt"), 
	row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

qv_dat <- qv_dat[, lapply(.SD, function(vec)-log10(vec + 10^(-100)))]

heatPlot(or_dat, qv_dat, com_dat, pdffout, pcol = "-qValueLog")

