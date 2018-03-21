#!/usr/bin/env Rscript

###
# LOLA on input bed
###

libfs <- c("LOLA", "data.table", "GenomicRanges", "argparse")
invisible(sapply(libfs, function(f)suppressPackageStartupMessages(require(f, character.only = T))))

parser <- ArgumentParser()
parser$add_argument("--bedin", type = "character", required = TRUE,
  help = "An input bed file.")
parser$add_argument("--dout", type = "character", required = TRUE,
  help = "Output LOLA result directory.")
parser$add_argument("--bgbed", type = "character", required = FALSE,
  help = "Background bed for LOLA (optional).")
args <- parser$parse_args()
attach(args)

bedin_dat <- fread(bedin, header = FALSE)
bedin_dat[, V2 := V2 + 1]
userSets <- makeGRangesFromDataFrame(bedin_dat,
	keep.extra.columns = TRUE,
	seqnames.field = "V1", 
	start.field = "V2",
	end.field = "V3",
)
userSets <- split(userSets, mcols(userSets))

if (!is.null(bgbed)) {
	message("Background bed supplied: ", bgbed)

	bedbg_dat <- fread(bedbg, header = FALSE)
	bedbg_dat[, V2 := V2 +1] 
	userUniverse <- makeGRangesFromDataFrame(bedbg_dat,
		keep.extra.columns = TRUE,
		seqnames.field = "V1", 
		start.field = "V2",
		end.field = "V3",
	)

} else {
	message("Using tiled genome as background")

	chrlenf <- "/home/liuyiyua/athena/Gencode/mm10/sequence/chrNameLengthKnownChr.sizes"
	chrlen_dat <- fread(chrlenf, header = FALSE)
	seqlengths <- chrlen_dat[, V2]
	names(seqlengths) <- chrlen_dat[, V1]
	userUniverse <- tileGenome(seqlengths, tilewidth = 10000, cut.last.tile.in.chrom = TRUE)

}

dbPath <- file.path("~/athena/LOLA/mm10")
regionDB <- loadRegionDB(dbPath)

res <- runLOLA(userSets, userUniverse, regionDB, cores=1)
res <- res[order(userSet, pValueLog, decreasing=TRUE),]

dir.create(dout, showWarnings = FALSE, recursive = TRUE)
writeCombinedEnrichment(res, outFolder= dout, includeSplits = TRUE)