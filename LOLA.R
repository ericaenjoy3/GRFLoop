#!/usr/bin/env Rscript

###
# LOLA on input bed
# Process:
###

libfs <- c("LOLA", "data.table", "GenomicRanges", "argparse")
invisible(sapply(libfs, function(f)suppressPackageStartupMessages(require(f, character.only = T))))

parser <- ArgumentParser()
parser$add_argument("--bedin", type = "character", required = FALSE,
  help = "An input bed file.")
parser$add_argument("--fout", type = "character", required = FALSE,
  help = "Output LOLA results.")
# parser$add_argument("--peak", type = "character", nargs='+', required = FALSE,
#   help = "A peak file of peaks.")
# parser$add_argument("--config", type = "character", required = FALSE,
#   help = "Bigwig configuration files.")
# parser$add_argument("--matf", type = "character", required = FALSE,
#   help = "Output signal matrix file.")
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

# bedbg <- file.path("~/athena/PlathChrom", "KPlath_chrom_mm10.bed")
# bedbg_dat <- fread(bedbg, header = FALSE)
# userUniverse <- makeGRangesFromDataFrame(bedbg_dat,
# 	keep.extra.columns = TRUE,
# 	seqnames.field = "V1", 
# 	start.field = "V2",
# 	end.field = "V3",
# )

chrlenf <- "/home/liuyiyua/athena/Gencode/mm10/sequence/chrNameLengthKnownChr.sizes"
chrlen_dat <- fread(chrlenf, header = FALSE)
seqlengths <- chrlen_dat[, V2]
names(seqlengths) <- chrlen_dat[, V1]
userUniverse <- tileGenome(seqlengths, tilewidth = 10000, cut.last.tile.in.chrom = TRUE)

dbPath <- file.path("~/athena/LOLA/mm10")
regionDB <- loadRegionDB(dbPath)

res <- runLOLA(userSets, userUniverse, regionDB, cores=1)
class(res$userSet) <- "character"
res[, userSet := gsub(".bed", "", basename(bedin))]
res <- res[order(pValueLog, decreasing=TRUE),]
#res <- res[pValueLog!="Inf"]
dir.create(dirname(fout), showWarnings = FALSE, recursive = TRUE)
fwrite(res, file = fout, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")