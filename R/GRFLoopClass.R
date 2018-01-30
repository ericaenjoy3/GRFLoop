#####rlang .data prevents R CMD check from giving a NOTE about undefined global variables
#' @import ModClusterR
#' @import RColorBrewer
#' @import ggplot2
#' @import dplyr
#' @import methods
#' @importFrom data.table fread rbindlist setnames melt data.table .I .N .SD
#' @importFrom rlang .data
#' @importFrom grDevices dev.off pdf png colorRampPalette
#' @importFrom stats na.omit quantile median
#' @importFrom graphics plot axis abline text par mtext
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom EnrichedHeatmap EnrichedHeatmap '+.AdditiveUnit'
#' @importFrom circlize colorRamp2
#' @importFrom utils read.table write.table setTxtProgressBar

#' @title An S4 class to represent loops.
#' @name loop-class
#' @rdname loop-class
#' @description Store cluster sets in a bed object.
#' @slot igraph data.frame object
#' @exportClass fet
inlibs <- c("data.table", "RColorBrewer", "ggplot2", "ComplexHeatmap", "circlize", "igraph", "multiplot", "gtools", "RNA", "ggpubr")
sapply(inlibs, library, character.only = TRUE)
setOldClass("igraph")
setClassUnion("FactorOrNULL", c("factor", "NULL"))
setClass(
  Class = "loop",
  representation = representation(g = "igraph", loop = "data.table", split = "FactorOrNULL"),
  validity=function(object) {
    if (nrow(object@loop) < 1) {
      return("No loop was given.")
    }
    if (any(!colnames(object@loop) %in% c("loop", "PromGene", "rowid"))) {
      return("Column names of the hash slot do not match to 'loop', 'PromGene', 'rowid'.")
    }
    if (length(E(object@g)) < 1) {
      return("No edge in the g slot.")
    }
    if (length(V(object@g)) < 1) {
      return("No vertices in the g slot.")
    }
    if (any(!E(object@g)$loop %in% object@loop[["loop"]])) {
      return("The loop attribute of the edges of the g slot does not match to the loop column in the loop slot.")
    }
    if (length(E(object@g)) != length(unique(object@loop[["loop"]])) {
      return("The number of edges must equal to the number of unique loops in loop slot.")
    }
    if (!is.null(object@split) && any(is.na(object@split))) {
      return("The split slot can not contain NA.")
    }
    if (!is.null(object@split) && length(object@split) != nrow(object@loop)) {
      return("The split slot does not have the same number as the row number of the loop slot.")
    }
    return(TRUE)
  }
)

# ginfo includes info for all genes.
setClassUnion("datatableOrNULL", c("data.table", "NULL"))
setClass(
  Class = "info",
  representation = representation(gene = "data.table", tad = "datatableOrNULL"),
  validity=function(object) {
    if (nrow(object@gene) < 1) {
      return("No gene information was given.")
    }
    if (any(!c("chr", "start", "end", "gene", "type") %in% colnames(object@gene))) {
      return("Column names of the gene slot do not contain all of 'chr', 'start', 'end', 'gene', 'type'.")
    }
    if (any(!c("chr", "start", "end") %in% key(object@gene))) {
      return("The gene slot must be keyed on 'chr', 'start' and 'end' columns.")
    }
    if (!is.null(object@tad) && nrow(object@tad) < 1) {
      return("No tad information was given.")
    }
    if (!is.null(object@tad) && any(!colnames(object@tad) %in% c("chr", "start", "end", "tadid"))) {
      return("Column names of the hash slot do not match to 'chr', 'start', 'end' and 'tadid'.")
    }
    if (!is.null(object@tad) && any(!c("chr", "start", "end") %in% key(object@tad))) {
      return("The tad slot must be keyed on 'chr', 'start' and 'end' columns.")
    }
    return(TRUE)
  }
)

setClass(
  Class = "fet",
  representation = representation(dat_list = "list", hash = "data.table"),
  validity=function(object) {
    if (length(object@dat_list) < 1) {
      return("No feture was given.")
    }
    for (i in seq_along(object@dat_list)) {
      if (!all(sapply(object@dat_list[[i]], class) %in% c("numeric", "integer"))) {
        return(paste0("Data in ", i, "th matrix does not have all columns as numeric or integer."))
      }
    }
    if (is.null(names(object@dat_list))) {
      return("dat_list slot does not names.")
    }
    if (nrow(object@hash) < 1) {
      return("No hash was given.")
    }
    if (nrow(object@hash) != length(object@dat_list)) {
      return("The row number of the hash slot does not equal to the length of dat_list.")
    }
    if (any(colnames(object@hash) != c("sms", "grps", "enhs", "grps_enhs"))) {
      return("Column names of the hash slot do not match to 'sms', 'grps', 'enhs', 'grps_enhs'.")
    }
    if (any(sapply(object@hash, class) != "factor")) {
      return("Columns of the hash slot of the object class must be of 'factor' class.")
    }
    return(TRUE)
  }
)
