#' @include GRFLoopClass.R

#' @export conPlot
setGeneric(name = "conPlot",
  def = function(loop.obj, dout){
    standardGeneric("conPlot")
  }
)

#' @rdname conPlot-methods
setMethod(f = "conPlot",
  signature = c("loop"),
  definition = function(loop.obj, dout){
    dir.create(dout, showWarnings = FALSE, recursive = TRUE)
    # loop through promoter/enhancer vertices
    type <- c("Prom", "Enh")
    for (j in seq_along(type)) {
      # identify incident loop
      ve <- V(loop.obj@g)$name[V(loop.obj@g)$type == type[j]]
      ed <- incident_edges(loop.obj@g, ve)
      # extract rowid from loop slot of loop.obj
      rowid_list <- lapply(ed, function(es, loop_hash){
        lp <- gsub("|", "_", as_ids(es), fixed = TRUE)
        rowid <- loop_hash[lp][["rowid"]]
        stopifnot(all(!is.na(rowid)))
        return(rowid)
      }, loop_hash = loop_hash)
      # glob index for indexing different connectivites
      centric.cnt <- sapply(rowid_list, length)
      tabulate <- as.data.frame(table(centric.cnt))
      tabulate[, 1] <- as.integer(tabulate[, 1])
      cumsum.idx <- nrow(tabulate) + 1 - which.max(cumsum(rev(tabulate[,2]))>10)
      ceilling.idx <- max(which(tabulate[,2] > 10)) + 1
      centric.cnt <- replace(centric.cnt,
        which(centric.cnt > tabulate[min(cumsum.idx, ceilling.idx), 1]),
        tabulate[min(cumsum.idx, ceilling.idx), 1])
      tabulate <- as.data.frame(table(centric.cnt))
      colnames(tabulate) <- c("variable", "value")
      levels(tabulate[["variable"]])[nlevels(tabulate[["variable"]])] <- paste0(levels(tabulate[["variable"]])[nlevels(tabulate[["variable"]])], "+")
      # boxplot of connectivites
      pdffout <- file.path(dout, 
        paste(type[j], "conBarPlot.pdf", sep = "_"))
      theme_set(theme_grey(base_size = 15))
      p1 <- ggplot(tabulate, aes(x = variable, y = value)) +
      geom_bar(stat = "identity", fill = "#FF9999") +
      labs(x = "", y = "Counts") +
      theme(legend.position = "none")
      ggsave(filename = pdffout, p1)
    }
  }
)