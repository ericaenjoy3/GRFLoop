#' @include GRFLoopClass.R

#' @export multiHeatPlot
setGeneric(name = "multiHeatPlot",
  def = function(loop.obj, fet.obj, rng_list, pdffout, scheme = c("wr", "br")){
    standardGeneric("multiHeatPlot")
  }
)

#' @rdname multiHeatPlot-methods
setMethod(f = "multiHeatPlot",
  signature = c("loop", "fet"),
  definition = function(loop.obj, fet.obj, rng_list, pdffout, scheme = c("wr", "br")) {
    scheme <- match.arg(scheme)
    if (is.null(loop.obj@split)) {
      ht_list <- NULL
    } else {
      col1 <- if(length(levels(loop.obj@split)) < 3) {
        brewer.pal(3, "Dark2")[1:2]
      } else if (length(levels(loop.obj@split)) <= 8) {
        brewer.pal(length(levels(loop.obj@split)), "Dark2")
      } else if (length(levels(loop.obj@split)) <= 12) {
        brewer.pal(length(levels(loop.obj@split)), "Paired")
      } else if (length(levels(loop.obj@split)) > 12) {
        colorRampPalette(brewer.pal(12, "Paired"))(length(levels(loop.obj@split)))
      }
      ht_list <- Heatmap(as.character(loop.obj@split), col = structure(col1, names = levels(loop.obj@split)), name = "Groups",
        show_row_names = FALSE, show_column_names = FALSE, width = unit(3, "mm"),
        use_raster = FALSE, show_row_dend = FALSE, show_column_dend = FALSE, split = loop.obj@split,
        combined_name_fun = NULL, gap = unit(5, "mm"))
    }
    for (i in seq_along(fet.obj@dat_list)) {
      message("processing ",i, " in heatMulti")
      rng <- rng_list[[as.numeric(fet.obj@hash$sms)[i]]]
      colvec <- c("white", "red")
      message("evaluating data range")
      if (fet.obj@hash$sms[i] == "RNA" & scheme == "br") {
        rng <- c(rng, 0)
        rng <- sort(rng)
        colvec <- c("blue", "white", "red")
      }
      message("heatmap")
      ht_list <- ht_list + Heatmap(fet.obj@dat_list[[i]], colorRamp2(rng, colvec),
        show_row_names = FALSE, show_column_names = TRUE, cluster_rows = FALSE,
        show_row_dend = FALSE,  cluster_columns = FALSE, show_column_dend = FALSE,
        split = loop.obj@split, gap = unit(5, "mm"),
        heatmap_legend_param = list(title = fet.obj@hash[["sms"]][i], color_bar = "continuous",
        legend_direction = "horizontal", legend_width = unit(5, "cm"), title_position = "lefttop"))
    }
    pdf(pdffout, width = 7 * length(fet.obj@dat_list)/2)
    draw(ht_list, annotation_legend_side = "bottom")
    dev.off()
  }
)
