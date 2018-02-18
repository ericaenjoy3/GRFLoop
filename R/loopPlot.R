#' @include GRFLoopClass.R

#' @export loopTypePlot
setGeneric(name = "loopTypePlot",
  def = function(loop.obj, pdffout){
    standardGeneric("loopTypePlot")
  }
)

#' @rdname loopTypePlot-methods
setMethod(f = "loopTypePlot",
  signature = c("loop"),
  definition = function(loop.obj, pdffout) {
  	dat <- data.table(etype = E(loop.obj@g)$etype)
  	stats <- dat[, .(N = .N, Pct = round(100*.N/nrow(dat), digits = 2)), by = etype]
  	theme_set(theme_grey(base_size=15))
	  p1 <- ggplot(stats, aes(x = etype, y = N, fill = etype))+
		geom_bar(stat = "identity") +
		labs(x = "",y = "Counts") +
		theme(legend.title = element_blank(),panel.spacing = unit(2, "lines"), legend.position = "top", 
			axis.text.x = element_text(angle = 45, hjust = 1)) +
		guides(fill = guide_legend(reverse = TRUE), colour = guide_legend(reverse = TRUE))+
		geom_text(aes(x = etype, y = N, label = Pct), hjust = 0.5, vjust = -0.5,size = 4, data = stats, inherit.aes = FALSE)
	ggsave(pdffout, p1)
  }
)

#' @export loopTypePlot
setGeneric(name = "loopDistPlot",
  def = function(loop.obj, pdffout){
    standardGeneric("loopDistPlot")
  }
)

#' @rdname loopTypePlot-methods
setMethod(f = "loopDistPlot",
  signature = c("loop"),
  definition = function(loop.obj, pdffout) {
    if (!is.null(E(loop.obj@g)$cluster)) {
      dat <- data.table(etype = E(loop.obj@g)$etype, dist = E(loop.obj@g)$dist, cluster = E(loop.obj@g)$cluster)
      dat[, etype := factor(etype, levels = sort(unique(etype)), ordered = TRUE)]
      cmp <- data.table(combn(unique(E(loop.obj@g)$cluster), 2))
      p1 <- ggviolin(dat, x = "cluster", y = "dist", fill = "cluster", 
        add = "boxplot", add.params = list(fill = "white"),
        xlab = "", ylab = "",
        legend.title = "") +
        stat_compare_means(comparisons = cmp)
      p1 <- facet(p1, scales = "free", facet.by = c("etype"), ncol = 1, panel.labs.background = list(fill = "transparent"))
      ggsave(filename = pdffout, p1)
    } else {
      dat <- data.table(etype = E(loop.obj@g)$etype, dist = E(loop.obj@g)$dist)
      dat[, etype := factor(etype, levels = c("Enh|Enh", "Prom|Enh", "Prom|Prom"), ordered = TRUE)]
      cmp <- data.table(combn(unique(E(loop.obj@g)$etype), 2))
      p1 <- ggviolin(dat, x = "etype", y = "dist", fill = "etype", 
        add = "boxplot", add.params = list(fill = "white"),
        xlab = "", ylab = "",
        legend.title = "") +
        stat_compare_means(comparisons = cmp)
      ggsave(filename = pdffout, p1)      
    }
  }
)