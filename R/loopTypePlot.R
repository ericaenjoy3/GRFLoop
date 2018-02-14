#' @include GRFLoopClass.R

#' @export loopTypePlot
setGeneric(name = "loopTypePlot",
  def = function(loop.obj, pdffout){
    standardGeneric("loopTypePlot")
  }
)

#' @rdname infoFilter-methods
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