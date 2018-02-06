#' @include GRFLoopClass.R

#' @export mViolinViolin
setGeneric(name = "mViolinPlot",
  def = function(loop.obj, fet.obj, dout){
    standardGeneric("mViolinPlot")
  }
)

#' @rdname mViolinPlot-methods
setMethod(f = "mViolinPlot",
  signature = c("loop", "fet"),
  definition = function(loop.obj, fet.obj, dout){
    dir.create(dout, showWarnings = FALSE, recursive = TRUE)
    kpt.idx <- !duplicated(loop.obj@loop[["loop"]])
    # fetures
    dat <- rbindlist(
      lapply(seq_along(fet.obj@dat_list), function(i, fet.obj, kpt.idx){
        message(i)
        dd <- melt(fet.obj@dat_list[[i]][kpt.idx])
        nd <- data.table(sms = rep(fet.obj@hash[["sms"]][i], nrow(dd)),
          grps = rep(fet.obj@hash[["grps"]][i], nrow(dd)),
          split = loop.obj@split[kpt.idx],
          dd);
        return(nd)
      }, fet.obj = fet.obj, kpt.idx = kpt.idx))
    string <- gsub("DAY", "D", levels(dat[["variable"]]))
    lel <- unique(gsub(".*(MEF|D3|D6|D9|ESC).*", "\\1", string))
    dat[, c("variable") := gsub("DAY", "D", variable)]
    dat[, c("variable") := factor(gsub(".*(MEF|D3|D6|D9|ESC).*", "\\1", variable), levels = lel)]
    ndat <- subset(dat, dat[["variable"]] == "ESC")
    # plot
    # theme_set(theme_grey(base_size=15))
    # p1 <- ggplot(ndat, aes(x = variable, y = value, fill = split)) +
    # geom_violin(position = position_dodge(width = 0.9)) +
    # geom_boxplot(width = 0.1, position = position_dodge(width = 0.9)) +
    # facet_grid(sms ~ grps, scale = "free") +
    # labs(x = "", y = "") +
    # theme(legend.title = element_blank(), panel.spacing = unit(2, "lines"),
    #   legend.position = "top",
    #   axis.text.x = element_text(angle = 90, hjust = 1))
    # alternative plot
    cmp <- data.table(combn(unique(ndat[["split"]]), 2))
    ndat_list <- split(ndat, ndat[["sms"]])
    plist <- lapply(seq_along(ndat_list), function(j, ndat_list, dout) {
      pdffout <- file.path(dout, paste0(names(dat_list)[j], "_medianViolin.pdf"))
      p1 <- ggviolin(ndat_list[[j]], x = "split", y = "value", fill = "split", 
      add = "boxplot", add.params = list(fill = "white"), 
      facet.by = c("sms", "grps"), xlab = "", ylab = "", 
      legend.title = "", title = names(ndat_list)[j]) +
      stat_compare_means(comparisons = cmp)
      p1 <- facet(p1, scales = "free", facet.by = c("grps"))
      ggsave(filename = pdffout, p1)
      return(p1)
    }, ndat_list = ndat_list, dout = dout)   
  }
)
