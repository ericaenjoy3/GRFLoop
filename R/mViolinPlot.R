#' @include GRFLoopClass.R

#' @export mViolinViolin
setGeneric(name = "mViolinPlot",
  def = function(loop.obj, fet.obj, pdffout){
    standardGeneric("mViolinPlot")
  }
)

#' @rdname mViolinPlot-methods
setMethod(f = "mViolinPlot",
  signature = c("loop", "fet"),
  definition = function(loop.obj, fet.obj, pdffout){
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
    p1 <- ggviolin(ndat, x = "split", y = "value", fill = "split", 
      add = "boxplot", add.params = list(fill = "white"), 
      facet.by = c("sms", "grps"), xlab = "", ylab = "", legend.title = "")+
    stat_compare_means(comparisons = cmp)
    p1 <- facet(p1, scales = "free", facet.by = c("sms", "grps"))
    ggsave(filename = pdffout, p1)
  }
)
