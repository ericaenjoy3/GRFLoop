#' @include GRFLoopClass.R GRFLoopGeneric.R

#' @export linePlot
setGeneric(name = "linePlot",
  def = function(loop.obj, info.obj, linepdf, uniqueLoopGene = FALSE){
    standardGeneric("linePlot")
  }
)

#' @rdname linePlot-methods
setMethod(f = "linePlot",
  signature = c("loop", "info"),
  definition = function(loop.obj, info.obj, linepdf, uniqueLoopGene) {

    con_num <- 1:3

    enh_list <- lapply(con_num, function(j){
      message(j);
      unlist(geneListProc(loop.obj, info.obj, nmin = j, nmax = j, type = "Enh", uniqueLoopGene = uniqueLoopGene))
      })

    rna_list <- lapply(dat_list, function(gid){
        dd <- info.obj@gene[chmatch(gid, gene), grep("tpm", colnames(info.obj@gene)), with = FALSE]
        return(dd)
      })

    rna_dat <- rbindlist(lapply(seq_along(rna_list), function(j){
       data.table(grp = con_num[j], rna_list[[j]])
      }))
    col_nm <- colnames(rna_dat)[grep("tpm", colnames(rna_dat))]
    setnames(rna_dat, colnames(rna_dat), gsub("tpm_", "", colnames(rna_dat)))
    rna_ldat <- melt(rna_dat, id.vars = "grp")
    rna_ldat[, variable := factor(variable, levels = c("MEF", "D3", "D6", "D9", "ESC"))]
    # log2
    rna_ldat[, ("value") := lapply(.SD, function(x)log2(x + 0.05)), .SD = "value"]
    # z-score
    # rna_ldat[, value := (value - mean(value))/sd(value), by = .(grp, variable)]

    plot_dat <- rna_ldat[, list(
      .SD[,.N],
      sapply(.SD, function(x){mean(x)}),
      sapply(.SD, function(x){mean(x) + sd(x)/sqrt(length(x))}),
      sapply(.SD, function(x){mean(x) - sd(x)/sqrt(length(x))})), by = c("grp", "variable")]
    setnames(plot_dat, paste0("V", 1:4), c("N", "mean", "upside", "dnside"))
    plot_dat[, c("grp") := paste(grp, "(", N, ")")]

    # ggline(rna_ldat, x = "variable", y = "value", color = "variable", palette = "jco", add = "median_iqr", facet.by = "grp")
    theme_set(theme_grey(base_size=15))

    # separate line plots per connectivity
    p1 <- ggplot(plot_dat, aes(x = variable, y = mean, group = 1)) + 
      geom_line(aes(color = factor(grp))) +
      facet_grid(.~grp) +
      geom_ribbon(data = plot_dat, aes(ymin = dnside, ymax = upside), fill = "grey70", alpha = 0.4) +
      labs(x = "", y = "log2(TPM)") +
      theme(legend.title = element_blank(), panel.spacing = unit(2, "lines"), legend.position = "top")

    # plot all lines of all connectivities in one
    p2 <- ggplot(plot_dat, aes(x = variable, y = mean, group = grp)) + 
      geom_line(aes(color = factor(grp))) +
      geom_ribbon(data = plot_dat, aes(ymin = dnside, ymax = upside), fill = "grey70", alpha = 0.4) +
      labs(x = "", y = "log2(TPM)") +
      theme(legend.title = element_blank(), panel.spacing = unit(2, "lines"), legend.position = "top")

    # 
    pdf(pdffout)
    grid.arrange(p1, p2, ncol=2)
    dev.off()
  }
)

