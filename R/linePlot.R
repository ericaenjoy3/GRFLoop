#' @include GRFLoopClass.R GRFLoopGeneric.R

#' @export linePlot
setGeneric(name = "linePlot",
  def = function(loop.obj, info.obj, pdffout, vtype, uniqueLoopGene = FALSE){
    standardGeneric("linePlot")
  }
)

#' @rdname linePlot-methods
setMethod(f = "linePlot",
  signature = c("loop", "info"),
  definition = function(loop.obj, info.obj, pdffout, vtype, uniqueLoopGene) {

    proc <- function(dd) {
      nd <- dd[, list(
        .SD[,.N],
        sapply(.SD, function(x){mean(x)}),
        sapply(.SD, function(x){mean(x) + sd(x)/sqrt(length(x))}),
        sapply(.SD, function(x){mean(x) - sd(x)/sqrt(length(x))})), by = c("grp", "variable")]
      setnames(nd, paste0("V", 1:4), c("N", "mean", "upside", "dnside"))
      nd[, c("grp") := paste(grp, "(", N, ")")]
      return(nd)
    }

    con_num <- 1:3

    dat_list <- lapply(con_num, function(j){
      message(j);
      unlist(geneListProc(loop.obj, info.obj, nmin = j, nmax = j, type = vtype, uniqueLoopGene = uniqueLoopGene))
      })

    rna_list <- lapply(dat_list, function(gid){
        dd <- info.obj@gene[chmatch(gid, gene), grep("tpm", colnames(info.obj@gene)), with = FALSE]
        return(dd)
      })

    rna_dat <- rbindlist(lapply(seq_along(rna_list), function(j){
       data.table(grp = con_num[j], rna_list[[j]])
      }))
    col_nm <- colnames(rna_dat)[grep("tpm", colnames(rna_dat))]
    ncol_nm <- gsub("tpm_", "", col_nm)
    setnames(rna_dat, col_nm, ncol_nm)
    setcolorder(rna_dat, c("grp", "MEF", "D3", "D6", "D9", "ESC"))
    ncol_nm <- colnames(rna_dat)[-1]
    rna_dat[, (ncol_nm) := lapply(.SD, function(x)log2(x + 0.05)), .SD = ncol_nm]
    mat <- combn(ncol_nm, 2)
    diff_dat <- copy(rna_dat[, .(grp)])
    for (j in 1:ncol(mat)) {
      sms <- ncol_nm[ncol_nm %in% mat[, j]]
      psm <- paste(sms, collapse="_")
      set(diff_dat, i = 1:nrow(rna_dat), j = psm, value = rna_dat[[sms[2]]] - rna_dat[[sms[1]]])
    }

    rna_ldat <- melt(rna_dat, id.vars = "grp")
    rna_ldat[, variable := factor(variable, levels = c("MEF", "D3", "D6", "D9", "ESC"))]
    # z-score
    # rna_ldat[, value := (value - mean(value))/sd(value), by = .(grp)]
    plot_dat <- proc(rna_ldat)

    # ggline(rna_ldat, x = "variable", y = "value", color = "variable", palette = "jco", add = "median_iqr", facet.by = "grp")
    theme_set(theme_grey(base_size=15))

    # separate line plots per connectivity
    p1 <- ggplot(plot_dat, aes(x = variable, y = mean, group = 1)) + 
      geom_line(aes(color = factor(grp))) +
      facet_grid(.~grp) +
      geom_ribbon(data = plot_dat, aes(ymin = dnside, ymax = upside), fill = "grey70", alpha = 0.4) +
      labs(x = "", y = expression(paste(log[2],"(TPM)"))) +
      theme(legend.title = element_blank(), panel.spacing = unit(2, "lines"), legend.position = "top")

    # plot all lines of all connectivities in one
    p2 <- ggplot(plot_dat, aes(x = variable, y = mean, group = grp)) + 
      geom_line(aes(color = factor(grp))) +
      geom_ribbon(data = plot_dat, aes(ymin = dnside, ymax = upside), fill = "grey70", alpha = 0.4) +
      labs(x = "", y = expression(paste(log[2],"(TPM)"))) +
      theme(legend.title = element_blank(), panel.spacing = unit(2, "lines"), legend.position = "top")

    pdf(pdffout)
    grid.arrange(p1, p2, ncol = 2)
    dev.off()

    # parwise comparison
    diff_ldat <- melt(diff_dat, id.vars = "grp")

    p3 <- ggboxplot(diff_ldat, x = "variable", y = "value", color = "variable", palette = "jco", xlab = "", 
      ylab = expression(paste(log[2],"(FC)")), legend.title = "")
    p3 <- facet(p3, scales = "free", facet.by = c("grp")) + rotate_x_text(45)

    pdf(gsub(".pdf", "_FC.pdf", pdffout), width = 7*1.5)
    grid.arrange(p3, ncol = 1)
    dev.off()
  }
)

