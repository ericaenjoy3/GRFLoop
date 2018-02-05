#' @include GRFLoopClass.R GRFLoopGeneric.R

#' @export linePlot
setGeneric(name = "linePlot",
  def = function(loop.obj, info.obj, fet.obj, linepdf, uniqueLoopGene = FALSE){
    standardGeneric("linePlot")
  }
)

#' @rdname linePlot-methods
setMethod(f = "linePlot",
  signature = c("loop", "info"),
  definition = function(loop.obj, info.obj, fet.obj, linepdf, uniqueLoopGene) {

    con_num <- 1:3
    dat_list <- lapply(con_num, function(j)unlist(geneListProc(loop.obj, info.obj, nmin = j, nmax = j, type = "Prom", uniqueLoopGene = uniqueLoopGene)))

    row_dat <- loop.obj@loop[, .SD[1, c("rowid"), with = FALSE], by = PromGene]
    stopifnot(all(sapply(dat_list, function(vec)all(vec %in% row_dat[, PromGene]))))

    idx <- grep("RNA", fet.obj@hash[["sms"]])[1]
    rna_dat <- fet.obj@dat_list[[idx]]

    rna_list <- lapply(dat_list, function(gid){
        ridx <- row_dat[chmatch(gid, row_dat[,PromGene]), rowid]
        dd <- rna_dat[ridx]
        return(dd)
      })

    rna_dat <- rbindlist(lapply(seq_along(rna_list), function(j){
       data.table(grp = con_num[j], rna_list[[j]])
      }))
    rna_ldat <- melt(rna_dat, id.vars = "grp")

    plot_dat <- rna_ldat[, list(
      .SD[,.N],
      sapply(.SD, function(x){mean(x)}),
      sapply(.SD, function(x){mean(x) + sd(x)/sqrt(length(x))}),
      sapply(.SD, function(x){mean(x) - sd(x)/sqrt(length(x))})), by = c("grp", "variable")]
    setnames(plot_dat, paste0("V", 1:4), c("N", "mean", "upside", "dnside"))
    plot_dat[, c("grp") := paste(grp, "(", N, ")")]

    # ggline(rna_ldat, x = "variable", y = "value", color = "variable", palette = "jco", add = "median_iqr", facet.by = "grp")
    theme_set(theme_grey(base_size=15))
    p1 <- ggplot(plot_dat, aes(x = variable, y = mean, group = 1)) + 
    geom_line(aes(color = factor(grp))) +
    facet_grid(.~grp) +
    geom_ribbon(data = plot_dat, aes(ymin = dnside, ymax = upside), fill = "grey70", alpha = 0.4) +
    labs(x = "", y = "Z scores") +
    coord_cartesian(ylim = c(-1, 1)) + 
    theme(legend.title = element_blank(), panel.spacing = unit(2, "lines"), legend.position = "top")

    ggsave(linepdf, p1)
  }
)

