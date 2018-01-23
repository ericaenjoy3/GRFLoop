#' @include GRFLoopClass.R GRFLoopGeneric.R
setMethod(f = "multiHeat",
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

#' @rdname quantRm-methods
setMethod(f = "medianViolin",
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
    theme_set(theme_grey(base_size=15))
    p1 <- ggplot(ndat, aes(x = variable, y = value, fill = split)) +
    geom_violin(position = position_dodge(width = 0.9)) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.9)) +
    facet_grid(sms ~ grps, scale = "free") +
    labs(x = "", y = "") +
    theme(legend.title = element_blank(), panel.spacing = unit(2, "lines"),
      legend.position = "top",
      axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(filename = pdffout, p1)
  }
)

setMethod(f = "conPlot",
  signature = c("loop", "fet"),
  definition = function(loop.obj, fet.obj, dout){
    # dedup loop in the loop slot of loop.obj
    kpt.idx <- !duplicated(loop.obj@loop[["loop"]])
    loop_hash <- loop.obj@loop[kpt.idx]
    setkeyv(loop_hash, "loop")
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
      pdffout <- file.path(dout, fet.obj@hash[["enhs"]][1],
        paste(fet.obj@hash[["enhs"]][1], type[j], "conBarPlot.pdf", sep = "_"))
      theme_set(theme_grey(base_size=15))
      p1 <- ggplot(tabulate, aes(x = variable, y = value)) +
      geom_bar(stat = "identity", fill = "#FF9999") +
      labs(x = "", y = "Counts") +
      theme(legend.position = "none")
      ggsave(filename = pdffout, p1)
      # loop through index
      stopifnot(type[j] %in% fet.obj@hash[["grps"]])
      fouts <- copy(fet.obj@hash)
      fouts[, "pdffout" := file.path(dout, enhs, paste(enhs, grps, sms, "vbox.pdf", sep = "_"))]
      sm_idx_vec <- which(fet.obj@hash[["grps"]] == type[j])
      stopifnot(length(sm_idx_vec) > 0)
      for (sm_idx in sm_idx_vec) {
        message("processing", fouts[["pdffout"]][sm_idx])
        dat <- rbindlist(lapply(unique(sort(centric.cnt)), function(k, fet.obj, rowid_list, centric.cnt) {
          idx_vec <- unlist(rowid_list[centric.cnt == k])
          stopifnot(max(idx_vec) <= nrow(fet.obj@dat_list[[sm_idx]]))
          con <- ifelse(max(centric.cnt) == k, " + ", " ")
          con <- paste0(k, con, "(", length(rowid_list[centric.cnt == k]), ")")
          dd <- data.table(con = rep(con, length(idx_vec)), fet.obj@dat_list[[sm_idx]][idx_vec])
          return(dd)
        }, fet.obj = fet.obj, rowid_list =rowid_list, centric.cnt = centric.cnt))
        ldat <- melt(dat, id.vars = "con")
        ldat[["con"]] <- factor(ldat[["con"]], levels = unique(ldat[["con"]]))
        ndat <- subset(ldat, grepl("ESC", ldat[["variable"]]), ignore.case = TRUE)
        if (grepl("RNA", fet.obj@hash[["sms"]][sm_idx], ignore.case = TRUE)) {
          ylab <- expression(paste(log[2], (TPM)))
        } else if (all(sapply(ndat[["value"]], function(num)num %%1 == 0))) {
          ylab <- "Peak Counts"
        } else {
          ylab <- "RPKM"
        }
        pdffout <- fouts[["pdffout"]][sm_idx]
        theme_set(theme_grey(base_size=15))
        p1 <- ggplot(ndat, aes(x = con, y = value, fill = con)) +
          geom_violin(position = position_dodge(width = 0.9)) +
          geom_boxplot(width = 0.1, position = position_dodge(width = 0.9)) +
          labs(x = "", y = ylab) +
          theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 1))
        ggsave(pdffout, p1)
      }
    }
  }
)
