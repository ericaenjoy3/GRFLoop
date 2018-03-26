#' @include GRFLoopClass.R GRFLoopGeneric.R

#' @export maxCoreg
setGeneric(name = "selCoreg",
  def = function(loop.obj, info.obj){
    standardGeneric("selCoreg")
  }
)

#' @rdname geneList-methods
setMethod(f = "selCoreg",
  signature = c("loop", "info"),
  definition = function(loop.obj, info.obj) {

    type <- "Enh"

    # dedup loop in the loop slot of loop.obj
    kpt.idx <- copy(!duplicated(loop.obj@loop[["loop"]]))
    loop_hash <- copy(loop.obj@loop[kpt.idx])
    setkeyv(loop_hash, "loop")

    # identify incident loop
    stopifnot(sum(V(loop.obj@g)$vtype == type) > 0)
    ve <- V(loop.obj@g)$name[V(loop.obj@g)$vtype == type]
    ed <- incident_edges(loop.obj@g, ve)
    kpt_num <- sapply(ed, length) >= 2

    # (genunine) extract PromGene from loop slot of loop.obj
    gene_list <- list()
    gene_list[[1]] <- lapply(ed[which(kpt_num)], function(es, loop_hash){
      # potentially multiple loops => mutiple fragments
      lp <- es$loop
      dd <- loop_hash[loop %in% lp]
      stopifnot(all(dd[, all(!is.na(gene1)), by = loop][, V1]))
      dd <- dd[, unique(gene1), by = loop]
      gs <- split(dd[, V1], dd[, loop])
      # keep gene names the same order as loop names
      gs <- gs[lp]
      names(gs) <- NULL
      return(gs)
    }, loop_hash = loop.obj@loop)
    names(gene_list[[1]]) <- NULL
    shadow_gen_list <- lapply(gene_list[[1]], function(glist)sapply(glist, function(vec)mixedsort(vec)))

    # (genunine) unique gene sets across connectomes
    kpt_agset <- !duplicated(shadow_gen_list)

    # (genuine) remove duplicated gene sets within connectomes
    kpt_wgset <- sapply(gene_list[[1]], function(v_list){
      nv_list <- unique(v_list);
      if (length(nv_list) == 1 ) {
        return(FALSE)
      } else {
        return(TRUE)
      }
    })

    stopifnot(length(kpt_agset) == length(kpt_wgset))
    message(sum(!(kpt_agset & kpt_wgset)), " hubs removed due to duplication in idnetical genes contacted within or across hubs.")

    # (genunine) final gene_list
    gene_list[[1]] <- gene_list[[1]][which(kpt_agset & kpt_wgset)]
    ve <- ve[which(kpt_num)][which(kpt_agset & kpt_wgset)]
    ed <- ed[which(kpt_num)][which(kpt_agset & kpt_wgset)]

    # (genuine) gene pairs
    gpair_list <- list()
    gpair_list[[1]] <- lapply(gene_list[[1]], function(glist) {
      gpair <- glist2gpair(glist)$gpair_dat
      return(gpair)
    })

    gidx_list <- list()
    gidx_list[[1]] <- lapply(gene_list[[1]], function(glist){
      gidx <- glist2gpair(glist)$gidx_dat
      return(gidx)
    })

    deg_list <- list()
    deg_list[[1]] <- lapply(gpair_list[[1]], function(dd){
      gene2pairwiseLab(dd, info.obj)
    })

    # keep co-regulation edges
    ed_list <- lapply(seq(deg_list[[1]]), function(i){
      
      idx <- which(deg_list[[1]][[i]] == "2")
      if (length(idx) == 0) {
        return(NA)
      }

      # promoter fragments
      dd <- gidx_list[[1]][[i]]
      dd <- separate(dd, "Var1", c("prom_idx1", "gene_idx1"), 
          convert = TRUE, remove = TRUE) %>% separate("Var2", c("prom_idx2", "gene_idx2"), 
          convert = TRUE, remove = TRUE);
      dd <- dd[idx]
      ed_vec <- { 
        ed_idx <- unique(c(dd[, prom_idx1], dd[, prom_idx2]))
        stopifnot(length(ed_idx) > 0)
        ed[[i]][ed_idx]
      }
      return(ed_vec)
    })
    kpt_coreg <- !is.na(ed_list)
    eid <- ed_list[kpt_coreg][[1]]
    for(i in 2:length(ed_list[kpt_coreg])) {
      eid <- c(eid, ed_list[kpt_coreg][[i]])
    }
    
    loop.obj@g <- subgraph.edges(loop.obj@g, eid = eid, delete.vertices = TRUE)
    loop.obj@loop <- loop.obj@loop[loop %in% E(loop.obj@g)$loop]
    validObject(loop.obj)
    return(list(loop.obj = loop.obj))

  }
)