#' @include GRFLoopClass.R

#' @export gene2pairwiseLab
setGeneric(name = "gene2pairwiseLab",
  def = function(gene_list, info.obj){
    standardGeneric("gene2pairwiseLab")
  }
)

#' @rdname gene2pairwiseLab-methods
setMethod(f = "gene2pairwiseLab",
  signature = c("list", "info"),
  definition = function(gene_list, info.obj){
    idx <- sapply(gene_list, function(x)is.character(x))
    message(sum(idx), " intervals overlap with genes")
    message(sum(!idx), " intervals do not overlap with genes")
    gene_list <- gene_list[idx]
    # deg labels for these genes
    col_nm <- colnames(info.obj@gene)[grep("^DEG", colnames(info.obj@gene))]
    deg_list <- lapply(gene_list, function(gs, info.obj, col_nm){
      idx <- chmatch(gs, info.obj@gene[["gene"]])
      stopifnot(all(!is.na(idx)))
      deg_dat <- info.obj@gene[idx, col_nm, with = FALSE]
      stopifnot(nrow(deg_dat) == length(gs))
      return(deg_dat)
    }, info.obj = info.obj, col_nm = col_nm)
    lab_list <- lapply(deg_list, function(dd){
      row_idmat <- combn(nrow(dd), 2)
      sapply(1:ncol(row_idmat), function(j){
        vec <- dd[row_idmat[,j]][[ncol(dd)]]
        if (any(is.na(vec))) {
          return(0)
        } else if (all(vec == "Up")) {
          return(1)
        } else if (all(vec == "Down")) {
          return(-1)
        } else {
          return(0)
        }
      }) 
    })
    lab <- unlist(lab_list)
    return(lab)
  }
)