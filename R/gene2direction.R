#' @include GRFLoopClass.R

#' @export gene2direction
setGeneric(name = "gene2direction",
  def = function(gene_list, info.obj){
    standardGeneric("gene2direction")
  }
)

#' @rdname gene2direction-methods
setMethod(f = "gene2direction",
  signature = c("list", "info"),
  definition = function(gene_list, info.obj){
    idx <- sapply(gene_list, function(x)is.character(x))
    message(sum(idx), " intervals overlap with genes")
    message(sum(!idx), " intervals do not overlap with genes")
    gene_list <- gene_list[idx]
    # deg labels for these genes
    col_nm <- copy(colnames(info.obj@gene)[grep("^DEG", colnames(info.obj@gene))])
    deg_list <- lapply(gene_list, function(gs, info.obj, col_nm){
      idx <- chmatch(gs, info.obj@gene[["gene"]])
      stopifnot(all(!is.na(idx)))
      deg_dat <- info.obj@gene[idx, col_nm, with = FALSE]
      stopifnot(nrow(deg_dat) == length(gs))
      return(deg_dat)
    }, info.obj = info.obj, col_nm = col_nm)
    up_list <- lapply(deg_list, function(dd){
      data.table(direction = "up", dd[, lapply(.SD, function(x){
        up <- 100 * sum(tolower(x) == "up", na.rm = TRUE)/length(x);
        dn <- 100 * sum(tolower(x) == "down", na.rm = TRUE)/length(x);
        return(up)
      }), .SDcols = colnames(dd)])
    })
    dn_list <- lapply(deg_list, function(dd){
      data.table(direction = "dn", dd[, lapply(.SD, function(x){
        dn <- 100 * sum(tolower(x) == "down", na.rm = TRUE)/length(x);
        return(dn)
      }), .SDcols = colnames(dd)])
    })
    stopifnot(length(up_list) == length(dn_list))
    deg_pct_list <- lapply(seq_along(up_list), function(j){
      rbindlist(list(up_list[[j]], dn_list[[j]]))
    })
    return(deg_pct_list)
  }
)
