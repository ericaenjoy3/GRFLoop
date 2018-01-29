#' @include GRFLoopClass.R

#' @export ceilFet
setGeneric(name = "ceilFet",
  def = function(fet.obj){
    standardGeneric("ceilFet")
  }
)

#' @rdname ceilFet-methods
setMethod(f = "ceilFet",
  signature = c("fet"),
  definition = function(fet.obj) {
    # update dat_list slot of fet object
    idx_list <- split(seq_len(nrow(fet.obj@hash)), fet.obj@hash[["sms"]])
    rng_list <- lapply(seq_along(idx_list), function(i, idx_list, fet.obj){
      idx_vec <- idx_list[[i]]
      if (!grepl("H3K27AC", names(idx_list)[i], ignore.case = TRUE)) {
        vec <- range(unlist(fet.obj@dat_list[idx_vec]))
      } else {
        vec <- quantile(unlist(fet.obj@dat_list[idx_vec]), probs = c(0, 0.95))
        names(vec) <- NULL
      }
      invisible(sapply(idx_vec, function(idx){
        col_nm <- colnames(fet.obj@dat_list[[idx]])
        fet.obj@dat_list[[idx]][, (col_nm) := lapply(.SD, function(x){x[x > vec[2] ] <- vec[2]; return(x)}), .SDcols = col_nm]}))
      return(rbind(vec))}, idx_list = idx_list, fet.obj = fet.obj)
    names(rng_list) <- names(idx_list)
    validObject(fet.obj)
    return(list(fet.obj = fet.obj, rng_list = rng_list))
  }
)
