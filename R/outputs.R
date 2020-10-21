#'@export 
as_array <- function(fit, ...) UseMethod("as_array")
as_array.dbtmb <- function(fit, name='median') {
  rp <- fit$obj$report()
  o <- array(rp[[name]], rp$rdims)
  dimnames(o)[[1]] <- sort(unique(all_cb$age))
  dimnames(o)[[2]] <- sort(unique(all_cb$yob))
  dimnames(o)[[3]] <- rownames(readRDS('data/nbmat_SSA.rds'))
  o
}
#'@export 
predict <- function(fit, ...) UseMethod("predict")
#' Predict from fitted dbtmb
#' 
#' Predict from fitted dbtmb
#' 
#' @param fit fitted dbtmb class
#' @param par parameter set to simulate
#' @param cutoff vector, can be length one
#' @export
predict.dbtmb <- function(fit, par, type=c("eversex"), cutoff=24, long=TRUE,...) {
  if (missing(par)) par <- fit$obj$env$last.par
  r <- fit$obj$report(par) # can check this and not redoing env$reportenv$...
  scale <- array(r$lambda, r$rdims)
  shape <- sweep(array(0, r$rdims), 3, r$alpha_vec, '+')
  skew <- sweep(array(0, r$rdims), 3, r$a_vec, '+')
  o <- sapply(cutoff, function(x) mapply(F_gllogisI, x, scale, shape, skew))
  dim(o) <- c(r$rdims, length(cutoff))
  dimnames(o)[[1]] <- sort(unique(all_cb$age))
  dimnames(o)[[2]] <- sort(unique(all_cb$yob))
  dimnames(o)[[3]] <- rownames(readRDS('data/nbmat_SSA.rds'))
  dimnames(o)[[4]] <- cutoff
  if (long) {
    o <- as.data.table(o)
    setnames(o, colnames(o), c('age', 'yob', 'ISO_A3', 'kut', type))
    ccols <- c('age', 'yob', 'kut')
    o[, ccols] <- lapply(o[, ccols], as.double)
  }
  o
}
