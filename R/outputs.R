#'@export 
as_array <- function(fit, ...) UseMethod("as_array")
#'@export 
as_array.dbtmb <- function(fit, name='median') {
  rp <- fit$obj$report()
  o <- array(rp[[name]], rp$rdims)
  dimnames(o) <- with(fit$meta, list(age=age, yob=yob, cc=cc_id$ISO_A3))
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
predict.dbtmb <- function(fit, par, type=c("eversex"), cutoff=24, long=TRUE, quick=FALSE,...) {
  if (missing(par)) par <- fit$obj$env$last.par
  r <- fit$obj$report(par) # can check this and not redoing env$reportenv$...
  scale <- array(r$lambda, r$rdims)
  shape <- sweep(array(0, r$rdims), 3, r$alpha_vec, '+')
  skew <- sweep(array(0, r$rdims), 3, r$a_vec, '+')
  o <- sapply(cutoff, function(x) mapply(F_gllogisI, x, scale, shape, skew))
  dim(o) <- c(r$rdims, length(cutoff))
  dimnames(o) <- with(fit$meta, list(age=age, yob=yob, cc=cc_id$ISO_A3, kut=cutoff))
  if (quick) 
    return(o)
  if (long) {
    o <- as.data.table(o)
    setnames(o, colnames(o), c('age', 'yob', 'ISO_A3', 'kut', type))
    ccols <- c('age', 'yob', 'kut')
    o[, ccols] <- lapply(o[, ccols], as.double)
  }
  o
}
