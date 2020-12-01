#'@export 
as_array <- function(fit, ...) UseMethod("as_array")
#'@export 
as_array.dbtmb <- function(fit, name='median') {
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
#' @param long long form
#' @param quick array form
#' @export
predict.dbtmb <- function(fit, par, type=c("eversex"), cutoff=24, long=TRUE, quick=FALSE,...) {
  if (missing(par)) par <- fit$obj$env$last.par
  options(tape.parallel=FALSE, DLL='dbtmb')
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
    o$age <- as.double(o$age)
    o$yob <- as.double(o$yob)
    o$kut <- as.double(o$kut)
  }
  o
}
#'@export
eversex_ui <- function(fit, smp, cutoff=15:20, n_cores=20) {
  message('predicting sample:\n')
  o <- parallel::mclapply(1:nrow(smp), function(x) predict(fit, smp[x, ], cutoff=cutoff, quick=TRUE), mc.cores = n_cores)
  o <- do.call('rbind', o)
  o <- parallel::mclapply(1:ncol(o), function(x) ktools::quantile95(o[, x]), mc.cores=n_cores)
  o <- do.call('cbind', o)
  dim(o) <- c(3, fit$obj$env$reportenv$rdims, length(cutoff))
  dimnames(o) <- with(fit$meta, list(c('lo', 'med', 'up'),  
    age=age, yob=yob, ISO_A3=cc_id$ISO_A3, kut=cutoff))
  o <- data.table::as.data.table(o) 
  o$age <- as.double(o$age)
  o$yob <- as.double(o$yob)
  o$kut <- as.double(o$kut)
  data.table::dcast(o, ... ~ V5, value.var='value')
}

#'@export
uncertainty <- function(obj, ...) UseMethod("uncertainty")
#'@export
uncertainty.dbtmb <- function(fit, smp, n_cores=20) {
  options(tape.parallel=FALSE, DLL='dbtmb')
  r <- fit$obj$report(smp[1, ]) # do you know why this is prefered

  tmp  <- parallel::mclapply(1:nrow(smp), function(x) {
    r <- fit$obj$report(smp[x, ])
    age_e <- r$ccxage + rep(r$age_rw2, r$rdims[3])
    yob_e <- r$ccxyob + rep(r$yob_rw2, r$rdims[3])
    list(itc=r$intercept, skew=r$a_vec, shape=r$alpha_vec, age=age_e, yob=yob_e, cc=r$cc_vec, median=r$median)
  }, mc.cores=n_cores)

  itc <- quantile95(sapply(tmp, function(x) x$itc))

  skew <- apply(sapply(tmp, function(x) x$skew), 1, quantile95)
  dimnames(skew) <- with(fit$meta, list(c('lo', 'med', 'up'), ISO_A3=cc_id$ISO_A3))
  skew <- data.table::as.data.table(t(skew), TRUE)
  setnames(skew, 'rn', 'ISO_A3')

  shape <- apply(sapply(tmp, function(x) x$shape), 1, quantile95)
  dimnames(shape) <- with(fit$meta, list(c('lo', 'med', 'up'), ISO_A3=cc_id$ISO_A3))
  shape <- data.table::as.data.table(t(shape), TRUE)
  setnames(shape, 'rn', 'ISO_A3')

  cc_e <- apply(sapply(tmp, function(x) x$cc), 1, quantile95)
  dimnames(cc_e) <- with(fit$meta, list(c('lo', 'med', 'up'), ISO_A3=cc_id$ISO_A3))
  cc_e <- data.table::as.data.table(t(cc_e), TRUE)
  setnames(cc_e, 'rn', 'ISO_A3')

  median_e <- sapply(tmp, function(x) x$median)
  median_e <- parallel::mclapply(1:nrow(median_e), function(x) quantile95(median_e[x, ]), mc.cores=n_cores)
  median_e <- do.call('cbind', median_e)
  dim(median_e) <- c(3, fit$obj$env$reportenv$rdims)
  dimnames(median_e) <- with(fit$meta, list(c('lo', 'med', 'up'),  
    age=age, yob=yob, ISO_A3=cc_id$ISO_A3))
  median_e <- data.table::as.data.table(median_e) 
  scol <- c('age', 'yob')
  median_e$age <- as.double(median_e$age)
  median_e$yob_e <- as.double(median_e$yob_e)

  age_e <- apply(sapply(tmp, function(x) x$age), 1, quantile95)
  dim(age_e) <- c(3, r$rdims[-2])
  dimnames(age_e) <- with(fit$meta, list(c('lo', 'med', 'up'),  
    age=age, ISO_A3=cc_id$ISO_A3))
  age_e <- data.table::as.data.table(age_e)
  age_e[, age := as.double(age)]

  yob_e <- apply(sapply(tmp, function(x) x$yob), 1, quantile95)
  dim(yob_e) <- c(3, r$rdims[-1])
  dimnames(yob_e) <- with(fit$meta, list(c('lo', 'med', 'up'),  
    yob=yob, ISO_A3=cc_id$ISO_A3))
  yob_e <- data.table::as.data.table(yob_e)
  yob_e[, yob := as.double(yob)]
  
  list(
    itc = itc, 
    skew = skew,
    shape = shape,
    cc  = cc_e,
    age = data.table::dcast(age_e, ... ~ V3, value.var='value'),
    yob = data.table::dcast(yob_e, ... ~ V3, value.var='value'),
    median = data.table::dcast(median_e, ... ~ V4, value.var='value')
  )
}