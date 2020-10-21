#' Sample TMB fit
#'
#' @param fit The TMB fit
#' @param nsample Number of samples
#' @param random_only Random only
#' @param verbose If TRUE prints additional information.
#'
#' @return Sampled fit.
#' @export
sample_tmb <- function(fit, nsample = 1000, random_only = TRUE, verbose = TRUE) {
  to_tape <- TMB:::isNullPointer(fit$obj$env$ADFun$ptr)
  if (to_tape) 
    fit$obj$retape(FALSE)
  par.full <- fit$obj$env$last.par.best  
  if(!random_only) {
    if(verbose) print("Calculating joint precision")
    hess <- TMB::sdreport(fit$obj, fit$fit$par, getJointPrecision=TRUE)
    if(verbose) print("Drawing sample")
    smp <- rmvnorm_sparseprec(nsample, par.full, hess)
  } else {
    r <- fit$obj$env$random
    par_f <- par.full[-r]
    par_r <- par.full[r]
    hess_r <- fit$obj$env$spHess(par.full, random = TRUE)
    smp_r <- rmvnorm_sparseprec(nsample, par_r, hess_r)
    smp <- matrix(0, nsample, length(par.full))
    smp[ , r] <- smp_r
    smp[ ,-r] <- matrix(par_f, nsample, length(par_f), byrow = TRUE)
    colnames(smp) <- names(par.full)
  }
  smp
}

rmvnorm_sparseprec <- function(n, mean = rep(0, nrow(prec)), prec = diag(lenth(mean))) {
  z = matrix(rnorm(n * length(mean)), ncol = n)
  # y=z
  # y = backsolve(chol(prec), z) 
  # L_inv = Matrix::Cholesky(prec, perm=FALSE, LDL=FALSE)
  L_inv = Matrix::Cholesky(prec)
  v <- mean + 
    Matrix::solve(as(L_inv, "pMatrix"), 
      Matrix::solve(
        Matrix::t(as(L_inv, "Matrix")), z) )
  as.matrix(Matrix::t(v))
}