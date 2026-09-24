## ---------------------------------------------------------------------------
## Fit diagnostics
##
## Two families of quantities are kept apart on purpose:
##
##  * "sampler" diagnostics describe how well the *estimation scheme* did its
##    job (mixing for MCMC, quality of the Gaussian approximation for Laplace,
##    convergence of the optimiser for MAP). They are by-products of fitting,
##    are cheap, and are therefore computed once in slgp() / .retrain_SLGP()
##    and stored in the 'diagnostics' slot.
##
##  * "predictive" diagnostics describe how well the *fitted field* describes
##    data (log predictive density, PIT calibration, and for Laplace the
##    importance-sampling efficiency of the Gaussian draws). They require
##    evaluating the model at every observation and are therefore computed on
##    request by summary(..., diagnostics = TRUE).
## ---------------------------------------------------------------------------

## Sampler diagnostics for a stanfit produced by rstan::sampling().
#' @keywords internal
#' @noRd
.slgp_mcmc_diag <- function(fit) {
  s <- rstan::summary(fit)$summary
  keep <- setdiff(rownames(s), "lp__")
  rhat <- s[keep, "Rhat"]
  ess  <- s[keep, "n_eff"]
  ndiv <- tryCatch(rstan::get_num_divergent(fit), error = function(e) NA_integer_)
  ntd  <- tryCatch(rstan::get_num_max_treedepth(fit), error = function(e) NA_integer_)
  bfmi <- tryCatch(rstan::get_bfmi(fit), error = function(e) NA_real_)
  list(
    scheme      = "MCMC",
    n_chains    = ncol(as.array(fit)[, , 1, drop = FALSE]),
    n_iter      = nrow(as.array(fit)[, , 1, drop = FALSE]),
    rhat_range  = range(rhat, na.rm = TRUE),
    rhat_max    = max(rhat, na.rm = TRUE),
    n_rhat_bad  = sum(rhat > 1.01, na.rm = TRUE),
    ess_range   = range(ess, na.rm = TRUE),
    ess_min     = min(ess, na.rm = TRUE),
    n_ess_bad   = sum(ess < 400, na.rm = TRUE),
    n_divergent = ndiv,
    n_max_treedepth = ntd,
    bfmi_min    = if (all(is.na(bfmi))) NA_real_ else min(bfmi, na.rm = TRUE),
    lp_mean     = unname(s["lp__", "mean"])
  )
}

## Sampler diagnostics for the Laplace approximation. 'logq' is the Gaussian
## log-density of each stored draw under the approximation; keeping it (a
## vector of length ndraws) is what later allows summary() to form
## self-normalised importance weights against the true posterior without
## storing the p x p covariance matrix.
#' @keywords internal
#' @noRd
.slgp_laplace_diag <- function(hessian, nugget, logq, ndraws) {
  ev <- tryCatch(eigen(hessian, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NA_real_)
  list(
    scheme      = "Laplace",
    ndraws      = ndraws,
    nugget      = nugget,
    hessian_min_eigen = if (all(is.na(ev))) NA_real_ else min(ev),
    hessian_cond      = if (all(is.na(ev))) NA_real_ else max(abs(ev)) / min(abs(ev)),
    logq        = logq
  )
}

## Optimiser diagnostics for MAP.
#' @keywords internal
#' @noRd
.slgp_map_diag <- function(fit) {
  list(
    scheme      = "MAP",
    return_code = tryCatch(fit$return_code, error = function(e) NA_integer_),
    converged   = isTRUE(tryCatch(fit$return_code == 0L, error = function(e) NA))
  )
}

## ---------------------------------------------------------------------------
## Predictive diagnostics
## ---------------------------------------------------------------------------

## Pointwise log predictive density: an n x ndraws matrix whose (i, s) entry is
## log f(t_i | x_i, epsilon_s). This is the object every information criterion
## below is built from.
#' @keywords internal
#' @noRd
.slgp_loglik <- function(object, data) {
  cols <- c(object@responseName, object@covariateName)
  nd <- data[, cols, drop = FALSE]
  pred <- predict(object, newdata = nd, type = "density")
  pdf_cols <- grep("^pdf_", colnames(pred), value = TRUE)
  if (!length(pdf_cols))
    stop("No density columns returned; cannot compute predictive diagnostics.")
  ll <- log(as.matrix(pred[, pdf_cols, drop = FALSE]))
  ll[!is.finite(ll)] <- NA_real_
  ll
}

## log(mean(exp(x))) computed stably.
#' @keywords internal
#' @noRd
.log_mean_exp <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  m <- max(x)
  m + log(mean(exp(x - m)))
}

## Probability integral transform of the observations under the posterior
## predictive CDF, with a Kolmogorov-Smirnov statistic against the uniform.
#' @keywords internal
#' @noRd
.slgp_pit <- function(object, data) {
  cols <- c(object@responseName, object@covariateName)
  nd <- data[, cols, drop = FALSE]
  pred <- predict(object, newdata = nd, type = "cdf")
  cdf_cols <- grep("^cdf_", colnames(pred), value = TRUE)
  if (!length(cdf_cols)) return(NULL)
  u <- rowMeans(as.matrix(pred[, cdf_cols, drop = FALSE]), na.rm = TRUE)
  u <- u[is.finite(u)]
  if (length(u) < 2L) return(NULL)
  n <- length(u)
  us <- sort(u)
  ## Two-sided KS statistic against U(0, 1).
  d <- max(pmax(seq_len(n) / n - us, us - (seq_len(n) - 1) / n))
  list(pit = u, ks = d,
       ks_pvalue = tryCatch(stats::ks.test(u, "punif", 0, 1)$p.value,
                            error = function(e) NA_real_),
       mean = mean(u), sd = stats::sd(u))
}

## Self-normalised importance-sampling efficiency of the Laplace approximation:
## how many of the 'ndraws' Gaussian draws are worth, in effective terms, when
## reweighted towards the true posterior. Low values mean the Gaussian is a poor
## match and that Laplace-based uncertainty should be treated with caution.
#' @keywords internal
#' @noRd
.slgp_laplace_is_ess <- function(object, ll) {
  d <- object@diagnostics
  if (is.null(d$logq) || ncol(ll) != length(d$logq)) return(NULL)
  sigma2 <- object@hyperparams$sigma2
  eps <- object@coefficients
  ## log prior: epsilon ~ N(0, sigma2 I)
  logprior <- -0.5 * rowSums(eps^2) / sigma2 -
    0.5 * ncol(eps) * log(2 * pi * sigma2)
  logpost <- colSums(ll, na.rm = TRUE) + logprior
  lw <- logpost - d$logq
  lw <- lw - max(lw[is.finite(lw)])
  w <- exp(lw)
  w[!is.finite(w)] <- 0
  if (sum(w) <= 0) return(NULL)
  w <- w / sum(w)
  list(ess = 1 / sum(w^2), ndraws = length(w))
}

## Wall-clock cost of a fit, split into the two phases that scale differently:
## 'setup' (normalisation, quadrature pre-computation, basis evaluation) grows
## with the number of distinct covariate values and the grid resolution, while
## 'estimation' grows with the rank and, for MCMC, the number of iterations.
## All times are wall-clock seconds: with cores > 1 the cumulated CPU time
## exceeds the elapsed time, and child-process times are not reported at all
## on Windows, so only 'elapsed' is comparable across platforms.
## For MCMC the per-chain warmup / sampling split reported by Stan is kept too.
#' @keywords internal
#' @noRd
.slgp_timing <- function(setup, estimation, fit = NULL, method = NA_character_) {
  out <- list(
    setup      = unname(setup["elapsed"]),
    estimation = unname(estimation["elapsed"]),
    total      = unname(setup["elapsed"] + estimation["elapsed"])
  )
  if (identical(method, "MCMC") && !is.null(fit)) {
    ct <- tryCatch(rstan::get_elapsed_time(fit), error = function(e) NULL)
    if (!is.null(ct)) {
      out$warmup  <- unname(ct[, "warmup"])
      out$sample  <- unname(ct[, "sample"])
      out$chain_total <- unname(rowSums(ct))
    }
  }
  out
}
