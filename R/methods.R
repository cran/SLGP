## =============================================================================
##  S4 methods for the SLGP class: print, show, summary, plot, predict,
##  simulate, update, and a handful of standard extractors.
## =============================================================================


## ----------------------------------------------------------------------------
## Internal helpers (not exported)
## ----------------------------------------------------------------------------

## Number of posterior/prior coefficient draws stored in the model.
.slgp_ndraws <- function(object) {
  if (length(object@coefficients) == 0L) return(0L)
  nrow(object@coefficients)
}

.slgp_discrete <- function(object) {
  if (methods::.hasSlot(object, "discrete") && length(object@discrete))
    isTRUE(object@discrete) else FALSE
}

## Compact wall-clock formatting: seconds below a minute, then m/s, then h/m.
.slgp_fmt_time <- function(x) {
  if (is.null(x) || !is.finite(x)) return("<unknown>")
  if (x < 60) return(paste0(format(x, digits = 3), "s"))
  if (x < 3600) return(paste0(floor(x / 60), "m ", round(x %% 60), "s"))
  paste0(floor(x / 3600), "h ", round((x %% 3600) / 60), "m")
}

.slgp_nIntegral <- function(object) {
  if (methods::.hasSlot(object, "nIntegral") && length(object@nIntegral))
    object@nIntegral else 101
}
## Warn when the caller overrides the response type or the support grid
.slgp_check_discrete <- function(object, discrete, nIntegral) {
  fitted_discrete <- .slgp_discrete(object)
  if (isTRUE(fitted_discrete) && !isTRUE(discrete))
    warning("Model was fitted with discrete = TRUE but discrete = FALSE was ",
            "requested. The response is treated as continuous.")
  if (!isTRUE(fitted_discrete) && isTRUE(discrete))
    warning("Model was fitted with discrete = FALSE but discrete = TRUE was ",
            "requested. The response is treated as discrete.")
  if (isTRUE(discrete) && !missing(nIntegral) &&
      length(nIntegral) == 1L && !is.na(nIntegral) &&
      nIntegral != .slgp_nIntegral(object))
    warning("'nIntegral' (", nIntegral, ") differs from the value used at ",
            "fitting (", .slgp_nIntegral(object), "). The support grid will ",
            "not match the one the model was estimated on.")
  invisible(NULL)
}

## Human-readable label for the basis family.
.slgp_basis_label <- function(object) {
  bf <- object@basisFunctionsUsed
  if (length(bf) == 0L || is.na(bf)) return("<unset>")
  bf
}

## Pretty one-line range, e.g. "[0, 50]".
.slgp_fmt_range <- function(x) {
  if (length(x) < 2L) return("<unset>")
  paste0("[", format(x[1], trim = TRUE), ", ", format(x[2], trim = TRUE), "]")
}

## Is the model fitted (as opposed to a pure prior, method = "none")?
.slgp_is_fitted <- function(object) {
  length(object@method) == 1L && !is.na(object@method) && object@method != "none"
}


## ----------------------------------------------------------------------------
## print / show
## ----------------------------------------------------------------------------

#' Print a brief description of an SLGP model
#'
#' Compact, one-screen summary of a fitted (or prior) \code{\link{SLGP-class}}
#' object: the model formula, the estimation method, the basis family and its
#' rank, the response range, the data dimensions, and the key hyperparameters.
#' No internal slots are dumped.
#'
#' @param x An object of class \code{\link{SLGP-class}}.
#' @param ... Ignored, for S3/S4 generic compatibility.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rep(seq(0, 1, length.out = 6), each = 5))
#' d$y <- rnorm(nrow(d), mean = sin(2 * pi * d$x), sd = 0.2)
#'
#' fit <- slgp(y ~ x,
#'   data = d, method = "MAP", basisFunctionsUsed = "RFF",
#'   predictorsLower = 0, predictorsUpper = 1,
#'   responseRange = range(d$y), seed = 1,
#'   opts_BasisFun = list(nFreq = 20, MatParam = 5 / 2))
#'
#' fit
#'
#' @export
setMethod("print", signature(x = "SLGP"), function(x, ...) {
  method <- if (length(x@method)) x@method else "<unset>"
  ndraws <- .slgp_ndraws(x)
  nobs   <- if (nrow(x@data)) nrow(x@data) else 0L
  ndistinct <- if (length(x@covariateName) && nobs > 0L) {
    nrow(unique(x@data[, x@covariateName, drop = FALSE]))
  } else NA_integer_

  sigma2 <- tryCatch(x@hyperparams$sigma2, error = function(e) NULL)
  ls     <- tryCatch(x@hyperparams$lengthscale, error = function(e) NULL)

  cat("Spatial Logistic Gaussian Process (SLGP) model\n")
  cat("  Formula        : ", deparse(x@formula), "\n", sep = "")
  cat("  Response       : ", x@responseName,
      "  in ", .slgp_fmt_range(x@responseRange),
      if (.slgp_discrete(x))
        paste0("  (discrete, ", .slgp_nIntegral(x), " support points)") else "",
      "\n", sep = "")
  cat("  Covariate(s)   : ", paste(x@covariateName, collapse = ", "), "\n", sep = "")
  cat("  Estimation     : ", method,
      if (!.slgp_is_fitted(x)) "  (no coefficients)" else "", "\n", sep = "")
  cat("  Basis          : ", .slgp_basis_label(x),
      "  (rank p = ", x@p, ")\n", sep = "")
  if (!is.null(ls))
    cat("  Length-scales  : ", paste(format(ls, digits = 3), collapse = ", "),
        "\n", sep = "")
  if (!is.null(sigma2))
    cat("  Variance sigma2: ", format(sigma2, digits = 4), "\n", sep = "")
  cat("  Data           : ", nobs, " obs",
      if (!is.na(ndistinct)) paste0(", ", ndistinct, " distinct covariate value(s)") else "",
      "\n", sep = "")
  .tm <- if (methods::.hasSlot(x, "diagnostics")) x@diagnostics$timing else NULL
  if (!is.null(.tm))
    cat("  Fitting time   : ", .slgp_fmt_time(.tm$total), "\n", sep = "")
  if (ndraws > 0L)
    cat("  Coefficient draws: ", ndraws, "\n", sep = "")
  if (length(x@logPost) == 1L && is.finite(x@logPost))
    cat("  Log-posterior  : ", format(x@logPost, digits = 6), "\n", sep = "")
  invisible(x)
})

#' @rdname print-SLGP-method
#' @param object An object of class \code{\link{SLGP-class}}.
#' @export
setMethod("show", signature(object = "SLGP"), function(object) {
  print(object)
})


## ----------------------------------------------------------------------------
## summary
## ----------------------------------------------------------------------------

#' Summarise a fitted SLGP model
#'
#' Extends \code{\link{print}} with diagnostics. Two kinds are reported.
#'
#' \emph{Sampler diagnostics} describe how well the estimation scheme  worked.
#' This is the by-products of fitting, stored in the model, and always shown:
#' R-hat, effective sample sizes, divergent transitions and BFMI for
#' \code{"MCMC"}; the nugget added to the Hessian and its conditioning for
#' \code{"Laplace"}; the optimiser return code for \code{"MAP"}.
#'
#' \emph{Predictive diagnostics} describe how well the fitted field describes
#' data: the mean log predictive density per observation, and the probability integral
#' transform (PIT) of the observations with a Kolmogorov-Smirnov statistic
#' against the uniform. These require evaluating the model at every observation
#' and are therefore computed only when \code{diagnostics = TRUE}.
#'
#' For a Laplace fit, the importance-sampling effective sample size of the
#' Gaussian draws relative to the true posterior is also reported: it says how
#' many of the stored draws are worth, in effective terms, once reweighted, and
#' hence reflects whether the Gaussian approximation is adequate.
#'
#' @param object An object of class \code{\link{SLGP-class}}.
#' @param diagnostics Logical; if \code{TRUE}, compute the predictive
#'   diagnostics described above. Default \code{FALSE}, since the cost grows
#'   with the number of observations and of coefficient draws.
#' @param newdata Optional \code{data.frame} of held-out observations on which
#'   to compute the predictive diagnostics. If \code{NULL} (default), the
#'   training data are used, in which case the log predictive density and the
#'   PIT are in-sample and therefore optimistic; WAIC and PSIS-LOO correct for
#'   this, the raw log predictive density does not.
#' @param ... Ignored.
#'
#' @return An object of class \code{"summary.SLGP"} (a list), returned invisibly.
#'
#' @export
setMethod("summary", signature(object = "SLGP"),
          function(object, diagnostics = FALSE, newdata = NULL, ...) {
            ndraws <- .slgp_ndraws(object)
            nobs   <- if (nrow(object@data)) nrow(object@data) else 0L
            ndistinct <- if (length(object@covariateName) && nobs > 0L) {
              nrow(unique(object@data[, object@covariateName, drop = FALSE]))
            } else NA_integer_

            out <- list(
              formula            = object@formula,
              responseName       = object@responseName,
              responseRange      = object@responseRange,
              discrete           = .slgp_discrete(object),
              nIntegral          = .slgp_nIntegral(object),
              covariateName      = object@covariateName,
              method             = if (length(object@method)) object@method else NA_character_,
              fitted             = .slgp_is_fitted(object),
              basisFunctionsUsed = .slgp_basis_label(object),
              p                  = object@p,
              hyperparams        = object@hyperparams,
              nobs               = nobs,
              ndistinct          = ndistinct,
              ndraws             = ndraws,
              logPost            = if (length(object@logPost) == 1L) object@logPost else NA_real_,
              timing             = if (methods::.hasSlot(object, "diagnostics"))
                object@diagnostics$timing else NULL,
              sampler            = if (methods::.hasSlot(object, "diagnostics"))
                object@diagnostics else list(),
              predictive         = NULL
            )

            if (isTRUE(diagnostics) && .slgp_is_fitted(object)) {
              dat <- if (is.null(newdata)) object@data else newdata
              ll  <- .slgp_loglik(object, dat)
              pit <- .slgp_pit(object, dat)
              out$predictive <- list(
                n        = nrow(ll),
                insample = is.null(newdata),
                elpd     = mean(apply(ll, 1L, .log_mean_exp), na.rm = TRUE),
                pit      = pit,
                is_ess   = if (identical(out$sampler$scheme, "Laplace"))
                  .slgp_laplace_is_ess(object, ll) else NULL
              )
            }

            class(out) <- "summary.SLGP"
            out
          })

#' Print method for SLGP summaries
#'
#' @param x An object of class \code{"summary.SLGP"}.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @exportS3Method base::print
#'
print.summary.SLGP <- function(x, ...) {
  cat("Summary of a Spatial Logistic Gaussian Process (SLGP) model\n")
  cat("  Formula        : ", deparse(x$formula), "\n", sep = "")
  cat("  Response       : ", x$responseName,
      "  in ", .slgp_fmt_range(x$responseRange),
      if (isTRUE(x$discrete))
        paste0("  (discrete, ", x$nIntegral, " support points)") else "",
      "\n", sep = "")
  cat("  Covariate(s)   : ", paste(x$covariateName, collapse = ", "), "\n", sep = "")
  cat("  Estimation     : ", x$method,
      if (!isTRUE(x$fitted)) "  (no coefficients)" else "", "\n", sep = "")
  cat("  Basis          : ", x$basisFunctionsUsed,
      "  (rank p = ", x$p, ")\n", sep = "")
  if (!is.null(x$hyperparams$lengthscale))
    cat("  Length-scales  : ",
        paste(format(x$hyperparams$lengthscale, digits = 3), collapse = ", "),
        "\n", sep = "")
  if (!is.null(x$hyperparams$sigma2))
    cat("  Variance sigma2: ", format(x$hyperparams$sigma2, digits = 4), "\n", sep = "")
  cat("  Data           : ", x$nobs, " obs",
      if (!is.na(x$ndistinct)) paste0(", ", x$ndistinct, " distinct covariate value(s)") else "",
      "\n", sep = "")

  ## -- sampler diagnostics ---------------------------------------------------
  d <- x$sampler
  cat("\nEstimation diagnostics\n")
  if (!is.null(x$timing)) {
    cat("  Fitting time      : ", .slgp_fmt_time(x$timing$total),
        "   (setup ", .slgp_fmt_time(x$timing$setup),
        ", estimation ", .slgp_fmt_time(x$timing$estimation), ")\n", sep = "")
    if (!is.null(x$timing$warmup))
      cat("  Per chain         : warmup ",
          paste(vapply(x$timing$warmup, .slgp_fmt_time, ""), collapse = ", "),
          " | sampling ",
          paste(vapply(x$timing$sample, .slgp_fmt_time, ""), collapse = ", "),
          "\n", sep = "")
  }
  if (x$ndraws > 0L)
    cat("  Coefficient draws : ", x$ndraws, "\n", sep = "")
  if (!is.na(x$logPost) && is.finite(x$logPost)) {
    cat("  Log-posterior     : ", format(x$logPost, digits = 6),
        if (identical(d$scheme, "MCMC")) "  (posterior mean)"
        else " (at mode for MAP/Laplace)", "\n", sep = "")
  } else {
    cat("  Log-posterior     : not available for method '", x$method, "'\n", sep = "")
  }
  if (identical(d$scheme, "MCMC")) {
    cat("  Chains / iter     : ", d$n_chains, " / ", d$n_iter, "\n", sep = "")
    cat("  R-hat             : max ", format(d$rhat_max, digits = 4),
        "  (", d$n_rhat_bad, " of ", x$p, " above 1.01)\n", sep = "")
    cat("  Eff. sample size  : min ", format(d$ess_min, digits = 4),
        "  (", d$n_ess_bad, " of ", x$p, " below 400)\n", sep = "")
    cat("  Divergences       : ", d$n_divergent,
        if (!is.na(d$n_max_treedepth))
          paste0("   max-treedepth hits: ", d$n_max_treedepth) else "",
        "\n", sep = "")
    if (!is.na(d$bfmi_min))
      cat("  Min BFMI          : ", format(d$bfmi_min, digits = 3), "\n", sep = "")
    if (isTRUE(d$n_rhat_bad > 0) || isTRUE(d$n_divergent > 0))
      cat("  ! Chains show signs of poor mixing; consider more iterations.\n")
  }
  if (identical(d$scheme, "Laplace")) {
    cat("  Hessian nugget    : ",
        if (isTRUE(d$nugget > 0)) format(d$nugget, digits = 3) else "none (invertible)",
        "\n", sep = "")
    if (!is.na(d$hessian_cond))
      cat("  Hessian condition : ", format(d$hessian_cond, digits = 4), "\n", sep = "")
  }
  if (identical(d$scheme, "MAP") && !is.na(d$return_code))
    cat("  Optimiser         : ",
        if (isTRUE(d$converged)) "converged" else
          paste0("return code ", d$return_code, " (not converged)"), "\n", sep = "")

  ## -- predictive diagnostics ------------------------------------------------
  pd <- x$predictive
  if (!is.null(pd)) {
    cat("\nPredictive diagnostics (",
        if (isTRUE(pd$insample)) "in-sample" else "held-out",
        ", n = ", pd$n, ")\n", sep = "")
    cat("  Mean log density  : ", format(pd$elpd, digits = 5), " per observation\n", sep = "")
    if (!is.null(pd$pit))
      cat("  PIT vs uniform    : KS = ", format(pd$pit$ks, digits = 4),
          if (!is.na(pd$pit$ks_pvalue))
            paste0(", p = ", format(pd$pit$ks_pvalue, digits = 3)) else "",
          "\n", sep = "")
    if (!is.null(pd$is_ess))
      cat("  Laplace IS-ESS    : ", format(pd$is_ess$ess, digits = 4),
          " of ", pd$is_ess$ndraws, " draws\n", sep = "")
    if (isTRUE(pd$insample))
      cat("  (in-sample; pass newdata= for an honest assessment)\n")
  }
  invisible(x)
}


## ----------------------------------------------------------------------------
## Standard extractors
## ----------------------------------------------------------------------------

#' Extract the coefficient draws of an SLGP model
#'
#' @param object An object of class \code{\link{SLGP-class}}.
#' @param ... Ignored.
#' @return The matrix of finite-rank GP coefficients (draws in rows, basis
#'   functions in columns).
#' @export
setMethod("coef", signature(object = "SLGP"), function(object, ...) {
  object@coefficients
})

#' Wall-clock cost of fitting an SLGP model
#'
#' Returns the time spent fitting, split into the two phases that scale
#' differently: \code{setup} (normalisation, quadrature pre-computation and
#' basis evaluation, which grow with the number of distinct covariate values and
#' with \code{nDiscret}) and \code{estimation} (the call to \pkg{rstan}, which
#' grows with the rank and, for MCMC, with the number of iterations).
#'
#' The result is a one-row \code{data.frame}, so that timings collected over a
#' set of fits can be stacked with \code{\link[base]{rbind}} and plotted
#' directly.
#'
#' All times are wall-clock seconds. With more than one core the cumulated CPU
#' time exceeds the elapsed time, so only elapsed times are comparable across
#' settings and platforms.
#'
#' @param object An object of class \code{\link{SLGP-class}}.
#' @param ... Ignored.
#'
#' @return A \code{data.frame} with one row and the columns \code{method},
#'   \code{p}, \code{nobs}, \code{ndraws}, \code{setup}, \code{estimation},
#'   \code{total} and, for MCMC, \code{n_chains}, \code{warmup} and
#'   \code{sample} (means over chains). Returns \code{NULL} for models fitted
#'   with a version of the package that did not record timings.
#'
#' @examples
#' \dontrun{
#' ## Compare estimation schemes over several replicates
#' fits <- lapply(1:10, function(i)
#'   slgp(depth ~ long, data = quakes, method = "MAP",
#'        basisFunctionsUsed = "RFF", seed = i,
#'        opts_BasisFun = list(nFreq = 200, MatParam = 5/2)))
#' tm <- do.call(rbind, lapply(fits, timing))
#' boxplot(total ~ method, data = tm, ylab = "Fitting time [s]")
#' }
#'
#' @export
setGeneric("timing", function(object, ...) standardGeneric("timing"))

#' @rdname timing
#' @export
setMethod("timing", signature(object = "SLGP"), function(object, ...) {
  tm <- if (methods::.hasSlot(object, "diagnostics"))
    object@diagnostics$timing else NULL
  if (is.null(tm)) return(NULL)
  out <- data.frame(
    method     = if (length(object@method)) object@method else NA_character_,
    p          = object@p,
    nobs       = nrow(object@data),
    ndraws     = .slgp_ndraws(object),
    setup      = tm$setup,
    estimation = tm$estimation,
    total      = tm$total,
    stringsAsFactors = FALSE
  )
  out$n_chains <- if (is.null(tm$warmup)) NA_integer_ else length(tm$warmup)
  out$warmup   <- if (is.null(tm$warmup)) NA_real_    else mean(tm$warmup)
  out$sample   <- if (is.null(tm$sample)) NA_real_    else mean(tm$sample)
  out
})

#' Number of observations used to fit an SLGP model
#'
#' @param object An object of class \code{\link{SLGP-class}}.
#' @param ... Ignored.
#' @return Integer number of observations.
#' @export
setMethod("nobs", signature(object = "SLGP"), function(object, ...) {
  nrow(object@data)
})

#' Model formula of an SLGP object
#'
#' @param x An object of class \code{\link{SLGP-class}}.
#' @param ... Ignored.
#' @return The model \code{formula}.
#' @export
setMethod("formula", signature(x = "SLGP"), function(x, ...) {
  x@formula
})


## ----------------------------------------------------------------------------
## predict
## ----------------------------------------------------------------------------

#' Predict from a fitted SLGP model
#'
#' Single entry point for all distributional summaries produced by an
#' \code{\link{SLGP-class}} model. The \code{type} argument selects what is
#' returned and dispatches to the corresponding workhorse function:
#' \describe{
#'   \item{\code{"density"}}{conditional density values, \code{newdata} must
#'   contain both the
#'     response and the covariate columns.}
#'   \item{\code{"cdf"}}{conditional CDF values, \code{newdata} as above.}
#'   \item{\code{"quantiles"}}{conditional quantiles at probabilities,
#'   \code{newdata} contains covariate columns only.}
#'   \item{\code{"moments"}}{conditional moments of order(s) \code{power},
#'   \code{newdata} contains covariate columns only.}
#' }
#'
#' @param object An object of class \code{\link{SLGP-class}}.
#' @param newdata A \code{data.frame} of evaluation points. For
#'   \code{type \%in\% c("density","cdf")} it must contain the response and
#'   covariate columns; for \code{type \%in\% c("quantiles","moments")} only the
#'   covariate columns are required.
#' @param type Character; one of \code{"density"} (default), \code{"cdf"},
#'   \code{"quantiles"}, \code{"moments"}.
#' @param probs Numeric vector of probabilities in \eqn{(0,1)}; required when
#'   \code{type = "quantiles"}.
#' @param power Numeric vector of moment orders; required when
#'   \code{type = "moments"}.
#' @param centered Logical; for \code{type = "moments"}, whether moments are
#'   centered. Default \code{FALSE}.
#' @param interpolateBasisFun Integral-approximation scheme, one of
#'   \code{"nothing"}, \code{"NN"}, \code{"WNN"} (default).
#' @param nIntegral,nDiscret Integration / discretisation resolutions passed
#'   through to the utilitarian functions. \code{nIntegral} defaults to the value
#'   recorded in \code{object} at fitting.
#' @param discrete Logical; treat the response as discrete. Defaults to the
#'   value recorded in \code{object} at fitting, so that a model fitted with
#'   \code{discrete = TRUE} is predicted on its support rather than being
#'   silently treated as continuous. Overriding either argument emits a warning.
#' @param ... Ignored.
#'
#' @return A \code{data.frame} as returned by the corresponding
#'   \code{predictSLGP_*} function.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rep(seq(0, 1, length.out = 6), each = 5))
#' d$y <- rnorm(nrow(d), mean = sin(2 * pi * d$x), sd = 0.2)
#'
#' fit <- slgp(y ~ x,
#'   data = d, method = "MAP", basisFunctionsUsed = "RFF",
#'   predictorsLower = 0, predictorsUpper = 1,
#'   responseRange = range(d$y), seed = 1,
#'   opts_BasisFun = list(nFreq = 20, MatParam = 5 / 2))
#'
#' ## Prediction grid for density and CDF evaluations:
#' grid <- expand.grid(y = seq(min(d$y), max(d$y), length.out = 100),
#'   x = c(0.25, 0.75))
#'
#' ## Predict conditional densities
#' pred_density <- predict(fit, newdata = grid, type = "density")
#'
#' ## Predict conditional cumulative distribution functions
#' pred_cdf <- predict(fit, newdata = grid, type = "cdf")
#'
#' ## Prediction locations for summaries depending only on the covariates
#' newx <- data.frame(x = c(0.25, 0.75))
#'
#' ## Predict conditional quantiles
#' pred_q <- predict(fit, newdata = newx, type = "quantiles",
#'                   probs = c(0.25, 0.5, 0.75))
#'
#' ## Predict the first two raw moments
#' pred_m <- predict(fit, newdata = newx, type = "moments", power = c(1, 2))
#'
#' @export
setMethod("predict", signature(object = "SLGP"),
          function(object, newdata,
                   type = c("density", "cdf", "quantiles", "moments"),
                   probs = NULL, power = NULL, centered = FALSE,
                   interpolateBasisFun = "WNN",
                   nIntegral = .slgp_nIntegral(object), nDiscret = 101,
                   discrete = .slgp_discrete(object), ...) {
            type <- match.arg(type)
            if (missing(newdata) || is.null(newdata))
              stop("'newdata' is required.")
            .slgp_check_discrete(object, discrete, nIntegral)
            switch(type,
                   density = .predict_density(SLGPmodel = object, newNodes = newdata,
                                              interpolateBasisFun = interpolateBasisFun,
                                              nIntegral = nIntegral, nDiscret = nDiscret,
                                              discrete = discrete),
                   cdf = .predict_cdf(SLGPmodel = object, newNodes = newdata,
                                      interpolateBasisFun = interpolateBasisFun,
                                      nIntegral = nIntegral, nDiscret = nDiscret,
                                      discrete = discrete),
                   quantiles = {
                     if (is.null(probs))
                       stop("'probs' is required when type = \"quantiles\".")
                     .predict_quantiles(SLGPmodel = object, newNodes = newdata,
                                        probs = probs,
                                        interpolateBasisFun = interpolateBasisFun,
                                        nIntegral = nIntegral, nDiscret = nDiscret,
                                        discrete = discrete)
                   },
                   moments = {
                     if (is.null(power))
                       stop("'power' is required when type = \"moments\".")
                     .predict_moments(SLGPmodel = object, newNodes = newdata,
                                      power = power, centered = centered,
                                      interpolateBasisFun = interpolateBasisFun,
                                      nIntegral = nIntegral, nDiscret = nDiscret,
                                      discrete = discrete)
                   }
            )
          })


## ----------------------------------------------------------------------------
## simulate
## ----------------------------------------------------------------------------

#' Simulate responses from a fitted SLGP model
#'
#' This method provides the standard user interface for simulating from a fitted
#' \code{\link{SLGP-class}} object.
#'
#' @param object An object of class \code{\link{SLGP-class}}.
#' @param nsim Number of samples to draw at each covariate value. Default 1.
#' @param seed Optional integer seed for reproducibility.
#' @param newdata A \code{data.frame} of covariate values.
#' @param interpolateBasisFun Integral-approximation scheme; default \code{"WNN"}.
#' @param nIntegral,nDiscret Integration / discretisation resolutions.
#' \code{nIntegral} defaults to the value recorded in \code{object} at fitting:
#' for a discrete response it is the size of the support.
#' @param type Character string; either \code{"predictive"} (default) or
#'   \code{"draws"}. With \code{"predictive"}, all responses are drawn from the
#'   posterior predictive distribution, i.e. the CDF averaged over the
#'   coefficient draws. With \code{"draws"}, each replicate is assigned one
#'   posterior draw of the SLGP and inverted against that draw's CDF, so the
#'   simulated sample also carries the between-draw variability of the fitted
#'   field. Use \code{"predictive"} to sample from the fitted model, and
#'   \code{"draws"} to propagate posterior uncertainty into downstream
#'   computations. The two coincide for models fitted with
#'   \code{method = "MAP"}, which carry a single coefficient vector.
#' @param discrete Logical: if \code{TRUE}, the response is treated as
#'   supported on the \code{nIntegral} nodes spanning \code{responseRange}, and
#'   sampling returns the first node whose CDF exceeds a uniform draw, so
#'   simulated values always lie on the support. If \code{FALSE}, the predictive
#'   CDF is inverted by linear interpolation. Defaults to the value recorded in
#'   \code{object} at fitting; overriding it emits a warning.
#' @param ... Ignored.
#'
#' @return A \code{data.frame} of sampled responses with the covariate columns
#'   from \code{newdata}.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rep(seq(0, 1, length.out = 6), each = 5))
#' d$y <- rnorm(nrow(d), mean = sin(2 * pi * d$x), sd = 0.2)
#'
#' fit <- slgp(y ~ x,
#'   data = d, method = "MAP", basisFunctionsUsed = "RFF",
#'   predictorsLower = 0, predictorsUpper = 1,
#'   responseRange = range(d$y), seed = 1,
#'   opts_BasisFun = list(nFreq = 20, MatParam = 5 / 2))
#'
#' ## Draw 10 samples from the conditional distributions at two locations
#' sim <- simulate(fit, type = "predictive",
#'                 newdata = data.frame(x = c(0.25, 0.75)), nsim = 10)
#'
#' head(sim)
#'
#' @export
setMethod("simulate", signature(object = "SLGP"),
          function(object, nsim = 1, seed = NULL, newdata,
                   type = "predictive",
                   interpolateBasisFun = "WNN",
                   nIntegral = .slgp_nIntegral(object), nDiscret = 101,
                   discrete = .slgp_discrete(object), ...) {
            type <- match.arg(type, c("predictive", "draws"))
            if (missing(newdata) || is.null(newdata))
              stop("'newdata' (covariate values) is required.")
            .slgp_check_discrete(object, discrete, nIntegral)
            .simulate_SLGP(SLGPmodel = object, newX = newdata, n = nsim,
                           type = type,
                           interpolateBasisFun = interpolateBasisFun,
                           nIntegral = nIntegral, nDiscret = nDiscret,
                           seed = seed, discrete = discrete)
          })


## ----------------------------------------------------------------------------
## update
## ----------------------------------------------------------------------------

#' Re-fit an SLGP model under a new method or with new data
#'
#' Standard-generic front end to \code{\link{retrainSLGP}}: re-estimates an
#' existing \code{\link{SLGP-class}} model, optionally with new data, under a
#' chosen estimation method, reusing the existing basis and ranges.
#'
#' @param object An object of class \code{\link{SLGP-class}}.
#' @param method Estimation method: one of \code{"none"}, \code{"Prior"},
#'   \code{"MCMC"}, \code{"MAP"}, \code{"Laplace"}.
#' @param newdata Optional new \code{data.frame}; if \code{NULL}, the model's
#'   stored data are reused.
#' @param epsilonStart Optional initial coefficient values.
#' @param interpolateBasisFun Integral-approximation scheme; default \code{"WNN"}.
#' @param nIntegral,nDiscret Integration / discretisation resolutions. If
#'   \code{nIntegral} is \code{NULL} (default), the value recorded in
#'   \code{object} is reused.
#' @param hyperparams Optional list of hyperparameters; if \code{NULL}, those of
#'   \code{object} are reused.
#' @param sigmaEstimationMethod Variance-selection method; default \code{"none"}.
#' @param seed Optional integer seed.
#' @param opts List of method-specific options (e.g. \code{stan_chains},
#'   \code{stan_iter} for MCMC, \code{ndraws} for Laplace).
#' @param discrete Logical; whether the response is treated as discrete.
#'   If \code{NULL} (default), the value recorded in \code{object} is reused, so
#'   that re-fitting a discrete model does not silently turn it into a continuous one.
#' @param trend Optional trend function.
#' @param verbose Logical; verbosity. Default \code{FALSE}.
#' @param ... Ignored.
#'
#' @return An updated object of class \code{\link{SLGP-class}}.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rep(seq(0, 1, length.out = 6), each = 5))
#' d$y <- rnorm(nrow(d), mean = sin(2 * pi * d$x), sd = 0.2)
#'
#' fit_prior <- slgp(y ~ x,
#'   data = d, method = "Prior", basisFunctionsUsed = "RFF",
#'   predictorsLower = 0, predictorsUpper = 1,
#'   responseRange = range(d$y), seed = 1,
#'   opts_BasisFun = list(nFreq = 20, MatParam = 5 / 2))
#'
#' ## Refit the same model structure by MAP
#' fit_map <- update(fit_prior, method = "MAP")
#'
#' @seealso \code{\link{retrainSLGP}} for the low-level routine.
#'
#' @export
setMethod("update", signature(object = "SLGP"),
          function(object, method, newdata = NULL, epsilonStart = NULL,
                   interpolateBasisFun = "WNN", nIntegral = NULL, nDiscret = 101,
                   discrete = NULL,
                   hyperparams = NULL, sigmaEstimationMethod = "none",
                   seed = NULL, opts = list(), trend = NULL,
                   verbose = FALSE, ...) {
            if (missing(method))
              stop("'method' is required (one of \"none\", \"Prior\", \"MCMC\", ",
                   "\"MAP\", \"Laplace\").")
            .retrain_SLGP(SLGPmodel = object, newdata = newdata,
                          epsilonStart = epsilonStart,
                          method = method, interpolateBasisFun = interpolateBasisFun,
                          nIntegral = nIntegral, nDiscret = nDiscret,
                          hyperparams = hyperparams,
                          sigmaEstimationMethod = sigmaEstimationMethod,
                          seed = seed, opts = opts, trend = trend,
                          discrete = discrete, verbose = verbose)
          })


## ----------------------------------------------------------------------------
## plot
## ----------------------------------------------------------------------------

#' Plot the estimated density field of an SLGP model
#'
#' Visualises the conditional densities implied by a fitted (or prior)
#' \code{\link{SLGP-class}} model. For a single covariate, the densities are
#' drawn as response curves at a sequence of covariate slices.
#'
#' This method draws with base graphics so that the package gains a \code{plot}
#' method without taking on a hard \pkg{ggplot2} dependency. It returns, the
#' computed prediction \code{data.frame} invisibly, so users who prefer
#' \pkg{ggplot2} can build their own figure from it.
#'
#' @param x An object of class \code{\link{SLGP-class}}.
#' @param y Ignored (present for generic compatibility).#'
#' @param newdata Optional data frame of covariate values at which to draw
#'   slices. If \code{NULL}, \code{n_slices} equally spaced covariate values
#'   spanning the model's predictor range are used.
#' @param n_slices Number of covariate slices to display (single-covariate
#'   models only). Default 8.
#' @param n_response Number of response grid points per slice. If \code{NULL}
#'   (default), 101 for a continuous response and, for a discrete one, the
#'   number of support points recorded in \code{x} at fitting, so that the
#'   plotted grid coincides with the support.
#' @param interpolateBasisFun Integral-approximation scheme; default \code{"WNN"}.
#' @param draw Integer vector selecting coefficient draws to display, or
#'   \code{"mean"} to display the pointwise mean over all stored draws.
#'   Use \code{NULL} to display all available draws.
#' @param panels Logical; if \code{TRUE} (default), draw one panel per covariate
#'   slice via \code{par(mfrow=)}. If \code{FALSE}, overlay all slices in a
#'   single plot.
#' @param discrete Logical; whether the response is treated as discrete.
#'   Defaults to the value recorded in \code{x} at fitting, in which case a step
#'   plot is drawn and the y-axis is labelled as a probability rather than a
#'   density.
#' @param ... Further graphical parameters passed to \code{matplot}.
#'
#' @return The prediction \code{data.frame}, invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(x = rep(seq(0, 1, length.out = 6), each = 5))
#' d$y <- rnorm(nrow(d), mean = sin(2 * pi * d$x), sd = 0.2)
#'
#' fit <- slgp(y ~ x,
#'   data = d, method = "MAP", basisFunctionsUsed = "RFF",
#'   predictorsLower = 0, predictorsUpper = 1,
#'   responseRange = range(d$y), seed = 1,
#'   opts_BasisFun = list(nFreq = 20, MatParam = 5 / 2))
#'
#' plot(fit, draw = "mean")
#' plot(fit, draw = c("mean", 1:5))
#' }
#'
#' @importFrom graphics matplot par
#'
#' @export
setMethod("plot", signature(x = "SLGP", y = "missing"),
          function(x, y, newdata = NULL,
                   n_slices = 6, n_response = NULL,
                   interpolateBasisFun = "WNN", draw = "mean",
                   panels = TRUE, discrete = .slgp_discrete(x), ...)  {
            if (length(x@covariateName) != 1L)
              stop("Automatic plotting is implemented for single-covariate models only; ",
                   "use predict(x, ..., type = \"density\") and plot manually.")

            ## For a discrete response the response grid is the support, so it
            ## must have exactly the resolution used at fitting.
            if (is.null(n_response))
              n_response <- if (isTRUE(discrete)) .slgp_nIntegral(x) else 101L
            .slgp_check_discrete(x, discrete, n_response)

            cov  <- x@covariateName
            resp <- x@responseName
            xr   <- c(x@predictorsRange$lower[[1]], x@predictorsRange$upper[[1]])
            rr   <- x@responseRange

            if (is.null(newdata)) {
              x_slices <- seq(xr[1], xr[2], length.out = n_slices)
            } else {
              if (!cov %in% names(newdata))
                stop("'newdata' must contain the covariate column '", cov, "'.")
              x_slices <- newdata[[cov]]
            }

            t_grid   <- seq(rr[1], rr[2], length.out = n_response)

            grid <- expand.grid(t = t_grid, xs = x_slices)
            newNodes <- data.frame(grid$t, grid$xs)
            colnames(newNodes) <- c(resp, cov)

            plot_type <- if (isTRUE(discrete)) "s" else "l"

            if (is.null(draw)) {
              use_mean <- FALSE
            } else {
              use_mean <- any(draw == "mean")
            }

            if (length(draw) == 0L && !use_mean)
              stop("'draw' does not select any available coefficient draw.")

            pred <- predict(x,
                            newdata = newNodes,
                            type = "density",
                            interpolateBasisFun = interpolateBasisFun,
                            discrete = discrete)

            pdf_cols <- grep("^pdf_", colnames(pred), value = TRUE)
            if (length(pdf_cols) == 0L)
              stop("No predicted density columns returned, cannot plot.")

            if (is.null(draw)) {
              draw <- seq_along(pdf_cols)
            } else {
              draw_num <- suppressWarnings(as.integer(draw[draw != "mean"]))
              draw_num <- draw_num[draw_num >= 1L & draw_num <= length(pdf_cols)]
              draw <- draw_num
            }

            if (length(draw) == 0L && !use_mean)
              stop("'draw' does not select any available coefficient draw.")

            dens <- list()

            if (use_mean) {
              dens[["mean"]] <- matrix(
                rowMeans(pred[, pdf_cols, drop = FALSE]),
                nrow = n_response,
                ncol = length(x_slices)
              )
            }

            if (length(draw) > 0L) {
              use_cols <- pdf_cols[draw]
              dens_draws <- lapply(use_cols, function(cc) {
                matrix(pred[[cc]], nrow = n_response, ncol = length(x_slices))
              })
              names(dens_draws) <- use_cols
              dens <- c(dens, dens_draws)
            }

            is_mean <- names(dens) == "mean"

            lty_vec <- ifelse(is_mean, 1, 2)
            lwd_vec <- ifelse(is_mean, 2.5, 1)

            if (isTRUE(panels)) {
              oldpar <- par(no.readonly = TRUE)
              on.exit(par(oldpar), add = TRUE)

              nc <- ceiling(sqrt(length(x_slices)))
              nr <- ceiling(length(x_slices) / nc)
              par(mfrow = c(nr, nc))

              ymax <- max(unlist(dens), na.rm = TRUE)

              for (j in seq_along(x_slices)) {
                mat <- do.call(cbind, lapply(dens, function(m) m[, j]))
                matplot(
                  t_grid, mat,
                  type = plot_type,
                  lty = lty_vec,
                  lwd = lwd_vec,
                  xlab = resp,
                  ylab = if (isTRUE(discrete)) "Conditional probability" else
                    "Conditional density",
                  ylim = c(0, ymax),
                  main = paste0(cov, " = ", format(x_slices[j], digits = 3)),
                  ...
                )
              }
            } else {
              ymax <- max(unlist(dens), na.rm = TRUE)
              mat <- do.call(cbind, lapply(seq_along(x_slices), function(j) {
                do.call(cbind, lapply(dens, function(m) m[, j]))
              }))

              matplot(
                t_grid, mat,
                type = plot_type, lty = lty_vec,
                lwd = lwd_vec,
                xlab = resp,
                ylab = if (isTRUE(discrete)) "Conditional probability" else
                  "Conditional density",
                ylim = c(0, ymax),
                main = paste0("SLGP density field of ", resp, " across ", cov),
                ...
              )
            }
            invisible(pred)
          })
