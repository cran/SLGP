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
#' d <- data.frame(
#'   x = rep(seq(0, 1, length.out = 6), each = 5)
#' )
#' d$y <- rnorm(nrow(d), mean = sin(2 * pi * d$x), sd = 0.2)
#'
#' fit <- slgp(
#'   y ~ x,
#'   data = d,
#'   method = "none",
#'   basisFunctionsUsed = "RFF",
#'   predictorsLower = 0,
#'   predictorsUpper = 1,
#'   responseRange = range(d$y),
#'   opts_BasisFun = list(nFreq = 20, MatParam = 5 / 2),
#'   seed = 1
#' )
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
      "  in ", .slgp_fmt_range(x@responseRange), "\n", sep = "")
  cat("  Covariate(s)   : ", paste(x@covariateName, collapse = ", "), "\n", sep = "")
  cat("  Estimation     : ", method,
      if (!.slgp_is_fitted(x)) "  (prior, not fitted)" else "", "\n", sep = "")
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
#' Extends \code{\link{print}} with fit diagnostics where available: the
#' log-posterior at the mode (MAP / Laplace), the number of coefficient draws,
#' and basic structural information. The returned object has class
#' \code{"summary.SLGP"} and its own print method.
#'
#' @param object An object of class \code{\link{SLGP-class}}.
#' @param ... Ignored.
#'
#' @return An object of class \code{"summary.SLGP"} (a list), returned invisibly.
#'
#' @export
setMethod("summary", signature(object = "SLGP"), function(object, ...) {
  ndraws <- .slgp_ndraws(object)
  nobs   <- if (nrow(object@data)) nrow(object@data) else 0L
  ndistinct <- if (length(object@covariateName) && nobs > 0L) {
    nrow(unique(object@data[, object@covariateName, drop = FALSE]))
  } else NA_integer_

  out <- list(
    formula            = object@formula,
    responseName       = object@responseName,
    responseRange      = object@responseRange,
    covariateName      = object@covariateName,
    method             = if (length(object@method)) object@method else NA_character_,
    fitted             = .slgp_is_fitted(object),
    basisFunctionsUsed = .slgp_basis_label(object),
    p                  = object@p,
    hyperparams        = object@hyperparams,
    nobs               = nobs,
    ndistinct          = ndistinct,
    ndraws             = ndraws,
    logPost            = if (length(object@logPost) == 1L) object@logPost else NA_real_
  )
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
      "  in ", .slgp_fmt_range(x$responseRange), "\n", sep = "")
  cat("  Covariate(s)   : ", paste(x$covariateName, collapse = ", "), "\n", sep = "")
  cat("  Estimation     : ", x$method,
      if (!isTRUE(x$fitted)) "  (prior, not fitted)" else "", "\n", sep = "")
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
  cat("\nDiagnostics\n")
  if (x$ndraws > 0L)
    cat("  Coefficient draws : ", x$ndraws, "\n", sep = "")
  if (!is.na(x$logPost) && is.finite(x$logPost)) {
    cat("  Log-posterior     : ", format(x$logPost, digits = 6),
        " (at mode for MAP/Laplace)\n", sep = "")
  } else {
    cat("  Log-posterior     : not available for method '", x$method, "'\n", sep = "")
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
#'   through to the workhorse functions.
#' @param discrete Logical; treat the response as discrete. Default \code{FALSE}.
#' @param ... Ignored.
#'
#' @return A \code{data.frame} as returned by the corresponding
#'   \code{predictSLGP_*} function.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(
#'   x = rep(seq(0, 1, length.out = 6), each = 5)
#' )
#' d$y <- rnorm(nrow(d), mean = sin(2 * pi * d$x), sd = 0.2)
#'
#' fit <- slgp(
#'   y ~ x,
#'   data = d,
#'   method = "MAP",
#'   basisFunctionsUsed = "RFF",
#'   predictorsLower = 0,
#'   predictorsUpper = 1,
#'   responseRange = range(d$y),
#'   opts_BasisFun = list(nFreq = 20, MatParam = 5 / 2),
#'   seed = 1
#' )
#'
#' ## Prediction grid for density and CDF evaluations:
#' grid <- expand.grid(
#'   y = seq(min(d$y), max(d$y), length.out = 100),
#'   x = c(0.25, 0.75)
#' )
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
                   nIntegral = 101, nDiscret = 101, discrete = FALSE, ...) {
            type <- match.arg(type)
            if (missing(newdata) || is.null(newdata))
              stop("'newdata' is required.")

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
#' @param discrete Logical; treat the response as discrete. Default \code{FALSE}.
#' @param ... Ignored.
#'
#' @return A \code{data.frame} of sampled responses with the covariate columns
#'   from \code{newdata}.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(
#'   x = rep(seq(0, 1, length.out = 6), each = 5)
#' )
#' d$y <- rnorm(nrow(d), mean = sin(2 * pi * d$x), sd = 0.2)
#'
#' fit <- slgp(
#'   y ~ x,
#'   data = d,
#'   method = "MAP",
#'   basisFunctionsUsed = "RFF",
#'   predictorsLower = 0,
#'   predictorsUpper = 1,
#'   responseRange = range(d$y),
#'   opts_BasisFun = list(nFreq = 20, MatParam = 5 / 2),
#'   seed = 1
#' )
#'
#' ## Draw 10 samples from the conditional distributions at two locations
#' sim <- simulate(
#'   fit,
#'   newdata = data.frame(x = c(0.25, 0.75)),
#'   nsim = 10
#' )
#'
#' head(sim)
#'
#' @export
setMethod("simulate", signature(object = "SLGP"),
          function(object, nsim = 1, seed = NULL, newdata,
                   interpolateBasisFun = "WNN",
                   nIntegral = 101, nDiscret = 101, discrete = FALSE, ...) {
            if (missing(newdata) || is.null(newdata))
              stop("'newdata' (covariate values) is required.")
            .simulate_SLGP(SLGPmodel = object, newX = newdata, n = nsim,
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
#' @param method Estimation method: one of \code{"MCMC"}, \code{"MAP"},
#'   \code{"Laplace"}.
#' @param newdata Optional new \code{data.frame}; if \code{NULL}, the model's
#'   stored data are reused.
#' @param epsilonStart Optional initial coefficient values.
#' @param interpolateBasisFun Integral-approximation scheme; default \code{"WNN"}.
#' @param nIntegral,nDiscret Integration / discretisation resolutions.
#' @param hyperparams Optional list of hyperparameters; if \code{NULL}, those of
#'   \code{object} are reused.
#' @param sigmaEstimationMethod Variance-selection method; default \code{"none"}.
#' @param seed Optional integer seed.
#' @param opts List of method-specific options (e.g. \code{stan_chains},
#'   \code{stan_iter} for MCMC, \code{ndraws} for Laplace).
#' @param discrete Logical; whether the response is treated as discrete.
#'   Default \code{FALSE}.
#' @param trend Optional trend function.
#' @param verbose Logical; verbosity. Default \code{FALSE}.
#' @param ... Ignored.
#'
#' @return An updated object of class \code{\link{SLGP-class}}.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(
#'   x = rep(seq(0, 1, length.out = 6), each = 5)
#' )
#' d$y <- rnorm(nrow(d), mean = sin(2 * pi * d$x), sd = 0.2)
#'
#' prior <- slgp(
#'   y ~ x,
#'   data = d,
#'   method = "none",
#'   basisFunctionsUsed = "RFF",
#'   predictorsLower = 0,
#'   predictorsUpper = 1,
#'   responseRange = range(d$y),
#'   opts_BasisFun = list(nFreq = 20, MatParam = 5 / 2),
#'   seed = 1
#' )
#'
#' ## Refit the same model structure by MAP
#' fit_map <- update(prior, method = "MAP")
#'
#' @seealso \code{\link{retrainSLGP}} for the low-level routine.
#'
#' @export
setMethod("update", signature(object = "SLGP"),
          function(object, method, newdata = NULL, epsilonStart = NULL,
                   interpolateBasisFun = "WNN", nIntegral = 101, nDiscret = 101,
                   discrete = FALSE,
                   hyperparams = NULL, sigmaEstimationMethod = "none",
                   seed = NULL, opts = list(), trend = NULL,
                   verbose = FALSE, ...) {
            if (missing(method))
              stop("'method' is required (one of \"MCMC\", \"MAP\", \"Laplace\").")
            .retrain_SLGP(SLGPmodel = object, newdata = newdata, epsilonStart = epsilonStart,
                        method = method, interpolateBasisFun = interpolateBasisFun,
                        nIntegral = nIntegral, nDiscret = nDiscret,
                        hyperparams = hyperparams,
                        sigmaEstimationMethod = sigmaEstimationMethod,
                        seed = seed, opts = opts, trend = trend, discrete=discrete,
                        verbose = verbose)
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
#' @param n_response Number of response grid points per slice. Default 101.
#' @param interpolateBasisFun Integral-approximation scheme; default \code{"WNN"}.
#' @param draw Integer vector selecting coefficient draws to display, or
#'   \code{"mean"} to display the pointwise mean over all stored draws.
#'   Use \code{NULL} to display all available draws.
#' @param panels Logical; if \code{TRUE} (default), draw one panel per covariate
#'   slice via \code{par(mfrow=)}. If \code{FALSE}, overlay all slices in a
#'   single plot.
#' @param discrete Logical; whether the response is treated as discrete.
#'   Default \code{FALSE}.
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
#' fit <- slgp(
#'   y ~ x, data = d, method = "Laplace", basisFunctionsUsed = "RFF",
#'   predictorsLower = 0, predictorsUpper = 1, responseRange = range(d$y),
#'   opts_BasisFun = list(nFreq = 20, MatParam = 5 / 2),
#'   seed = 1, opts=list(ndraws=5)
#' )
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
                   n_slices = 6, n_response = 101,
                   interpolateBasisFun = "WNN", draw = "mean",
                   panels = TRUE, discrete = FALSE, ...) {
            if (length(x@covariateName) != 1L)
              stop("Automatic plotting is implemented for single-covariate models only; ",
                   "use predict(x, ..., type = \"density\") and plot manually.")

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
                  ylab = if (isTRUE(discrete)) "Conditional probability" else "Conditional density",
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
                ylab = if (isTRUE(discrete)) "Conditional probability" else "Conditional density",
                ylim = c(0, ymax),
                main = paste0("SLGP density field of ", resp, " across ", cov),
                ...
              )
            }
            invisible(pred)
          })
