#' Fit a Spatial Logistic Gaussian Process model
#'
#' Builds a finite-rank Spatial Logistic Gaussian Process (SLGP) model for
#' conditional density estimation. The model can be fitted by MAP, MCMC, or
#' Laplace approximation, or initialised without fitting by setting
#' \code{method = "none"}.
#'
#'
#'
#' @param formula A formula specifying the model structure, with the response on the left-hand side and covariates on the right.
#' @param data A data frame containing the variables used in the formula.
#' @param epsilonStart Optional numeric vector of initial weights for the finite-rank GP:
#'   \eqn{Z(x,t) = \sum_{i=1}^p \epsilon_i f_i(x, t)}.
#' @param method Character string specifying the training method: one of
#'   \code{"none"}, \code{"MCMC"}, \code{"MAP"}, or \code{"Laplace"}.
#' @param basisFunctionsUsed Character string describing the basis function type:
#'   one of "inducing points", "RFF", "Discrete FF", "filling FF", or "custom cosines".
#' @param interpolateBasisFun Character string indicating how to evaluate basis functions:
#'   "nothing" (exact eval), "NN" (nearest-neighbor), or "WNN" (weighted inverse-distance). Default is "NN".
#' @param nDiscret Integer controlling the resolution of the interpolation grid (used only for "NN" or "WNN").
#' @param nIntegral Number of quadrature points used for numerical integration over the response domain.
#' @param hyperparams Optional list of hyperparameters. Should contain:
#'   \itemize{
#'     \item \code{sigma2}: signal variance
#'     \item \code{lengthscale}: vector of lengthscales (one per covariate)
#'   }
#' @param sigmaEstimationMethod Method to heuristically estimate the variance \code{sigma2}.
#'   Either "none" (default) or "heuristic".
#' @param predictorsUpper Optional numeric vector for the upper bounds of the covariates (used for scaling).
#' @param predictorsLower Optional numeric vector for the lower bounds of the covariates.
#' @param responseRange Optional numeric vector of length 2 with the lower and upper bounds of the response.
#' @param seed Optional integer for reproducibility.
#' @param opts_BasisFun List of optional configuration parameters passed to the basis function initializer.
#' @param BasisFunParam Optional list of precomputed basis function parameters.
#' @param opts Optional list of extra settings passed to inference routines (e.g., \code{stan_iter}, \code{stan_chains}, \code{ndraws}).
#' @param trend Optional function returning the trend of the transformed GP.
#'   If not provided, a zero trend is used.
#' @param discrete Logical; whether the response is treated as discrete,d efaults to \code{FALSE}.
#' When \code{TRUE}, the \code{nIntegral} quadrature nodes are taken to be the support of
#' the response and the normalising integral becomes an exact finite sum.
#' Both \code{discrete} and \code{nIntegral} are recorded on the fitted object and
#' reused as defaults by \code{\link[stats]{predict}}, \code{\link[stats]{simulate}}, \code{\link[graphics]{plot}} and \code{\link[stats]{update}}.
#' @param verbose Logical; if \code{TRUE}, print progress and diagnostic messages during computation.
#'   Defaults to \code{FALSE}.
#'
#' @return An object of S4 class \code{\link{SLGP-class}}. Standard methods are
#'   available for fitted objects, including \code{\link[base]{summary}},
#'   \code{\link[graphics]{plot}}, \code{\link[stats]{predict}},
#'   \code{\link[stats]{simulate}}, \code{\link[stats]{update}},
#'   \code{\link[stats]{coef}}, \code{\link[stats]{formula}}, and
#'   \code{\link[stats]{nobs}}.
#'
#' @importFrom stats rnorm median
#' @importFrom mvnfast rmvn
#' @export
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
#'   method = "Prior",
#'   basisFunctionsUsed = "RFF",
#'   predictorsLower = 0,
#'   predictorsUpper = 1,
#'   responseRange = range(d$y),
#'   opts_BasisFun = list(nFreq = 20, MatParam = 5 / 2),
#'   seed = 1
#' )
#'
#' fit
#' summary(fit)
#'
#' @references
#' Gautier, Athénaïs (2023). "Modelling and Predicting Distribution-Valued Fields with Applications to Inversion Under Uncertainty." Thesis, Universität Bern, Bern.
#' \url{https://boristheses.unibe.ch/4377/}
#'
slgp <- function(formula,
                 data,
                 epsilonStart = NULL,
                 method,
                 basisFunctionsUsed,
                 interpolateBasisFun = "NN",
                 nIntegral = 101,
                 nDiscret = 101,
                 hyperparams = NULL,
                 predictorsUpper = NULL,
                 predictorsLower = NULL,
                 responseRange = NULL,
                 sigmaEstimationMethod = "none",
                 seed = NULL,
                 opts_BasisFun = list(),
                 BasisFunParam = NULL,
                 opts = list(),
                 trend = NULL,
                 discrete=FALSE,
                 verbose = FALSE) {
  if(!is.null(seed)){
    set.seed(seed)
  }
  # If formula contains ".", extract all variables from the data
  if ("." %in% all.vars(formula)) {
    responseName <- as.character(formula[[2]])
    predictorNames <- names(data)[names(data) != responseName]  # Exclude the response variable
  } else {
    # Extract response and predictor variables from the formula
    responseName <- as.character(formula[[2]])
    predictorNames <- all.vars(formula)[-1]  # Exclude the response variable
  }
  # Check if all predictor variable names are in the data
  if (!all(predictorNames %in% names(data))) {
    stop("Not all predictor variables in the formula are present in the data.")
  }
  #Match arguments
  method <- match.arg(method, c("none", "Prior", "MCMC", "MAP", "Laplace"))
  basisFunctionsUsed <- match.arg(
    basisFunctionsUsed,
    c("inducing points", "RFF", "Discrete FF", "filling FF", "custom cosines")
  )
  interpolateBasisFun <- match.arg(interpolateBasisFun, c("nothing", "NN", "WNN"))
  sigmaEstimationMethod <- match.arg(sigmaEstimationMethod, c("none", "heuristic"))
  ## Bring the range of data to [0, 1]
  if(is.null(predictorsUpper)){
    predictorsUpper<- apply(data[, predictorNames, drop=FALSE], 2, max)
  }else{
    predictorsUpper<- pmax(predictorsUpper,
                           apply(data[, predictorNames, drop=FALSE], 2, max))
  }
  if(is.null(predictorsLower)){
    predictorsLower<- apply(data[, predictorNames, drop=FALSE], 2, min)
  }else{
    predictorsLower<- pmin(predictorsLower,
                           apply(data[, predictorNames, drop=FALSE], 2, min))
  }
  if(is.null(responseRange)){
    responseRange <- range(data[, responseName])
  }else{
    responseRange[1]<- min(responseRange[1], min(data[, responseName]))
    responseRange[2]<- max(responseRange[2], max(data[, responseName]))
  }
  .t_setup <- proc.time()
  normalizedData <- normalize_data(data=data, predictorNames = predictorNames, responseName = responseName,
                                   predictorsUpper = predictorsUpper, predictorsLower = predictorsLower,
                                   responseRange = responseRange)
  dimension <- ncol(normalizedData)
  if(is.null(hyperparams)){
    sigma2 <- 1
    lengthscale <- rep(0.15, dimension)
  }else{
    sigma2 <- hyperparams$sigma2
    lengthscale <- hyperparams$lengthscale
  }
  # Do we perform exact function evaluation, or we use a grid and interpolate it.
  if(interpolateBasisFun=="nothing"){
    intermediateQuantities <- pre_comput_nothing(normalizedData=normalizedData,
                                                 predictorNames=predictorNames,
                                                 responseName=responseName,
                                                 nIntegral=nIntegral)
  }
  if(interpolateBasisFun =="NN"){
    intermediateQuantities <- pre_comput_NN(normalizedData=normalizedData,
                                            predictorNames=predictorNames,
                                            responseName=responseName,
                                            nIntegral=nIntegral,
                                            nDiscret=nDiscret)
  }
  if(interpolateBasisFun == "WNN"){
    intermediateQuantities <- pre_comput_WNN(normalizedData=normalizedData,
                                             predictorNames=predictorNames,
                                             responseName=responseName,
                                             nIntegral=nIntegral,
                                             nDiscret=nDiscret)
  }
  ## Check if all options for the basis functions are provided, if not, set them to default
  opts_BasisFun <- check_basisfun_opts(basisFunctionsUsed=basisFunctionsUsed,
                                       dimension=dimension,
                                       opts_BasisFun=opts_BasisFun)
  if(is.null(BasisFunParam)){
    ## Initialise the basis functions to use
    initBasisFun <- initialize_basisfun(basisFunctionsUsed=basisFunctionsUsed,
                                        dimension=dimension,
                                        lengthscale = lengthscale,
                                        opts_BasisFun=opts_BasisFun)
  }else{
    initBasisFun <-BasisFunParam
  }

  ## Evaluate basis funs on nodes
  functionValues <- evaluate_basis_functions(parameters=initBasisFun,
                                             X=intermediateQuantities$nodes,
                                             lengthscale=lengthscale)
  if(is.null(trend)){
    trend <- function(df){return(rep(0, nrow(df)))}
    trendValues <- rep(0, nrow(functionValues))
  }else{
    dftrend <- as.data.frame(t(t(as.matrix(intermediateQuantities$nodes))*
                                 c(responseRange[2]-responseRange[1],
                                   predictorsUpper - predictorsLower)+
                                 c(responseRange[1], predictorsLower)))
    colnames(dftrend) <- c(responseName, predictorNames)
    trendValues <- trend(df=dftrend)
    rm(dftrend)
  }

  if(sigmaEstimationMethod=="heuristic"){
    Nsim <- 50
    resSim <- sapply(seq(Nsim), function(i){
      eps <- rnorm(ncol(functionValues))
      return(diff(range(functionValues%*%eps)))
    })
    sigma2 <- median(5/resSim)^2
  }


  if(discrete){
    weightQuadrature <- rep(1, nIntegral)
  }else{
    weightQuadrature <- c(1/nIntegral/2, rep(1/(nIntegral-1), nIntegral-2), 1/nIntegral/2)
  }


  # Create the data list required for the estimation method selected
  if(interpolateBasisFun == "WNN"){
    stan_model <- stanmodels$likelihoodComposed
    temp <- intermediateQuantities$indSamplesToNodes
    temp[is.na(temp)]<- 1
    temp2 <- intermediateQuantities$weightSamplesToNodes
    temp2[is.na(temp2)]<- 0

    if(!(method %in% c("none", "Prior"))){
      stan_data <- list(
        n = nrow(intermediateQuantities$indSamplesToNodes),
        nIntegral = nIntegral,
        nPredictors = nrow(intermediateQuantities$nodes)/nIntegral,
        nNeigh = ncol(intermediateQuantities$indSamplesToNodes),
        p = ncol(functionValues),
        functionValues = functionValues,
        weightMatrix = temp2,
        indMatrix =temp,
        weightQuadrature = weightQuadrature,
        Sigma = diag(sigma2, ncol(functionValues)),
        mean_x = rep(0, ncol(functionValues)),
        trendValues = trendValues)
    }
    else{stan_data <- list() }
  }else{
    stan_model <- stanmodels$likelihoodSimple
    if(interpolateBasisFun=="nothing"){
      if(!(method %in% c("none", "Prior"))){
        stan_data <- list(
          n = nrow(intermediateQuantities$indSamplesToNodes),
          nIntegral = nIntegral,
          nPredictors = as.integer(max(intermediateQuantities$indNodesToIntegral, na.rm = TRUE)),
          p = ncol(functionValues),
          meanFvalues = colMeans(functionValues[intermediateQuantities$indSamplesToNodes,]),
          functionValues = functionValues[!is.na(intermediateQuantities$indNodesToIntegral),],
          multiplicities=c(as.matrix(table(intermediateQuantities$indSamplesToPredictor))[, 1]),
          weightQuadrature = weightQuadrature,
          Sigma = diag(sigma2, ncol(functionValues)),
          mean_x = rep(0, ncol(functionValues)),
          trendValues = trendValues[!is.na(intermediateQuantities$indNodesToIntegral)])
      }
      else{stan_data <- list() }
    }
    if(interpolateBasisFun =="NN"){
      if(!(method %in% c("none", "Prior"))){
        stan_data <- list(
          n = nrow(intermediateQuantities$indSamplesToNodes),
          nIntegral = nIntegral,
          nPredictors = as.integer(max(intermediateQuantities$indNodesToIntegral, na.rm = TRUE)),
          p = ncol(functionValues),
          meanFvalues = colMeans(functionValues[intermediateQuantities$indSamplesToNodes,]),
          functionValues = functionValues[!is.na(intermediateQuantities$indNodesToIntegral),],
          weightQuadrature = weightQuadrature,
          multiplicities=c(as.matrix(table(intermediateQuantities$indNodesToIntegral[intermediateQuantities$indSamplesToNodes]))[, 1]),
          Sigma = diag(sigma2, ncol(functionValues)),
          mean_x = rep(0, ncol(functionValues)),
          trendValues = trendValues[!is.na(intermediateQuantities$indNodesToIntegral)])
      }
      else{stan_data <- list() }
    }
  }

  # Call the right estimation method
  if(is.null(epsilonStart)){
    epsilonStart <- rnorm(ncol(functionValues))
  }

  .setup_time <- proc.time() - .t_setup
  .t_est <- proc.time()
  # Estimation
  if(method=="MCMC"){
    stan_chains <- opts$stan_chains
    if(is.null(stan_chains)){
      stan_chains<-4
    }
    stan_iter <- opts$stan_iter
    if(is.null(stan_iter)){
      stan_iter<-2000
    }
    fit <- rstan::sampling(stan_model,
                           data = stan_data,
                           iter = stan_iter,
                           chains = stan_chains,
                           init = rep(list(list(epsilon = epsilonStart)),
                                      stan_chains))
    epsilon <-  rstan::extract(fit)$epsilon
    fit_summary <- rstan::summary(fit)
    rhats <- fit_summary$summary[,"Rhat"]
    if(verbose){
      cat("Convergence diagnostics:\n")
      cat(paste0("  * A R-hat close to 1 indicates a good convergence.\nHere, the dimension of the sampled values is ", stan_data$p, ".\nThe range of the R-hats is: ",
                 round(min(rhats[-c(length(rhats))]), 5), " - ",
                 round(max(rhats[-c(length(rhats))]), 5)))
    }
    ess <- fit_summary$summary[,"n_eff"]
    if(verbose){
      cat(paste0("  * Effective Sample Sizes estimates the number of independent draws from the posterior.\nHigher values are better.\nThe range of the ESS is: ",
                 round(min(ess[-c(length(ess))]), 1), " - ",
                 round(max(ess[-c(length(ess))]), 1)))
      cat(paste0("  * Checking the Bayesian Fraction of Missing Information is also a way to locate issues.\n"))
    }
    if(verbose){ rstan::check_energy(fit) }
    diagnostics <- .slgp_mcmc_diag(fit)
    logPost <- diagnostics$lp_mean
  }
  if(method=="Laplace"){
    ## The Hessian is available in closed form, so we do not ask Stan to
    ## build it by finite differences: that would cost O(p) extra gradient
    ## evaluations after convergence and dominates the fit for large p.
    fit <- rstan::optimizing(
      object = stan_model,
      data = stan_data,
      init = list(epsilon = epsilonStart),
      hessian = FALSE
    )
    if (is.null(fit$par))
      stop("Stan failed to optimise the posterior; see the message above.")
    # The MAP estimates
    mode <- fit$par
    hessian <- .slgp_hessian(stan_data, as.vector(mode),
                             interpolateBasisFun, sigma2)
    ndraws <- opts$ndraws
    if(is.null(ndraws)){
      ndraws <- 1000
    }
    nugget <- 0
    sigma <- try(solve(hessian), silent = TRUE)
    if(is.character(sigma)){
      nugget <-1e-10
      while(is.character(sigma)){
        sigma <- try(solve(hessian+nugget*diag(nrow(hessian))), silent = TRUE)
        nugget <- 10*nugget
      }
      warning("The Hessian at the mode was not invertible; a nugget of ",
              format(nugget, digits = 3), " was added. The Laplace ",
              "uncertainty should be interpreted with caution.")
    }

    epsilon <- mvnfast::rmvn(n = ndraws,
                             mu = mode,
                             sigma = sigma,
                             ncores=2)
    ## Gaussian log-density of each draw: kept so that summary() can measure
    ## how well the approximation matches the true posterior.
    logq <- mvnfast::dmvn(epsilon, mu = mode, sigma = sigma, log = TRUE)
    diagnostics <- .slgp_laplace_diag(hessian = hessian, nugget = nugget,
                                      logq = logq, ndraws = ndraws)

    logPost <- c(fit$value)
  }
  if(method=="MAP"){
    fit <- rstan::optimizing(
      object = stan_model,
      data = stan_data,
      init = list(epsilon = epsilonStart),
      hessian = FALSE,
      iter = 5000,
      tol_grad = 1e-12,
      tol_param = 1e-12,
      tol_obj = 1e-14)
    if (is.null(fit$par))
      stop("Stan failed to optimise the posterior; see the message above.")
    # The MAP estimates
    epsilon <- matrix(fit$par, nrow=1)
    diagnostics <- .slgp_map_diag(fit)
    # The log-posterior value
    logPost <- c(fit$value)
  }
  if(method=="none"){
    epsilon <- matrix(nrow=0, ncol=ncol(functionValues))
    logPost <- NaN
  }
  if (method == "Prior") {
    ndraws <- opts$ndraws
    if (is.null(ndraws)) ndraws <- 1000
    epsilon <- matrix(rnorm(ndraws * ncol(functionValues), sd = sqrt(sigma2)),
                      nrow = ndraws, ncol = ncol(functionValues))
    logPost <- NaN
  }
  if(!exists("diagnostics", inherits = FALSE)){ diagnostics <- list(scheme = method) }
  .est_time <- proc.time() - .t_est
  diagnostics$timing <- .slgp_timing(.setup_time, .est_time,
                                     if (exists("fit", inherits = FALSE)) fit else NULL,
                                     method)
  gc()
  return(SLGP(formula = formula,
              data = data,
              responseName = responseName,
              covariateName = predictorNames,
              trend=trend,
              method = method,
              predictorsRange = list(upper=predictorsUpper, lower = predictorsLower),
              responseRange = responseRange,
              p=ncol(epsilon),
              basisFunctionsUsed=basisFunctionsUsed,
              opts_BasisFun=opts_BasisFun,
              BasisFunParam=initBasisFun,
              coefficients = epsilon,
              hyperparams=list(sigma2=sigma2, lengthscale=lengthscale),
              discrete = discrete,
              nIntegral = nIntegral,
              diagnostics=diagnostics,
              logPost=logPost))
}

#'
#' Refit an SLGP model
#'
#' @description
#' `retrainSLGP()` is deprecated; use
#' \code{\link[=update,SLGP-method]{update}(object, ...)} instead.
#'
#' Re-estimates an existing SLGP model, optionally with new data, under a chosen
#' estimation method while reusing the existing basis representation and model
#' ranges.
#'
#' @param SLGPmodel An object of class \code{\link{SLGP-class}} to be retrained.
#' @param newdata Optional data frame containing new observations. If \code{NULL}, the original data is reused.
#' @param epsilonStart Optional numeric vector with initial values for the coefficients \eqn{\epsilon}.
#' @param method Character string specifying the training method: one of
#'   \code{"none"}, \code{"Prior"}, \code{"MCMC"}, \code{"MAP"}, or \code{"Laplace"}.
#' @param interpolateBasisFun Character string specifying how basis functions are evaluated:
#'   \itemize{
#'     \item \code{"nothing"} — evaluate directly at sample locations;
#'     \item \code{"NN"} — interpolate using nearest neighbor;
#'     \item \code{"WNN"} — interpolate using weighted nearest neighbors (default).
#'   }
#' @param nDiscret Integer specifying the discretisation grid size (used only if interpolation is enabled).
#' @param nIntegral Integer specifying the number of quadrature points used to approximate integrals over the response domain.
#' @param hyperparams Optional list with updated hyperparameters. Must include:
#'   \itemize{
#'     \item \code{sigma2}: signal variance;
#'     \item \code{lengthscale}: vector of lengthscales for the inputs.
#'   }
#' @param sigmaEstimationMethod Character string indicating how to estimate \code{sigma2}:
#'   either \code{"none"} (default) or \code{"heuristic"}.
#' @param seed Optional integer to set the random seed for reproducibility.
#' @param opts Optional list of additional options passed to inference routines:
#'   \code{stan_chains}, \code{stan_iter}, \code{ndraws}, etc.
#' @param trend Optional function returning the trend of the transformed GP.
#'   If not provided, a zero trend is used.
#' @param discrete Logical; whether the response is treated as discrete.
#' If \code{NULL} (default), the value recorded in \code{SLGPmodel} is reused,
#' so that refitting a discrete model does not silently turn it into a continuousone.
#' Supply \code{TRUE} or \code{FALSE} to override.
#' @param verbose Logical; if \code{TRUE}, print progress and diagnostic messages during computation.
#'   Defaults to \code{FALSE}.
#'
#' @return An updated object of class \code{\link{SLGP-class}} with retrained coefficients and updated posterior information.
#'
#' @importFrom stats rnorm
#' @importFrom mvnfast rmvn
#'
#' @seealso \code{\link[stats]{update}} for the recommended user interface.
#' @export
retrainSLGP <- function(SLGPmodel,
                        newdata=NULL,
                        epsilonStart =NULL,
                        method,
                        interpolateBasisFun="WNN",
                        nIntegral=NULL,
                        nDiscret=101,
                        hyperparams = NULL,
                        sigmaEstimationMethod = "none",
                        seed=NULL,
                        opts = list(),
                        trend=NULL,
                        discrete = NULL,
                        verbose = FALSE) {
  .Deprecated(msg = paste("retrainSLGP() is deprecated;",
                          "use update(object, ...)."))
  .retrain_SLGP(SLGPmodel = SLGPmodel,
                newdata = newdata,
                epsilonStart = epsilonStart,
                method = method,
                interpolateBasisFun = interpolateBasisFun,
                nIntegral = nIntegral,
                nDiscret = nDiscret,
                hyperparams = hyperparams,
                sigmaEstimationMethod = sigmaEstimationMethod,
                seed = seed,
                opts = opts,
                trend = trend,
                discrete = discrete,
                verbose = verbose)
}

#' Internal worker for retrainSLGP / update()
#' @noRd
.retrain_SLGP <- function(SLGPmodel,
                          newdata=NULL,
                          epsilonStart =NULL,
                          method,
                          interpolateBasisFun="WNN",
                          nDiscret=101,
                          nIntegral = NULL,
                          hyperparams = NULL,
                          sigmaEstimationMethod = "none",
                          seed=NULL,
                          opts = list(),
                          trend = NULL,
                          discrete=NULL,
                          verbose = FALSE) {
  if(!is.null(seed)){
    set.seed(seed)
  }
  responseName <- SLGPmodel@responseName
  predictorNames <- SLGPmodel@covariateName
  predictorsUpper<- SLGPmodel@predictorsRange$upper
  predictorsLower<-SLGPmodel@predictorsRange$lower
  method <- match.arg(method, c("none", "Prior", "MCMC", "MAP", "Laplace"))
  if(is.null(discrete)){
    discrete <- .slgp_discrete(SLGPmodel)
  }
  if(is.null(nIntegral)){
    nIntegral <- .slgp_nIntegral(SLGPmodel)
  }
  if(isTRUE(discrete) && nIntegral != .slgp_nIntegral(SLGPmodel)){
    warning("'nIntegral' (", nIntegral, ") differs from the value used at ",
            "fitting (", .slgp_nIntegral(SLGPmodel), "); the discrete support ",
            "grid will change.")
  }
  responseRange <-SLGPmodel@responseRange
  if(!is.null(newdata)){
    SLGPmodel@data <- newdata
  }
  .t_setup <- proc.time()
  normalizedData <- normalize_data(data=SLGPmodel@data,
                                   predictorNames = predictorNames,
                                   responseName = responseName,
                                   predictorsUpper = predictorsUpper,
                                   predictorsLower = predictorsLower,
                                   responseRange = responseRange)
  dimension <- ncol(normalizedData)

  if(!is.null(hyperparams)){
    SLGPmodel@hyperparams <- hyperparams
  }
  lengthscale <- SLGPmodel@hyperparams$lengthscale
  sigma2 <- SLGPmodel@hyperparams$sigma2

  # Do we perform exact function evaluation, or we use a grid and interpolate it.
  if(interpolateBasisFun=="nothing"){
    intermediateQuantities <- pre_comput_nothing(normalizedData=normalizedData,
                                                 predictorNames=predictorNames,
                                                 responseName=responseName,
                                                 nIntegral=nIntegral)
  }
  if(interpolateBasisFun =="NN"){
    intermediateQuantities <- pre_comput_NN(normalizedData=normalizedData,
                                            predictorNames=predictorNames,
                                            responseName=responseName,
                                            nIntegral=nIntegral,
                                            nDiscret=nDiscret)
  }
  if(interpolateBasisFun == "WNN"){
    intermediateQuantities <- pre_comput_WNN(normalizedData=normalizedData,
                                             predictorNames=predictorNames,
                                             responseName=responseName,
                                             nIntegral=nIntegral,
                                             nDiscret=nDiscret)
  }

  initBasisFun <-SLGPmodel@BasisFunParam


  ## Evaluate basis funs on nodes
  functionValues <- evaluate_basis_functions(parameters=initBasisFun,
                                             X=intermediateQuantities$nodes,
                                             lengthscale=lengthscale)
  if(sigmaEstimationMethod=="heuristic"){
    Nsim <- 50
    resSim <- sapply(seq(Nsim), function(i){
      eps <- rnorm(ncol(functionValues))
      return(diff(range(functionValues%*%eps)))
    })
    sigma2 <- median(5/resSim)^2
  }
  if(is.null(trend)){
    trend <- SLGPmodel@trend
  }
  if(is.null(trend)){
    trend <- function(df){return(rep(0, nrow(df)))}
    trendValues <- rep(0, nrow(functionValues))
  }else{
    dftrend <- as.data.frame(t(t(as.matrix(intermediateQuantities$nodes))*
                                 c(responseRange[2]-responseRange[1],
                                   predictorsUpper - predictorsLower)+
                                 c(responseRange[1], predictorsLower)))
    colnames(dftrend) <- c(responseName, predictorNames)
    trendValues <- trend(df=dftrend)
    rm(dftrend)
  }

  if(discrete){
    weightQuadrature <- rep(1, nIntegral)
  }else{
    weightQuadrature <- c(1/nIntegral/2, rep(1/(nIntegral-1), nIntegral-2), 1/nIntegral/2)
  }

  # Create the data list required for the estimation method selected
  if(interpolateBasisFun == "WNN"){
    stan_model <- stanmodels$likelihoodComposed

    temp <- intermediateQuantities$indSamplesToNodes
    temp[is.na(temp)]<- 1
    temp2 <- intermediateQuantities$weightSamplesToNodes
    temp2[is.na(temp2)]<- 0
    if(!(method %in% c("none", "Prior"))){
      stan_data <- list(
        n = nrow(intermediateQuantities$indSamplesToNodes),
        nIntegral = nIntegral,
        nPredictors = nrow(intermediateQuantities$nodes)/nIntegral,
        nNeigh = ncol(intermediateQuantities$indSamplesToNodes),
        p = ncol(functionValues),
        functionValues = functionValues,
        weightMatrix = temp2,
        indMatrix =temp,
        weightQuadrature = weightQuadrature,
        Sigma = diag(sigma2, ncol(functionValues)),
        mean_x = rep(0, ncol(functionValues)),
        trendValues = trendValues)
    }
    else{stan_data <- list() }
  }else{
    stan_model <- stanmodels$likelihoodSimple
    if(interpolateBasisFun=="nothing"){
      if(!(method %in% c("none", "Prior"))){
        stan_data <- list(
          n = nrow(intermediateQuantities$indSamplesToNodes),
          nIntegral = nIntegral,
          nPredictors = as.integer(max(intermediateQuantities$indNodesToIntegral, na.rm = TRUE)),
          p = ncol(functionValues),
          meanFvalues = colMeans(functionValues[intermediateQuantities$indSamplesToNodes,]),
          functionValues = functionValues[!is.na(intermediateQuantities$indNodesToIntegral),],
          weightQuadrature = weightQuadrature,
          multiplicities=c(as.matrix(table(intermediateQuantities$indSamplesToPredictor))[, 1]),
          Sigma = diag(sigma2, ncol(functionValues)),
          mean_x = rep(0, ncol(functionValues)),
          trendValues = trendValues[!is.na(intermediateQuantities$indNodesToIntegral)])
      }
      else{stan_data <- list() }
    }
    if(interpolateBasisFun =="NN"){
      if(!(method %in% c("none", "Prior"))){
        stan_data <- list(
          n = nrow(intermediateQuantities$indSamplesToNodes),
          nIntegral = nIntegral,
          nPredictors = as.integer(max(intermediateQuantities$indNodesToIntegral, na.rm = TRUE)),
          p = ncol(functionValues),
          meanFvalues = colMeans(functionValues[intermediateQuantities$indSamplesToNodes,]),
          functionValues = functionValues[!is.na(intermediateQuantities$indNodesToIntegral),],
          weightQuadrature = weightQuadrature,
          multiplicities=c(as.matrix(table(intermediateQuantities$indNodesToIntegral[intermediateQuantities$indSamplesToNodes]))[, 1]),
          Sigma = diag(sigma2, ncol(functionValues)),
          mean_x = rep(0, ncol(functionValues)),
          trendValues = trendValues[!is.na(intermediateQuantities$indNodesToIntegral)])
      }
      else{stan_data <- list() }
    }
  }

  # Call the right estimation method
  ## A refit starts from the coefficients already in the model when there is
  ## a single set of them (a MAP fit), which is what makes 'MAP as a starting
  ## point for MCMC' actually true.
  if(is.null(epsilonStart) && nrow(SLGPmodel@coefficients) == 1L){
    epsilonStart <- as.vector(SLGPmodel@coefficients)
  }
  if(is.null(epsilonStart)){
    epsilonStart <- rnorm(ncol(functionValues))
  }

  .setup_time <- proc.time() - .t_setup
  .t_est <- proc.time()
  # Estimation
  if(method=="MCMC"){
    stan_chains <- opts$stan_chains
    if(is.null(stan_chains)){
      stan_chains<-4
    }
    stan_iter <- opts$stan_iter
    if(is.null(stan_iter)){
      stan_iter<-2000
    }
    fit <- rstan::sampling(stan_model,
                           data = stan_data,
                           iter = stan_iter,
                           chains = stan_chains,
                           init = rep(list(list(epsilon = epsilonStart)),
                                      stan_chains))
    epsilon <-  rstan::extract(fit)$epsilon
    fit_summary <- rstan::summary(fit)
    rhats <- fit_summary$summary[,"Rhat"]
    if(verbose){
      cat("Convergence diagnostics:\n")
      cat(paste0("  * A R-hat close to 1 indicates a good convergence.\nHere, the dimension of the sampled values is ", stan_data$p, ".\nThe range of the R-hats is: ",
                 round(min(rhats[-c(length(rhats))]), 5), " - ",
                 round(max(rhats[-c(length(rhats))]), 5)))
    }
    ess <- fit_summary$summary[,"n_eff"]
    if(verbose){
      cat(paste0("  * Effective Sample Sizes estimates the number of independent draws from the posterior.\nHigher values are better.\nThe range of the ESS is: ",
                 round(min(ess[-c(length(ess))]), 1), " - ",
                 round(max(ess[-c(length(ess))]), 1)))
      cat(paste0("  * Checking the Bayesian Fraction of Missing Information is also a way to locate issues.\n"))
    }
    if(verbose){ rstan::check_energy(fit) }
    diagnostics <- .slgp_mcmc_diag(fit)
    logPost <- diagnostics$lp_mean
  }
  if(method=="Laplace"){
    ## The Hessian is available in closed form, so we do not ask Stan to
    ## build it by finite differences: that would cost O(p) extra gradient
    ## evaluations after convergence and dominates the fit for large p.
    fit <- rstan::optimizing(
      object = stan_model,
      data = stan_data,
      init = list(epsilon = epsilonStart),
      hessian = FALSE
    )
    if (is.null(fit$par))
      stop("Stan failed to optimise the posterior; see the message above.")
    # The MAP estimates
    mode <- fit$par
    hessian <- .slgp_hessian(stan_data, as.vector(mode),
                             interpolateBasisFun, sigma2)
    ndraws <- opts$ndraws
    if(is.null(ndraws)){
      ndraws <- 1000
    }
    nugget <- 0
    sigma <- try(solve(hessian), silent = TRUE)
    if(is.character(sigma)){
      nugget <-1e-10
      while(is.character(sigma)){
        sigma <- try(solve(hessian+nugget*diag(nrow(hessian))), silent = TRUE)
        nugget <- 10*nugget
      }
      warning("The Hessian at the mode was not invertible; a nugget of ",
              format(nugget, digits = 3), " was added. The Laplace ",
              "uncertainty should be interpreted with caution.")
    }
    epsilon <- mvnfast::rmvn(n = ndraws,
                             mu = mode,
                             sigma = sigma,
                             ncores=2)
    ## Gaussian log-density of each draw: kept so that summary() can measure
    ## how well the approximation matches the true posterior.
    logq <- mvnfast::dmvn(epsilon, mu = mode, sigma = sigma, log = TRUE)
    diagnostics <- .slgp_laplace_diag(hessian = hessian, nugget = nugget,
                                      logq = logq, ndraws = ndraws)
    logPost <- c(fit$value)

  }
  if(method=="MAP"){
    fit <- rstan::optimizing(
      object = stan_model,
      data = stan_data,
      init = list(epsilon = epsilonStart),
      hessian = FALSE,
      iter = 5000,
      tol_grad = 1e-12,
      tol_param = 1e-12,
      tol_obj = 1e-14)
    if (is.null(fit$par))
      stop("Stan failed to optimise the posterior; see the message above.")
    # The MAP estimates
    epsilon <- matrix(fit$par, nrow=1)
    diagnostics <- .slgp_map_diag(fit)
    logPost <- c(fit$value)

  }
  if(method=="none"){
    epsilon <- matrix(nrow=0, ncol=ncol(functionValues))
    logPost <- NaN
  }
  if(!exists("diagnostics", inherits = FALSE)){ diagnostics <- list(scheme = method) }
  if (method == "Prior") {
    ndraws <- opts$ndraws
    if (is.null(ndraws)) ndraws <- 1000
    epsilon <- matrix(rnorm(ndraws * ncol(functionValues), sd = sqrt(sigma2)),
                      nrow = ndraws, ncol = ncol(functionValues))
    logPost <- NaN
  }
  .est_time <- proc.time() - .t_est
  diagnostics$timing <- .slgp_timing(.setup_time, .est_time,
                                     if (exists("fit", inherits = FALSE)) fit else NULL,
                                     method)
  gc()
  SLGPmodel@coefficients <- epsilon
  SLGPmodel@hyperparams <- list(sigma2=sigma2, lengthscale=lengthscale)
  SLGPmodel@method <- method
  SLGPmodel@logPost <- logPost
  SLGPmodel@diagnostics <- diagnostics
  SLGPmodel@trend <- trend
  SLGPmodel@discrete <- discrete
  SLGPmodel@nIntegral <- nIntegral
  return(SLGPmodel)
}
