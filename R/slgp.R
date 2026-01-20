#' Define and can train a Spatial Logistic Gaussian Process (SLGP) model
#'
#' This function builds and trains an SLGP model based on a specified formula and data.
#' The SLGP is a finite-rank Gaussian process model for conditional density estimation,
#' trained using MAP, MCMC, Laplace approximation, or left untrained ("none").
#'
#'
#'
#' @param formula A formula specifying the model structure, with the response on the left-hand side and covariates on the right.
#' @param data A data frame containing the variables used in the formula.
#' @param epsilonStart Optional numeric vector of initial weights for the finite-rank GP:
#'   \eqn{Z(x,t) = \sum_{i=1}^p \epsilon_i f_i(x, t)}.
#' @param method Character string specifying the training method: one of \{"none", "MCMC", "MAP", "Laplace"\}.
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
#' @param trend Optional a function that returns the trend of the transformed GP (not to be estimated). If not provided, it defaults to 0.
#' @param verbose Logical; if \code{TRUE}, print progress and diagnostic messages during computation.
#'   Defaults to \code{FALSE}.
#'
#' @return An object of S4 class \code{\link{SLGP-class}}, containing:
#' \describe{
#'   \item{coefficients}{Matrix of posterior (or prior) draws of the SLGP coefficients \eqn{\epsilon_i}.}
#'   \item{hyperparams}{List of fitted or provided hyperparameters.}
#'   \item{logPost}{Log-posterior (if MAP or Laplace used).}
#'   \item{method}{Estimation method used.}
#'   \item{...}{Other internal information such as ranges, basis settings, and data.}
#' }
#'
#' @importFrom stats rnorm median
#' @importFrom mvnfast rmvn
#' @export
#'
#' @examples
#' \donttest{
#' # Load Boston housing dataset
#' library(MASS)
#' data("Boston")
#' # Set input and output ranges manually (you can also use range(Boston$age), etc.)
#' range_x <- c(0, 100)
#' range_response <- c(0, 50)
#'
#'#' #Create a SLGP model but don't fit it
#' modelPrior <- slgp(medv ~ age,        # Use a formula to specify response and covariates
#'                  data = Boston,     # Use the original Boston housing data
#'                  method = "none",    # No training
#'                  basisFunctionsUsed = "RFF",         # Random Fourier Features
#'                  sigmaEstimationMethod = "heuristic",  # Auto-tune sigma2 (more stable)
#'                  predictorsLower = range_x[1],         # Lower bound for 'age'
#'                  predictorsUpper = range_x[2],         # Upper bound for 'age'
#'                  responseRange = range_response,       # Range for 'medv'
#'                  opts_BasisFun = list(nFreq = 200,     # Use 200 Fourier features
#'                                       MatParam = 5/2), # Matern 5/2 kernel
#'                  seed = 1)                             # Reproducibility
#'
#' # Train an SLGP model using MAP estimation and RFF basis
#' modelMAP <- slgp(medv ~ age,        # Use a formula to specify response and covariates
#'                  data = Boston,     # Use the original Boston housing data
#'                  method = "MAP",    # Train using Maximum A Posteriori estimation
#'                  basisFunctionsUsed = "RFF",         # Random Fourier Features
#'                  sigmaEstimationMethod = "heuristic",  # Auto-tune sigma2 (more stable)
#'                  predictorsLower = range_x[1],         # Lower bound for 'age'
#'                  predictorsUpper = range_x[2],         # Upper bound for 'age'
#'                  responseRange = range_response,       # Range for 'medv'
#'                  opts_BasisFun = list(nFreq = 200,     # Use 200 Fourier features
#'                                       MatParam = 5/2),  # Matern 5/2 kernel
#'                  seed = 1)                             # Reproducibility
#' }
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
  #
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
    sigma2 <- median(5/resSim)
  }

  if(is.null(opts$discrete)){
    weightQuadrature <- c(1/nIntegral/2, rep(1/(nIntegral-1), nIntegral-2), 1/nIntegral/2)
  }else{
    if(opts$discrete){
      weightQuadrature <- rep(1/(nIntegral), nIntegral)
    }else{
      weightQuadrature <- c(1/nIntegral/2, rep(1/(nIntegral-1), nIntegral-2), 1/nIntegral/2)
    }
  }

  # Create the data list required for the estimation method selected
  if(interpolateBasisFun == "WNN"){
    if(file.exists("./inst/extdata/composed_model.rds")){
      stan_model <- readRDS(system.file("extdata", "composed_model.rds", package = "SLGP"))
    }else{
      stan_model_path <- system.file("stan", "likelihoodComposed.stan", package = "SLGP")
      # Load and compile the Stan model
      stan_model <- rstan::stan_model(stan_model_path, model_name = "SLGP_Likelihood_composed")
    }
    temp <- intermediateQuantities$indSamplesToNodes
    temp[is.na(temp)]<- 1
    temp2 <- intermediateQuantities$weightSamplesToNodes
    temp2[is.na(temp2)]<- 0
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
      trendValues = trendValues
    )
  }else{

    if(file.exists("./inst/extdata/simple_model.rds")){
      stan_model <- readRDS(system.file("extdata", "simple_model.rds", package = "SLGP"))
    }else{
      stan_model_path <- system.file("stan", "likelihoodSimple.stan", package = "SLGP")
      # Load and compile the Stan model
      stan_model <- rstan::stan_model(stan_model_path, model_name = "SLGP_Likelihood_simple")
    }
    if(interpolateBasisFun=="nothing"){
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
        trendValues = trendValues
      )
    }
    if(interpolateBasisFun =="NN"){
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
        trendValues = trendValues
      )
    }
  }

  # Call the right estimation method
  if(is.null(epsilonStart)){
    epsilonStart <- rnorm(ncol(functionValues))
  }

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
                           chains = stan_chains)
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
    rstan::check_energy(fit)
    logPost <- NaN # To implement later
  }
  if(method=="Laplace"){
    fit <- rstan::optimizing(
      object = stan_model,
      data = stan_data,
      hessian = TRUE  # Set to TRUE if you want to estimate the Hessian (optional)
    )
    # The MAP estimates
    mode <- fit$par
    hessian <- -fit$hessian
    ndraws <- opts$ndraws
    if(is.null(ndraws)){
      ndraws <- 1000
    }
    sigma <- try(solve(hessian), silent = TRUE)
    if(is.character(sigma)){
      nugget <-1e-10
      while(is.character(sigma)){
        sigma <- try(solve(hessian+nugget*diag(nrow(hessian))))
        nugget <- 10*nugget
      }
    }

    epsilon <- mvnfast::rmvn(n = ndraws,
                             mu = mode,
                             sigma = sigma,
                             ncores=2)

    logPost <- c(fit$value)
  }
  if(method=="MAP"){
    fit <- rstan::optimizing(
      object = stan_model,
      data = stan_data,
      hessian = FALSE,
      iter = 5000,
      tol_grad = 1e-12,
      tol_param = 1e-12,
      tol_obj = 1e-14)
    # The MAP estimates
    epsilon <- matrix(fit$par, nrow=1)
    # The log-posterior value
    logPost <- c(fit$value)
  }
  if(method=="none"){
    epsilon <- matrix(nrow=0, ncol=ncol(functionValues))
    logPost <- NaN
  }
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
              logPost=logPost))
}

#' Retrain a fitted SLGP model with new data and/or estimation method
#'
#' This function retrains an existing SLGP model using either a Bayesian MCMC estimation,
#' a Maximum A Posteriori (MAP) estimation, or a Laplace approximation. The model can be retrained
#' using new data, new inference settings, or updated hyperparameters. It reuses the structure and
#' basis functions from the original model.
#'
#' @param SLGPmodel An object of class \code{\link{SLGP-class}} to be retrained.
#' @param newdata Optional data frame containing new observations. If \code{NULL}, the original data is reused.
#' @param epsilonStart Optional numeric vector with initial values for the coefficients \eqn{\epsilon}.
#' @param method Character string specifying the estimation method: one of \{"MCMC", "MAP", "Laplace"\}.
#' @param interpolateBasisFun Character string specifying how basis functions are evaluated:
#'   \itemize{
#'     \item \code{"nothing"} — evaluate directly at sample locations;
#'     \item \code{"NN"} — interpolate using nearest neighbor;
#'     \item \code{"WNN"} — interpolate using weighted nearest neighbors (default).
#'   }
#' @param nDiscret Integer specifying the discretization grid size (used only if interpolation is enabled).
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
#' @param trend Optional a function that returns the trend of the transformed GP (not to be estimated). If not provided, it defaults to 0.
#' @param verbose Logical; if \code{TRUE}, print progress and diagnostic messages during computation.
#'   Defaults to \code{FALSE}.
#'
#' @return An updated object of class \code{\link{SLGP-class}} with retrained coefficients and updated posterior information.
#'
#' @importFrom stats rnorm
#' @importFrom mvnfast rmvn
#'
#' @export
#'
#' @references
#' Gautier, A. (2023). *Modelling and Predicting Distribution-Valued Fields with Applications to Inversion Under Uncertainty*.
#' PhD Thesis, Universität Bern. \url{https://boristheses.unibe.ch/4377/}
#'
#' @examples
#' \donttest{
#' # Load Boston housing dataset
#' library(MASS)
#' data("Boston")
#' range_x <- c(0, 100)
#' range_response <- c(0, 50)
#'
#' #Create a SLGP model but don't fit it
#' modelPrior <- slgp(medv ~ age,        # Use a formula to specify response and covariates
#'                  data = Boston,     # Use the original Boston housing data
#'                  method = "none",    # No training
#'                  basisFunctionsUsed = "RFF",         # Random Fourier Features
#'                  sigmaEstimationMethod = "heuristic",  # Auto-tune sigma2 (more stable)
#'                  predictorsLower = range_x[1],         # Lower bound for 'age'
#'                  predictorsUpper = range_x[2],         # Upper bound for 'age'
#'                  responseRange = range_response,       # Range for 'medv'
#'                  opts_BasisFun = list(nFreq = 200,     # Use 200 Fourier features
#'                                       MatParam = 5/2), # Matern 5/2 kernel
#'                  seed = 1)                             # Reproducibility
#' #Retrain using the Boston Housing dataset and a Laplace approximation scheme
#' modelLaplace <- retrainSLGP(SLGPmodel=modelPrior,
#'                             newdata = Boston,
#'                             method="Laplace")
#' }
#'
retrainSLGP <- function(SLGPmodel,
                        newdata=NULL,
                        epsilonStart =NULL,
                        method,
                        interpolateBasisFun="WNN",
                        nIntegral=101,
                        nDiscret=101,
                        hyperparams = NULL,
                        sigmaEstimationMethod = "none",
                        seed=NULL,
                        opts = list(),
                        trend=NULL,
                        verbose = FALSE) {
  if(!is.null(seed)){
    set.seed(seed)
  }
  responseName <- SLGPmodel@responseName
  predictorNames <- SLGPmodel@covariateName
  predictorsUpper<- SLGPmodel@predictorsRange$upper
  predictorsLower<-SLGPmodel@predictorsRange$lower

  responseRange <-SLGPmodel@responseRange

  if(!is.null(newdata)){
    SLGPmodel@data <- newdata
  }

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
    sigma2 <- median(5/resSim)
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
  if(is.null(opts$discrete)){
    weightQuadrature <- c(1/nIntegral/2, rep(1/(nIntegral-1), nIntegral-2), 1/nIntegral/2)
  }else{
    if(opts$discrete){
      weightQuadrature <- rep(1/(nIntegral), nIntegral)
    }else{
      weightQuadrature <- c(1/nIntegral/2, rep(1/(nIntegral-1), nIntegral-2), 1/nIntegral/2)
    }
  }

  # Create the data list required for the estimation method selected
  if(interpolateBasisFun == "WNN"){
    stan_model <- stanmodels$likelihoodComposed

    temp <- intermediateQuantities$indSamplesToNodes
    temp[is.na(temp)]<- 1
    temp2 <- intermediateQuantities$weightSamplesToNodes
    temp2[is.na(temp2)]<- 0
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
      trendValues = trendValues
    )
  }else{
    stan_model <- stanmodels$likelihoodSimple

    if(interpolateBasisFun=="nothing"){
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
        trendValues = trendValues[!is.na(intermediateQuantities$indNodesToIntegral)]
      )
    }
    if(interpolateBasisFun =="NN"){
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
        trendValues = trendValues[!is.na(intermediateQuantities$indNodesToIntegral)]
      )
    }
  }


  # Call the right estimation method
  if(is.null(epsilonStart)){
    epsilonStart <- rnorm(ncol(functionValues))
  }

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
                           chains = stan_chains)
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
    rstan::check_energy(fit)
    logPost <- NaN
  }
  if(method=="Laplace"){
    fit <- rstan::optimizing(
      object = stan_model,
      data = stan_data,
      hessian = TRUE  # Set to TRUE if you want to estimate the Hessian (optional)
    )
    # The MAP estimates
    mode <- fit$par
    hessian <- -fit$hessian
    ndraws <- opts$ndraws
    if(is.null(ndraws)){
      ndraws <- 1000
    }
    sigma <- try(solve(hessian), silent = TRUE)
    if(is.character(sigma)){
      nugget <-1e-10
      while(is.character(sigma)){
        sigma <- try(solve(hessian+nugget*diag(nrow(hessian))))
        nugget <- 10*nugget
      }
    }

    epsilon <- mvnfast::rmvn(n = ndraws,
                             mu = mode,
                             sigma = sigma,
                             ncores=2)

    logPost <- c(fit$value)

  }
  if(method=="MAP"){
    fit <- rstan::optimizing(
      object = stan_model,
      data = stan_data,
      hessian = FALSE,
      iter = 5000,
      tol_grad = 1e-12,
      tol_param = 1e-12,
      tol_obj = 1e-14)
    # The MAP estimates
    epsilon <- matrix(fit$par, nrow=1)
    logPost <- c(fit$value)

  }
  if(method=="none"){
    epsilon <- matrix(nrow=0, ncol=ncol(functionValues))
    logPost <- NaN
  }
  gc()
  SLGPmodel@coefficients <- epsilon
  SLGPmodel@hyperparams <- list(sigma2=sigma2, lengthscale=lengthscale)
  SLGPmodel@method <- method
  SLGPmodel@logPost <- logPost
  SLGPmodel@trend <- trend
  return(SLGPmodel)
}
