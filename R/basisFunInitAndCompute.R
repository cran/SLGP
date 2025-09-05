#' Check basis function parameters
#'
#' Checks and completes the parameter list for a given basis function type.
#'
#' @param basisFunctionsUsed Character. Type of basis function to use.
#'   One of: "inducing points", "RFF", "Discrete FF", "filling FF", "custom cosines".
#' @param dimension Integer. The dimension of the input space (typically \eqn{[\mathbf{x}, t]}).
#' @param opts_BasisFun List. Options specific to the chosen basis function.
#' Users can refer to the documentation of specific basis function initialization functions
#' (e.g., \code{\link{initialize_basisfun_inducingpt}}, \code{\link{initialize_basisfun_RFF}},
#'  \code{\link{initialize_basisfun_fillingRFF}}, \code{\link{initialize_basisfun_discreteFF}}, etc.) for details on the available options.
#'
#' @return A completed list of options specific to the chosen basis function.
#'
#' @keywords internal
#'
#' @importFrom stats runif
#'
check_basisfun_opts <- function(basisFunctionsUsed,
                                dimension,
                                opts_BasisFun=list()) {
  # Check if basisFunctionsUsed is valid
  valid_types <- c("inducing points", "RFF", "Discrete FF", "filling FF", "custom cosines")
  if (!basisFunctionsUsed %in% valid_types) {
    stop("Invalid basisFunctionsUsed. Choose from: 'inducing points', 'RFF', 'Discrete FF', 'filling FF', 'custom sines'")
  }

  opts_BasisFunClean <- list()
  if (basisFunctionsUsed == "inducing points") {
    # Add custom parameters
    if(is.null(opts_BasisFun$numberPoints)){
      if(!is.null(opts_BasisFun$pointscoord)){
        opts_BasisFunClean$numberPoints <- nrow(opts_BasisFun$pointscoord)
      }else{
        opts_BasisFunClean$numberPoints <- 10
      }
      warning("You did not specify a number of anchoring points 'numberPoints' in 'opts_BasisFun', defaulted to 10")
    }else{
      opts_BasisFunClean$numberPoints <- opts_BasisFun$numberPoints
    }
    if(is.null(opts_BasisFun$kernel)){
      opts_BasisFunClean$kernel <- "Mat52"
      warning("You did not specify a kernel type 'kernel' in 'opts_BasisFun', defaulted to Matern 5/2")
    }else{
      opts_BasisFunClean$kernel <- opts_BasisFun$kernel
    }
    if(is.null(opts_BasisFun$pointscoord)){
      opts_BasisFunClean$pointscoord <- matrix(runif(opts_BasisFun$numberPoints*dimension), ncol=dimension)
    }else{
      opts_BasisFunClean$pointscoord <- opts_BasisFun$pointscoord
    }
  }
  if (basisFunctionsUsed == "RFF"){
    if(is.null(opts_BasisFun$MatParam)){
      opts_BasisFunClean$MatParam <- 5/2
      warning("You did not specify a valid Matern parameter 'MatParam' in 'opts_BasisFun', defaulted to Matern 5/2.")
    }else{
      if(!is.infinite(opts_BasisFun$MatParam)){
        if(0!=((opts_BasisFun$MatParam-1/2) %%1)){
          opts_BasisFunClean$MatParam <- 5/2
          warning("You did not specify a valid Matern parameter 'MatParam' in 'opts_BasisFun', defaulted to Matern 5/2.")
        }
      }
      opts_BasisFunClean$MatParam <- opts_BasisFun$MatParam
    }
    if(is.null(opts_BasisFun$nFreq)){
      opts_BasisFunClean$nFreq <- 5
      warning("You did not specify a valid number of frequencies 'nFreq' in 'opts_BasisFun', defaulted to 5 (gives 10 basis functions).")
    }else{
      opts_BasisFunClean$nFreq <- opts_BasisFun$nFreq
    }
  }
  if (basisFunctionsUsed == "filling FF") {
    opts_BasisFunClean$seed <- opts_BasisFun$seed
    if(is.null(opts_BasisFun$MatParam) | 0!=((opts_BasisFun$MatParam-1/2) %%1)){
      opts_BasisFunClean$MatParam <- 5/2
      warning("You did not specify a valid Matern parameter 'MatParam' in 'opts_BasisFun', defaulted to Matern 5/2.")
    }else{
      opts_BasisFunClean$MatParam <- opts_BasisFun$MatParam
    }
    if(is.null(opts_BasisFun$nFreq)){
      opts_BasisFunClean$nFreq <- 5
      warning("You did not specify a valid number of frequencies 'nFreq' in 'opts_BasisFun', defaulted to 5 (gives 10 basis functions).")
    }else{
      opts_BasisFunClean$nFreq <- opts_BasisFun$nFreq
    }
  }
  if (basisFunctionsUsed == "Discrete FF") {
    if(is.null(opts_BasisFun$maxOrdert)){
      opts_BasisFunClean$maxOrdert <- 2
      warning("You did not specify a maximum order for t 'maxOrdert' in 'opts_BasisFun',  defaulted to 2.")
    }else{
      opts_BasisFunClean$maxOrdert <- opts_BasisFun$maxOrdert
    }
    if(is.null(opts_BasisFun$maxOrderx)){
      opts_BasisFunClean$maxOrderx <- 2
      warning("You did not specify a maximum order for x 'maxOrderx' in 'opts_BasisFun',  defaulted to 2.")
    }else{
      opts_BasisFunClean$maxOrderx <- opts_BasisFun$maxOrderx
    }
  }

  if (basisFunctionsUsed == "custom cosines") {
    if(is.null(opts_BasisFun$freq)){
      stop("You did not specify a vector of frequencies 'freq' in 'opts_BasisFun' for your custom cosines.")
    }else{
      opts_BasisFunClean$freq <- opts_BasisFun$freq
    }
    if(is.null(opts_BasisFun$coef)){
      opts_BasisFunClean$coef <- 1
    }else{
      opts_BasisFunClean$coef <- opts_BasisFun$coef
    }
    if(is.null(opts_BasisFun$offset)){
      stop("You did not specify a vector of offsets 'offset' in 'opts_BasisFun' for your custom cosines.")
    }else{
      opts_BasisFunClean$offset <- opts_BasisFun$offset
    }
    # Add custom sines-specific parameters
    temp <- list(freq=opts_BasisFun$freq,
                 coef=opts_BasisFun$coef,
                 offset=opts_BasisFun$offset)
  }
  opts_BasisFunClean$lengthscale <- opts_BasisFun$lengthscale
  return(opts_BasisFunClean)
}


#' Initialize basis function parameters
#'
#' Initializes the parameter list needed for a basis function.
#'
#' @param basisFunctionsUsed Character. The type of basis function to use.
#' One of: "inducing points", "RFF", "Discrete FF", "filling FF", "custom cosines".
#'
#' @param dimension Integer. Dimension of the input space \eqn{[\mathbf{x},\,t]}{[x, t]}.
#' @param lengthscale Numeric vector. Lengthscales used for scaling the input space.
#' @param opts_BasisFun List. Optional. Additional options specific to the chosen basis function.
#' If the type is "custom cosines", the basis functions considered are \eqn{ coef\cos(freq^\top [x, t] + offset) }
#' and the user must provide three vectors: \code{opts_BasisFun$freq}, \code{opts_BasisFun$offset} and \code{opts_BasisFun$coef}.
#' Users can refer to the documentation of specific basis function initialization functions
#' (e.g., \code{\link{initialize_basisfun_inducingpt}}, \code{\link{initialize_basisfun_RFF}},
#' \code{\link{initialize_basisfun_fillingRFF}}, \code{\link{initialize_basisfun_discreteFF}}, etc.) for details on the available options.
#'
#' @return A list of initialized basis function parameters.
#'
#' @keywords internal
#'
initialize_basisfun <- function(basisFunctionsUsed, dimension, lengthscale, opts_BasisFun=list()) {
  # Initialize parameters based on basisFunctionsUsed
  parameters <- list(
    basisFunctionsUsed = basisFunctionsUsed,
    dimension = dimension
  )

  # Additional initialization based on basisFunctionsUsed
  if (basisFunctionsUsed == "inducing points") {
    temp <- initialize_basisfun_inducingpt(dimension=dimension,
                                           kernel = opts_BasisFun$kernel,
                                           lengthscale= lengthscale,
                                           pointscoord = opts_BasisFun$pointscoord,
                                           numberPoints = opts_BasisFun$numberPoints)
    parameters[["kernel"]] <-opts_BasisFun$kernel
    parameters[["sqrtInv_covMat"]] <- temp$sqrtInv_covMat
    parameters[["sqrt_covMat"]] <- temp$sqrt_covMat
    parameters[["lengthCoord"]] <- temp$lengthCoord
  }
  if (basisFunctionsUsed %in% c("RFF", "filling FF", "Discrete FF", "custom cosines")) {
    if (basisFunctionsUsed == "RFF") {
      temp <- initialize_basisfun_RFF(dimension=dimension,
                                      nFreq=opts_BasisFun$nFreq,
                                      MatParam = opts_BasisFun$MatParam,
                                      lengthscale=lengthscale)
    }
    if (basisFunctionsUsed == "filling FF") {
      temp <- initialize_basisfun_fillingRFF(dimension=dimension,
                                             nFreq=opts_BasisFun$nFreq,
                                             MatParam = opts_BasisFun$MatParam,
                                             lengthscale=lengthscale,
                                             seed=opts_BasisFun$seed)
    }
    if (basisFunctionsUsed == "Discrete FF") {
      temp <- initialize_basisfun_discreteFF(dimension=dimension,
                                             maxOrdert=opts_BasisFun$maxOrdert,
                                             maxOrderx=opts_BasisFun$maxOrderx)
    }

    if (basisFunctionsUsed == "custom cosines") {
      # Add custom sines-specific parameters
      temp <- list(freq=opts_BasisFun$freq,
                   coef=opts_BasisFun$coef,
                   offset=opts_BasisFun$offset)
    }
    parameters[["freq"]] <- temp$freq
    parameters[["coef"]] <- temp$coef
    parameters[["offset"]] <- temp$offset
  }
  return(parameters)
}

#' Initialize parameters for inducing-point basis functions
#'
#' Computes kernel matrix and its decompositions for use in inducing-point basis functions.
#'
#' @param dimension Integer. Input (\eqn{[\mathbf{x},\,t]}{[x, t]}) dimension.
#' @param kernel Character. Kernel type ("Exp", "Mat32", "Mat52", "Gaussian").
#' @param lengthscale Numeric vector. Lengthscales used for scaling the input space.
#' @param pointscoord Optional matrix of inducing point coordinates.
#' If none is provided, we sample them uniformly in the unit hypercube.
#' @param numberPoints Integer. Number of inducing points
#' (used if `pointscoord` is NULL).
#'
#' @return List with kernel square root and inverse root matrices, and scaled coordinates.
#'
#' @keywords internal
#'
initialize_basisfun_inducingpt <- function(dimension, kernel = "Mat52", lengthscale, pointscoord = NULL, numberPoints =NULL){

  temp <- t(t(pointscoord)/lengthscale)
  distances <- crossdist(temp, temp)
  if(kernel=="Exp"){
    K <- exp(-distances)
  }
  if(kernel=="Mat32"){
    K <- (1+sqrt(3)*distances)*exp(-sqrt(3)*distances)
  }
  if(kernel=="Mat52"){
    K <- (1+sqrt(5)*distances+5*distances^2/3)*exp(-sqrt(5)*distances)
  }
  if(kernel=="Gaussian"){
    K <- exp(-distances^2/2)
  }
  eig <- eigen(K)
  eig$values <- abs(eig$values)/2+eig$values/2+1e-20
  if(any(eig$values==0)){
    stop("The inducing points/kernel provided make the covariance matrix ill-conditioned. Select different points or regularise, please.")
  }
  Ksqrt <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  Kinvsqrt <- eig$vectors %*% diag(1/sqrt(eig$values)) %*% t(eig$vectors)

  # for Y ~ N(0, K), Y %*% Kinvsqrt ~ N(0, I)
  return(list(sqrtInv_covMat=Kinvsqrt, sqrt_covMat=Ksqrt, lengthCoord=temp))
}


#' Initialize parameters basis functions based on Random Fourier Features
#'
#' Draws parameters for standard RFF approximating a Matérn kernel.
#'
#' @param dimension Integer. Input (\eqn{[\mathbf{x},\,t]}{[x, t]}) dimension.
#' @param nFreq Integer. Number of frequency vectors to be considered.
#' @param MatParam Numeric. Matérn smoothness parameter (default = 5/2).
#' @param lengthscale Numeric vector. Lengthscales used for scaling the input space.
#'
#' @return List with frequency, offset, and coefficient parameters.
#'
#' @keywords internal
#'
#' @importFrom mvnfast rmvt
#'
initialize_basisfun_RFF <- function(dimension, nFreq, MatParam = 5/2, lengthscale) {
  # Check if basisFunctionsUsed is valid
  if (!requireNamespace("mvnfast", quietly = TRUE)) {
    stop("Package 'mvnfast' could not be used")
  }
  if(!is.infinite(MatParam)){
    freq <- mvnfast::rmvt(n=nFreq, mu=rep(0, dimension),
                          sigma=diag(dimension), df=2*MatParam)
  }else{
    freq <-  mvnfast::rmvn(n=nFreq, mu=rep(0, dimension),
                           sigma=4*diag(dimension))
  }
  freq <- rbind(freq, freq)
  offset <- c(rep(0, nFreq), rep(-pi/2, nFreq))
  coef <- rep(1/sqrt(nFreq), 2*nFreq)
  return(list(freq=freq, offset=offset, coef=coef))
}


#' Initialize space-filling Random Fourier Features
#'
#' Initializes RFF parameters with LHS-optimized frequency directions.
#'
#' @param dimension Integer. Input (\eqn{[\mathbf{x},\,t]}{[x, t]}) dimension.
#' @param nFreq Integer. Number of frequency vectors to be considered.
#' @param MatParam Numeric. Matérn smoothness parameter (default = 5/2).
#' @param lengthscale Numeric vector. Lengthscales used for scaling the input space.
#' @param seed Integer. Random seed.
#'
#' @return List with frequency, offset, and coefficient parameters.
#'
#' @keywords internal
#'
#' @importFrom DiceDesign lhsDesign maximinSA_LHS
#'
initialize_basisfun_fillingRFF <- function(dimension, nFreq, MatParam = 5/2, lengthscale, seed=0) {
  # Check if basisFunctionsUsed is valid
  message("You selected space-filling Fourier Features, in the current implementation, this is only available for Mat\ffffffc3rn (2k+1)/2 kernels.")
  if (!requireNamespace("DiceDesign", quietly = TRUE)) {
    stop("Package 'DiceDesign' could not be used")
  }
  if(is.null(seed)){
    warning("Recall that DiceDesign overrides the local seed...")
  }
  if(is.null(MatParam)){
    MatParam<- 5/2
  }
  X <- lhsDesign(nFreq, dimension, seed=seed)$design
  Xopt <- maximinSA_LHS(X, T0=10, c=0.99, it=2000)$design

  freq <- rosenblatt_transform_multivarStudent(Xopt, dimension = dimension, MatParam=MatParam)
  freq <- rbind(freq, freq)
  offset <- c(rep(0, nFreq), rep(-pi/2, nFreq))
  coef <- rep(1/sqrt(nFreq), 2*nFreq)
  return(list(freq=freq, offset=offset, coef=coef))
}


#' Initialize discrete Fourier features
#'
#' Generates basis using discrete cosine/sine terms for each input dimension.
#'
#' @param dimension Integer. Input (\eqn{[\mathbf{x},\,t]}{[x, t]}) dimension.
#' @param maxOrdert Integer. Maximum frequency in t.
#' @param maxOrderx Integer. Maximum frequency in each x.
#'
#' @return List with frequency, offset, and coefficient parameters.
#'
#' @keywords internal
#'
initialize_basisfun_discreteFF <- function(dimension, maxOrdert, maxOrderx) {
  dim_x <- dimension - 1
  # Init
  ai_temp <- as.matrix(rbind(expand.grid(seq(1, maxOrdert),
                                         seq(0, maxOrderx)),
                             expand.grid(seq(1, maxOrdert),
                                         -seq(1, maxOrderx))))
  if(dimension > 2){
    for(i in seq(dim_x-1)){
      comb <- as.matrix(expand.grid(seq(nrow(ai_temp)),
                                    seq(-maxOrderx, maxOrderx)))
      ai_temp <- cbind(ai_temp[comb[, 1], ], comb[, 2])
    }
  }
  for(i in seq(dimension, 1)){
    ai_temp <- ai_temp[order(ai_temp[, i]), ]
  }
  ai_temp <- ai_temp*pi*2
  colnames(ai_temp)<- c("t", paste0("x", seq(dim_x)))
  nFreq <- 2*nrow(ai_temp)
  freq <- rbind(ai_temp, ai_temp)
  offset <- c(rep(0, nFreq/2), rep(-pi/2, nFreq/2))
  coef <- rep(1/sqrt(nFreq), nFreq)

  return(list(freq=freq, offset=offset, coef=coef))
}

#' Evaluate basis functions at given locations.
#'
#' Evaluates all basis functions defined by a parameter list at new locations.
#'
#' @param parameters List of basis function parameters.
#' @param X Matrix or dataframe of evaluation locations.
#' @param lengthscale Numeric vector. Lengthscales used for scaling the input space.
#'
#' @return A matrix of basis function values.
#'
#' @keywords internal
#'
evaluate_basis_functions <- function(parameters, X, lengthscale) {
  type <- parameters$basisFunctionsUsed
  if(type=="inducing points"){
    X <- t(t(X)/lengthscale)
    distances <- crossdist(X, parameters$lengthCoord)
    kernel<- parameters$kernel
    if(kernel=="Exp"){
      kxX <- exp(-distances)
    }
    if(kernel=="Mat32"){
      kxX <- (1+sqrt(3)*distances)*exp(-sqrt(3)*distances)
    }
    if(kernel=="Mat52"){
      kxX <- (1+sqrt(5)*distances+5*distances^2/3)*exp(-sqrt(5)*distances)
    }
    if(kernel=="Gaussian"){
      kxX <- exp(-distances^2/2)
    }
    functions <- kxX %*% t(parameters$sqrtInv_covMat)
    # GP <- functions %*% epsilon
    return(functions)
  }else{
    if(is.data.frame(X)){
      X <- as.matrix(X)
    }
    X <- t(t(X)/lengthscale)
    Xfreq <- X%*%t(parameters$freq)
    basis_fun <- cos(t(Xfreq) + parameters$offset)
    basis_fun <- t(parameters$coef*basis_fun)
    return(basis_fun)
  }
}

