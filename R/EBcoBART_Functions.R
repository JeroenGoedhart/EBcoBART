#' Convenience function to correctly specify co-data matrix if X contains
#' factor variables.
#'
#' The R package dbarts uses dummy encoding for factor variables so
#' the co-data matrix should contain co-data information for each dummy.
#' If co-data #' is only available for the factor as a whole
#' (e.g. factor belongs to a group), #' use this function to set-up the co-data
#' in the right-format #' for the EBcoBART function.
#'
#' @param X Explanatory variables. Should be a data.frame. The function is only
#' useful when X contains factor variables.
#' @param CoData The co-data model matrix with co-data information on
#' explanatory variables in X. Should be a matrix, so not a data.frame.
#' If grouping information is present, please encode this yourself using dummies
#' with dummies representing which group a explanatory variable belongs to.
#' The number of rows of the co-data matrix should equal the number
#' of columns of X.
#'
#' @return A list object with X: the explanatory variables with factors encoded
#' as dummies and CoData: the co-data matrix with now co-data for all dummies.
#' @export
#'
#' @examples p <- 15
#' n <- 30
#' X <- matrix(runif(n*p),nrow = n, ncol = p) #all continuous variables
#' Fact <- factor(sample(1:3,n,replace = TRUE)) # factor variables
#' X <- cbind.data.frame(X,Fact)

#' G <- 4   #number of groups for co-data
#' CoDat <- rep(1:G, rep(ncol(X)/G,G)) # first 4 covariates in group 1,
#' #2nd 4 covariates in group 2, etc..
#' CoDat <- data.frame(factor(CoDat))
#' CoDat <- stats::model.matrix(~0+., CoDat) # encode the grouping structure
#' # with dummies
#' Dat <- Dat_EBcoBART(X = X, CoData = CoDat) #
#' X <- Dat$X
#' CoData <- Dat$CoData
#'
#'@author Jeroen M. Goedhart, \email{j.m.goedhart@@amsterdamumc.nl}
#'
Dat_EBcoBART <- function(X,CoData){

  ## control statements ##
  if (ncol(X) == 0 || nrow(X) == 0){stop("X not specified.")}
  if (ncol(CoData) == 0 || nrow(CoData) == 0){stop("CoData  not specified.")}
  if (!(is.data.frame(X))){stop("X should be specified as data.frame.)")}

  if(ncol(X) != nrow(CoData)){
    stop("number of columns of X should equal number of rows of CoData.")}

  idVars <- which(unlist(lapply(X, is.factor)))
  replication_times <- rep(1,nrow(CoData))

  if (length(idVars) > 0){
    for (i in idVars) {
      reps <- length(unique(X[,i]))
      replication_times[i] <- reps
      remove(reps)
    }
  }
  CoDat <- CoData[rep(seq_len(nrow(CoDat)), times = replication_times), ]
  X <- stats::model.matrix(~ . + 0, X)
  res <- list(X = X, CoData = CoDat)
  return(res)
}

#' Learning prior covariate weights for BART models with empirical Bayes
#' and co-data.
#'
#' Function that estimates the prior probabilities of variables being selected
#' in the splitting rules of Bayesian Additive Regression Trees (BART).
#' Estimation is performed using empirical Bayes and co-data, i.e. external
#' information on the explanatory variables.
#'
#' @param Y Response variable that can be either continuous or binary.
#' Should be a numeric.
#' @param X Explanatory variables. Should be a matrix. If X is a data.frame and
#' contains factors, you may consider the function Dat_EBcoBART
#' @param CoData The co-data model matrix with co-data information on
#' explanatory variables in X. Should be a matrix, so not a data.frame.
#' If grouping information is present, please encode this yourself using dummies
#' with dummies representing which group a explanatory variable belongs to.
#' The number of rows of the co-data matrix should equal the number of columns
#' of X. If no CoData is available, but one aims to estimate either prior para-
#' meter k, alpha or sigma, please specify CoData == NULL.
#' @param model What type of response variable Y. Can be either continuous or
#' binary
#' @param nIter Number of iterations of the EM algorithm
#' @param EB_k Logical (T/F). If true, the EM algorithm also estimates prior
#' parameter k (of leaf node parameter prior). Defaults to False.
#' Setting to true increases computational time.
#' @param EB_alpha Logical (T/F). If true, the EM algorithm also estimates prior
#' parameter alpha (of tree structure prior). Defaults to False.
#' Setting to true increases computational time.
#' @param EB_sigma Logical (T/F). If true, the EM algorithm also estimates prior
#' parameters of the error variance. To do so, the algorithm estimates
#' the degrees of freedom (sigdf) and the quantile (sigest) at which sigquant
#' of the probability mass is placed. Thus, the specified sigquant is kept fixed
#' and sigdf and sigest are updated. Defaults to False.
#' @param Prob_Init Initial vector of splitting probabilities for
#' explanatory variables X. #' Length should equal number of columns of X
#' (and number of rows in CoData).
#' Defaults to 1/p, i.e. equal weight for each variable.
#' @param verbose Logical. Asks whether algorithm progress
#' should be printed. Defaults to FALSE.
#' @param ndpost Number of posterior samples returned by dbarts after burn-in.
#' Same as in dbarts. Defaults to 5000.
#' @param nskip Number of burn-in samples. Same as in dbarts. Defaults to 5000.
#' @param nchain Number of independent mcmc chains. Same as in dbarts.
#' Defaults to 5.
#' @param keepevery Thinning. Same as in dbarts. Defaults to 1.
#' @param ntree Number of trees in the BART model. Same as in dbarts.
#' Defaults to 50.
#' @param alpha Alpha parameter of tree structure prior. Called base in dbarts.
#' Defaults to 0.95.
#' If EB_alpha is TRUE, this parameter will be the starting value.
#' @param beta Beta parameter of tree structure prior. Called power in dbarts.
#' Defaults to 2.
#' @param k Parameter for leaf node parameter prior. Same as in dbarts.
#' Defaults to 2. If EB_k is TRUE, this parameter will be the starting value.
#' @param sigest Only for continuous response. Estimate of error variance
#' used to set scaled inverse Chi^2 prior on error variance. Same as in dbarts.
#' Defaults to 0.667*var(Y). #' If EB_sigma is TRUE, this parameter will be the
#' starting value.
#' @param sigdf Only for continuous response. Degrees of freedom for error
#' variance prior. Same as in dbarts. Defaults to 10. If EB_sigma is TRUE,
#' this parameter will be the starting value.
#' @param sigquant Only for continuous response. Quantile at which sigest is
#' placed Same as in dbarts. Defaults to 0.75. If EB_sigma is TRUE,
#' this parameter will be fixed, only sigdf and sigest will be updated.
#'
#' @return A list object with the estimated variable weights,
#' i.e the probabilities #' that variables are selected in the splitting rules.
#' Additionally, the final co-data model is returned. If EB is set to TRUE,
#' estimates of k and/or alpha and/or (sigdf, sigest) are also returned.
#' The prior parameter estimates can then be used in your favorite BART R
#' package that supports #' manually setting the splitting variable
#' probability vector (dbarts and BARTMachine).
#' @export
#'
#' @examples
#' ###################################
#' ### Binary response example ######
#' ###################################
#' # For continuous response example, see README.
#' # Use data set provided in R package
#' # We set EB=T indicating that we also estimate
#' # tree structure prior parameter alpha
#' # and leaf node prior parameter k
#'
#' data(dat)
#' Xtr <- as.matrix(dat$Xtrain) # Xtr should be matrix object
#' Ytr <- dat$Ytrain
#' Xte <- as.matrix(dat$Xtest) # Xte should be matrix object
#' Yte <- dat$Ytest
#' CoDat <- dat$CoData
#' CoDat <- stats::model.matrix(~., CoDat) # encode grouping by dummies
#' #(include intercept)
#'
#' set.seed(4) # for reproducible results
#' Fit <- EBcoBART(Y = Ytr, X = Xtr, CoData = CoDat,
#'                 nIter = 2,         # Low! Only for illustration
#'                 model = "binary",
#'                 EB_k = TRUE, EB_alpha = TRUE,
#'                 EB_sigma = FALSE,
#'                 verbose = TRUE,
#'                 ntree = 5,         # Low! Only for illustration
#'                 nchain = 3,
#'                 nskip = 500,       # Low! Only for illustration
#'                 ndpost = 500,      # Low! Only for illustration
#'                 Prob_Init = rep(1/ncol(Xtr), ncol(Xtr)),
#'                 k = 2, alpha = .95, beta = 2)
#' EstProbs <- Fit$SplitProbs # estimated prior weights of variables
#' alpha_EB <- Fit$alpha_est
#' k_EB <- Fit$k_est
#' print(Fit)
#' summary(Fit)
#'
#' # The prior parameter estimates EstProbs, alpha_EB,
#' # and k_EB can then be used in your favorite BART fitting package
#' # We use dbarts:
#'
#' FinalFit <- dbarts::bart(x.train = Xtr, y.train = Ytr,
#'                          x.test = Xte,
#'                          ntree = 5,         # Low! Only for illustration
#'                          nchain = 3,        # Low! Only for illustration
#'                          nskip = 200,       # Low! Only for illustration
#'                          ndpost = 200,      # Low! Only for illustration
#'                          k = k_EB, base = alpha_EB, power = 2,
#'                          splitprobs = EstProbs,
#'                          combinechains = TRUE, verbose = FALSE)
#'
#' @references
#' \CRANpkg{dbarts}
#'
#' Jerome H. Friedman.
#' "Multivariate Adaptive Regression Splines."
#' The Annals of Statistics, 19(1) 1-67 March, 1991.
#'
#' Hugh A. Chipman, Edward I. George, Robert E. McCulloch.
#' "BART: Bayesian additive regression trees."
#' The Annals of Applied Statistics, 4(1) 266-298 March 2010.
#'
#' Jeroen M. Goedhart, Thomas Klausch, Jurriaan Janssen, Mark A. van de Wiel.
#' "Co-data Learning for Bayesian Additive Regression Trees."
#' arXiv preprint arXiv:2311.09997. 2023 Nov 16.
#' @author Jeroen M. Goedhart, \email{j.m.goedhart@@amsterdamumc.nl}

EBcoBART <- function(Y,X, model,
                     CoData,
                     nIter = 10,
                     EB_k = FALSE,
                     EB_alpha = FALSE,
                     EB_sigma = FALSE,
                     Prob_Init = c(rep(1 / ncol(X), ncol(X))),
                     verbose = FALSE,
                     ndpost = 5000,
                     nskip = 5000,
                     nchain = 5,
                     keepevery = 1,
                     ntree = 50,
                     alpha = .95, beta = 2, k = 2,
                     sigest = stats::sd(Y)*0.667, sigdf = 10, sigquant = .75
) {

  ## control statements ##
  if(!(model == "continuous" || model == "binary")) {
    stop("model should be specified as continuous or binary.")}
  if (!(is.logical(EB_k))) {
    stop("EB_k is not logical, specify as either TRUE or FALSE.")}
  if (!(is.logical(EB_alpha))) {
    stop("EB_alpha is not logical, specify as either TRUE or FALSE.")}
  if (!(is.logical(EB_sigma))) {
    stop("EB_sigma is not logical, specify as either TRUE or FALSE.")}
  if (model == "binary") {
    EB_sigma <- FALSE} #error variance only for continuous outcome


  if(!(is.numeric(Y))) {stop("Y is not a numeric. If Y is binary please specify
                             it as numeric vector coded with 0 and 1.")}
  if (!(is.matrix(X))) {stop("X should be specified as matrix.")}
  if (ncol(X) < 1 || nrow(X) < 1) {stop("X not specified.")}
  if (length(Y) == 0) {stop("Y vector is empty.")}

  ### CoData control statements ###
  if (!is.null(CoData)) {
    if(ncol(X) != nrow(CoData)) {
      stop("number of columns of X should equal number of rows of CoData.")}

    if(!(is.matrix(CoData))) {
      stop("CoData should be specified as a matrix")}
    if (ncol(CoData) < 1 || nrow(CoData) == 0) {
      stop("CoData  not specified")}
    if(sum(is.na(CoData)) > 0){
      stop("CoData has missing values")}
    CoDatFlag <- TRUE
  }

  if (is.null(CoData)) {
    message("CoData is NULL. Prior covariate weights are not estimated")
    if(!any(c(EB_k, EB_alpha, EB_sigma))) {
      stop("CoData is NULL and all EB estimators are set to FALSE.")
    }
    CoDatFlag <- FALSE
  }

  if (model == "continuous" & length(unique(Y)) < 3) {
    stop("Y has less than 3 distinct values
         while model = continuous is specified.")}
  if (model == "binary" & !all(Y == 1 | Y == 0)) {
    stop("Binary model, specify binary response as numeric coded with 0 and 1.")}

  if(!all(Prob_Init > 0 & Prob_Init < 1)) {
    stop("All prior splitting probabilities in Prob_Init
         should be between 0 and 1.")}
  if(base::sum(Prob_Init) != 1) {stop("Sum of Prob_Init should equal 1.")}


  if(nchain < 3) {stop("Use at least 3 independent chains")}
  if(!all(c(alpha, beta, k, nchain, ndpost, nskip, nIter, keepevery, ntree) > 0)) {
    stop("Check if input for bart are all positive numerics")
  }

  # Initialization
  p <- ncol(X)
  probs <- Prob_Init  # initial probs var is selected in the splitting rules
  CoData <- data.frame(CoData)  #required for glm fit

  # storage containers
  EstimatedProbs <- matrix(NA, nrow = nIter+1, ncol = ncol(X))
  Codatamodels <- vector("list", length = nIter)
  EstimatedProbs[1,] <- probs
  row.names(EstimatedProbs) <- c(paste("Iter",0:(nIter), sep = " "))
  WAICVector <- c()
  WAIC_Old <- 10e8
  EB <- FALSE

  if (EB_k == TRUE){
    k_Update <- c()
    k_Update[1] <- k
    EB <- TRUE
  }
  if (EB_alpha == TRUE){
    alpha_Update <- c()
    alpha_Update[1] <- alpha
    EB <- TRUE
  }

  if (EB_sigma == TRUE){
    sigdf_Update <- c()
    sigest_Update <- c()
    sigtau_Update <- c() # scale parameter of (scaled) inverse chi^2
    sigdf_Update[1] <- sigdf
    sigest_Update[1] <- sigest
    sigtau_Update[1] <- 0
  }

  for (i in seq_len(nIter)) {

    if (verbose == TRUE){
      cat("EM iteration ",i)
    }

    ### step 1: Fit BART model ###
    ##############################

    if(model == "continuous"){

      fit <- dbarts::bart(x.train = X, y.train = Y,
                          ndpost = ndpost,
                          nskip = nskip,
                          nchain = nchain,
                          keepevery = keepevery,
                          ntree = ntree,
                          keeptrees = EB,
                          verbose = FALSE,
                          k = k, base = alpha, power = beta,
                          sigest = sigest, sigdf = sigdf, sigquant = sigquant,
                          splitprobs = probs,
                          combinechains = TRUE)

      ## Estimate WAIC ##
      Ypred <- fit$yhat.train
      LogLikMatrix <- .LikelihoodCont(Ypred = Ypred, Y = Y, sigma = fit$sigma)
      WAICVector[i] <- suppressWarnings(loo::waic(LogLikMatrix)$estimates[3,1])

      if (verbose == TRUE){
        cat("   WAIC equals: ",WAICVector[i], "\n")
      }



      ## MCMC Convergence check ##
      if (i == 1){

        if (verbose == TRUE){
          cat("\n","Check convergence of mcmc chains")
        }

        samps<-matrix(fit$sigma, nrow = ndpost, ncol = nchain, byrow = TRUE)
        Rhat <- .MCMC_convergence(samps)

        if(Rhat < 1.1){
          if (verbose == TRUE){
            cat(" <- convergence okay", "\n","\n")
          }
        } else {
          stop("MCMC not converged, change mcmc sampling settings.")
        }
      }
      if (EB_sigma == TRUE){
        HypEsts <- .EstSigma(sigma = fit$sigma, quant = sigquant)
        sigdf <- HypEsts[1]
        sigest <- HypEsts[2]
        tau <- HypEsts[3]
        sigdf_Update[i+1] <- sigdf
        sigest_Update[i+1] <- sigest
        sigtau_Update[i+1] <- tau
      }
    }

    if(model == "binary"){

      fit <- dbarts::bart(x.train = X, y.train = Y,
                          ndpost = ndpost,
                          nskip = nskip,
                          nchain = nchain,
                          keepevery = keepevery,
                          ntree = ntree,
                          verbose = FALSE,
                          usequants = FALSE,
                          k = k, base = alpha, power = beta,
                          splitprobs = probs,
                          keeptrees = EB,
                          combinechains = TRUE)
      ## Estimate WAIC
      Ypred <- stats::pnorm(fit$yhat.train)
      Ypred[which(Ypred == 0)] <- .0000000000000001
      Ypred[which(Ypred == 1)] <- .9999999999999999
      LogLikMatrix <- .LikelihoodBin(Ypred = Ypred, Y = Y)
      WAICVector[i] <- suppressWarnings(loo::waic(LogLikMatrix)$estimates[3,1])

      if (verbose == TRUE){
        cat("    WAIC = ", WAICVector[i], "\n")
      }

      ## MCMC Convergence check
      if (i == 1){

        if (verbose == TRUE){
          cat("\n","Check convergence of mcmc chains")
        }

        samps <- fit$yhat.train[,sample(seq_len(length(Y)),1)]
        samps <- matrix(samps, nrow = ndpost,ncol = nchain, byrow = TRUE)
        Rhat <- .MCMC_convergence(samps)
        remove(samps)
        if(Rhat < 1.1){
          if (verbose == TRUE){
            cat(" <- convergence okay", "\n","\n")
          }
        } else {
          stop("Not converged yet, please change mcmc sampling settings.")
        }
      }
    }

    #### convergence check of EM algorithm
    if (WAICVector[i] > WAIC_Old) {
      EstProb <- EstimatedProbs[i-1,]
      EstWAIC <- WAIC_Old
      CodataModel <-  Codatamodels[[i-1]]
      Converged <- TRUE
      Iter <- i-1
      if (EB_k == TRUE) {
        Estk <- k_Update[i-1]
      }
      if (EB_alpha == TRUE) {
        EstAlpha <- alpha_Update[i-1]
      }
      if (EB_sigma == TRUE) {
        Est_sigdf <- sigdf_Update[i-1]
        Est_sigest <- sigest_Update [i-1]
        Est_tau <- sigtau_Update[i-1]
        Est_sigquant <- sigquant
      }

      if (verbose == TRUE){message("Minimum WAIC at iteration ", i-1)}
      break
    } else {
      WAIC_Old <- WAICVector[i]
    }



    if (CoDatFlag == TRUE) {
    # obtain average number of times each variable occurs in the splitting rules
    VarsUsed <- base::colSums(fit$varcount)
    # count of each variable occurring in the splitting rules
    VarsUsed <- VarsUsed / base::sum(VarsUsed)
    # normalize count of each variable to probabilities = pure EB updates
    # of hyperparameter S

    ### STEP 2: Fit co-data model ###
    coDataModel <- stats::glm(VarsUsed ~.-1,
                       data = CoData,family = stats::quasibinomial) # the model

    Codatamodels[[i]] <- coDataModel
    probs <- stats::predict(coDataModel, type = "response", newdata = CoData)
    # estimating the co-data moderated estimates of hyperparameter S
    probs[is.na(probs)] <- 0
    probs <- unname(probs)

    }

    EstimatedProbs[i+1,] <- probs

    ## Optional step: update other hyperparameters of BART using EB ##
    if (EB_k == TRUE || EB_alpha == TRUE) {

      trees <- dbarts::extract(fit, "trees", chainNum = c(1:nchain),
                               sampleNum=c(base::sample(1 : (ndpost / keepevery),
                                                        0.25 * (ndpost/keepevery),
                                                        replace = FALSE)))
      #for computation, we only randomly select 25% of the posterior samples

      # Update leaf node parameter k
      if (EB_k == TRUE){
        k <- .EstimateLeafNode(Trees = trees, ntree = ntree, model = model)[2]
        k_Update[i+1] <- k
      }
      if (EB_alpha == TRUE){
        trees <- trees[c("tree", "sample", "chain", "n","var","value" )]
        trees$depth <- unname(unlist(
          by(trees, trees[,c("tree", "sample", "chain")], .getDepth)))
        alpha <- stats::optim(alpha, .LikelihoodTreeStructure,
                              beta = beta, Trees = trees,
                              method = 'Brent',
                              lower = .00001, upper = .9999999)$par
        alpha_Update[i+1] <- alpha
      }
      remove(trees)
    }



    if (i == nIter){
      warning("EM algorithm not converged yet, consider increasing nIter.
              Return estimates at last iteration.")
      Converged <- FALSE
      EstProb <- EstimatedProbs[i,]
      EstWAIC <- WAICVector[i]
      CodataModel <-  Codatamodels[[i]]
      Iter <- nIter
      if (EB_k == TRUE) {
        Estk <- k_Update[i]
      }
      if (EB_alpha == TRUE){
        EstAlpha <- alpha_Update[i]
      }
      if (EB_sigma == TRUE) {
        Est_sigdf <- sigdf_Update[i]
        Est_sigest <- sigest_Update[i]
        Est_tau <- sigtau_Update[i]
        Est_sigquant <- sigquant
      }
    }

  }
  # collect results
  res <- list()

  if(CoDatFlag == TRUE) {
    res$SplitProbs <- EstProb
    res$CoDataModel <- CodataModel
  }
  if (EB_k == TRUE) {
    res$k_est <- Estk
  }
  if (EB_alpha == TRUE) {
    res$alpha_est <- EstAlpha
  }
  if (EB_sigma == TRUE) {
    res$sigma_est <- c("sigdf" = Est_sigdf, "sigest" = Est_sigest,
                       "sigquant" = Est_sigquant, "sigtau" = Est_tau)
  }
  res$WAIC <- WAICVector
  res$Convergence <- Converged
  res$iteration <- Iter
  class(res) <- "EBcoBART"

  return(res)
}


#####################################
#### Methods for EB-coBART class ####
#####################################



#' @export
print.EBcoBART  <- function(x, ...) {
  WAIC <- x$WAIC
  for (i in seq_along(WAIC)) {
    cat("Iteration ",i,"    WAIC = ", WAIC[i], "\n")

  }
}

#'@export
summary.EBcoBART  <- function(object, ...) {

  if (!object$Convergence){
    cat("EB-coBART not converged. Increase nIter.", "\n")
  }

  if (object$Convergence){
    cat("Convergence okay. Minumum WAIC reached at iteration ",object$iteration, "\n")
    cat("\nEB_coBART estimates:\n")
    if (!is.null(object$SplitProbs)) {
      cat("Prior covariate weights S\n")
    }
    if (!is.null(object$k_est)) {
      cat("Leaf node parameter k =", object$k_est,"\n")
    }
    if (!is.null(object$alpha_est)) {
      cat("Tree structure parameter alpha =", object$alpha_est,"\n")
    }
    if (!is.null(object$sigma_est)) {
      cat("Sigma prior parameters:\n", object$sigma_est,"\n")
    }
    if (!is.null(object$CoDataModel)) {
      cat("\nThe estimated co-data model is:\n\n")
      print(object$CoDataModel$coefficients)
    }
  }
}

###############################################################################
################################### Auxiliary Functions #######################
###############################################################################
.FiniteSum <- function(x) {
  base::sum(x[is.finite(x)])
}

.MCMC_convergence <- function(Samples){

  ## ---------------------------------------------------------------------
  ## Compute improved Rhat from Vehtari and Gelman to assess mcmc convergence
  ## for BART samples
  ## ---------------------------------------------------------------------

  if (!(is.matrix(Samples))) {
    stop("Samples should be specified as matrix with
         nsample rows and nchain columns")}
  if(base::nrow(Samples) < base::ncol(Samples)){
    warning("Are you sure Samples has nsample rows and nchain columns")}
  Rhat <- posterior::rhat(Samples)
  return(Rhat)
}

.LikelihoodBin <- function(Ypred,Y) {

  ## ---------------------------------------------------------------------
  ## Compute likelihood for mcmc samples for binary response
  ## ---------------------------------------------------------------------

  result <- apply(Ypred, 1, function(x) Y * log(x) + (1 - Y) * log(1 - x))
  return(t(result))
}

.LikelihoodCont <- function(Ypred, Y,sigma){

  ## ---------------------------------------------------------------------
  ## Compute likelihood for mcmc samples for continuous response
  ## ---------------------------------------------------------------------

  loglik <- -(0.5 * (1 / sigma ^ 2))*
    (base::sweep(Ypred, 2, Y) ^ 2) - .5 * log(sigma ^ 2)- .5 * log(2 * pi)
  return(loglik)
}

.getDepth <- function(tree) {

  ## ---------------------------------------------------------------------
  ## Compute detph for all nodes, required for .LikelihoodTreeStructure
  ## function
  ## This function is coded by Vincent Dorie (author of dbarts R package)
  ## ---------------------------------------------------------------------

  getDepthRecurse <- function(tree, depth) {
    node <- list(
      depth = depth
    )
    if (tree$var[1] == -1) {
      node$n_nodes <- 1
      return(node)
    }

    headOfLeftBranch <- tree[-1,]
    left <- getDepthRecurse(headOfLeftBranch, depth + 1)
    n_nodes.left <- left$n_nodes
    left$n_nodes <- NULL

    headOfRightBranch <- tree[seq.int(2 + n_nodes.left, nrow(tree)),]
    right <- getDepthRecurse(headOfRightBranch, depth + 1)
    n_nodes.right <- right$n_nodes
    right$n_nodes <- NULL

    node$n_nodes <- 1 + n_nodes.left + n_nodes.right
    node$depth <- c(node$depth, left$depth, right$depth)
    return(node)
  }
  result <- getDepthRecurse(tree, 0)

  return(result$depth)
}

.LikelihoodTreeStructure <- function(alpha,beta, Trees) {

  ## ---------------------------------------------------------------------
  ## likelihood for optimization of tree structure parameter alpha
  ## ---------------------------------------------------------------------

  LogLike <- ifelse(Trees$var == -1,log(1-alpha * (1 + Trees$depth) ^ (-beta)),
                    log(alpha * (1 + Trees$depth) ^ (-beta)))
  S <- .FiniteSum(LogLike)
  return(-S)
}

.EstimateLeafNode <- function(Trees, ntree, model) {

  ## ---------------------------------------------------------------------
  ## Estimate leaf node prior parameter k from tree output
  ## ---------------------------------------------------------------------

  if(!(model == "continuous" || model == "binary")) {
    stop("model should be specified as continuous or binary")
  }

  ids <- which(Trees$var == -1) #check which rows correspond to leaf nodes
  samples <- Trees$value[ids] #obtain samples of leaf nodes
  varhat <- (1 / length(samples)) * .FiniteSum(samples^2) #maximum likelihood
  # estimate of variance for known mean (equals 0)

  if (model == "continuous") {cnst <- 0.5}
  if (model == "binary") {cnst <- 3}

  k_hat <- cnst / (sqrt(varhat) * sqrt(ntree))
  return(c(varhat = varhat,k_hat = k_hat))
}

.EstSigma <- function(sigma, quant) {

  ## -----------------------------------------------------------------------
  ## Estimate error variance parameters df and sigest from posterior samples
  ## -----------------------------------------------------------------------


  HypEsts <- univariateML::mlinvgamma(sigma)
  shape <- unname(HypEsts[1])
  scale <- HypEsts[2]
  nu <- 2 * shape #transform parameters of invgamma to invchi^2
  tau <- scale / shape
  sigest <- extraDistr::qinvchisq(quant, nu = nu, tau = tau)
  return(c("df" = nu, "sigest" = sigest, "tau" = tau))
}
