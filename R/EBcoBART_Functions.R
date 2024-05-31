#' Convenience function to correctly specify co-data matrix if X contains factor variables.
#'
#' The R package dbarts uses dummy encoding for factor variables so
#' the co-data matrix should contain co-data information for each dummy. If co-data
#' is only available for the factor as a whole (e.g. factor belongs to a group),
#' use this function #' to set-up the co-data in the right-format
#' for the EBcoBART function
#'
#' @param X Explanatory variables. Should be a data.frame. The function is only
#' useful when X contains factor variables.
#' @param CoData The co-data model matrix with co-data information on explanatory variables in X.
#' Should be a matrix, so not a data.frame.
#' If grouping information is present, please encode this yourself using dummies
#' with dummies representing which group a explanatory variable belongs to.
#' The number of rows of the co-data matrix should equal the number of columns of X
#'
#' @return A list object with X: the explanatory variables with factors encoded as dummies
#' and CoData: the co-data matrix with now co-data for all dummies.
#' @export
#'
#' @examples p <- 15
#' n <- 30
#' X <- matrix(runif(n*p),nrow = n, ncol = p) #all continuous variables
#' Fact <- factor(sample(1:3,n,replace = T)) # factor variables
#' X <- cbind.data.frame(X,Fact)

#' G <- 4   #number of groups for co-data
#' CoDat <- rep(1:G, rep(ncol(X)/G,G)) # first 4 covariates in group 1, 2nd 4 covariates in group 2, etc..
#' CoDat <- data.frame(factor(CoDat))
#' CoDat <- model.matrix(~0+., CoDat) # encode the grouping structure with dummies
#' Dat <- Dat_EBcoBART(X = X, CoData = CoDat) #
#' X <- Dat$X
#' CoData <- Dat$CoData
#'
#'@author Jeroen M. Goedhart, \email{j.m.goedhart@@amsterdamumc.nl}
#'
Dat_EBcoBART <- function(X,CoData){

  ## ---------------------------------------------------------------------
  ## Convenience function to recreate co-data matrix if factor variables
  ## in X are present. dbarts uses dummy encoding for factor variables so
  ## CoData matrix should have co-data for each dummy. This function does
  ## that for you.
  ## ---------------------------------------------------------------------

  ## Arguments ##
  # X: covariate model matrix, should be a data.frame.
  #
  # CoData: The codata model matrix with co-data information on the p covariates.
  # If grouping information is present, please encode this yourself using dummies
  # The number of rows of the co-data matrix should equal the number of columns of X

  if (ncol(X) == 0 | nrow(X) == 0){stop("X not specified")}
  if (ncol(CoData) == 0 | nrow(CoData) == 0){stop("CoData  not specified")}
  if (class(X) != "data.frame"){stop("X should be specified as data.frame)")}

  if(ncol(X) != nrow(CoData)){stop("number of columns of X should equal number of rows of CoData")}

  idVars <- which(sapply(X, is.factor))
  replication_times <- rep(1,nrow(CoData))

  if (length(idVars) > 0){
    for (i in 1:length(idVars)) {
      id <- idVars[i]
      reps <- length(unique(X[,id]))
      replication_times[id] <- reps
      remove(id,reps)
    }
  }
  CoDat <- CoData[rep(seq(1:nrow(CoData)), times = replication_times), ]
  X <- model.matrix(~ . + 0, X)
  res <- list(X = X, CoData = CoDat)
  return(res)
}

#' Learning prior covariate weights for BART models using empirical Bayes and co-data.
#'
#' Function that estimates the prior probabilities of variables getting selected in the splitting rules
#' of Bayesian Additive Regression Trees (BART). Estimation is performed using empirical Bayes and co-data,
#' i.e. external information on the explanatory variables.
#'
#' @param Y Response variable that can be either continuous or binary. Should be a numeric.
#' @param X Explanatory variables. Should be a matrix. If X is a data.frame and contains
#' factors, you may consider the function Dat_EBcoBART
#' @param CoData The co-data model matrix with co-data information on explanatory variables in X.
#' Should be a matrix, so not a data.frame.
#' If grouping information is present, please encode this yourself using dummies
#' with dummies representing which group a explanatory variable belongs to.
#' The number of rows of the co-data matrix should equal the number of columns of X
#' @param model What type of response variable Y. Can be either continuous or binary
#' @param nIter Number of iterations of the EM algorithm
#' @param EB Logical (T/F). If true the EM algorithm also estimates prior parameters
#' alpha (of tree structure prior) and k (of leaf node parameter prior). Defaults to False.
#' Setting to true increases computational time.
#' @param Prob_Init Initial vector of splitting probabilities for explanatory variables X.
#' Lenght should equal number of columns of X (and number of rows in CoData).
#' Defaults to 1/p, i.e. equal weight for each variable.
#' @param Info Logical. Asks whether information about the algorithm progress should be printed.
#' Defaults to FALSE.
#' @param Seed Logical asking whether a seed should be set for reproduciblity. Defaults to TRUE.
#' @param ndpost Number of posterior samples returned by dbarts after burn-in. Same as in dbarts.
#' Defaults to 5000.
#' @param nskip Number of burn-in samples. Same as in dbarts. Defaults to 5000.
#' @param nchain Number of independent mcmc chains. Same as in dbarts. Defaults to 5.
#' @param keepevery Thinning. Same as in dbarts. Defaults to 1.
#' @param ntree Number of trees in the BART model. Same as in dbarts. Defaults to 50.
#' @param alpha Alpha parameter of tree structure prior. Called base in dbarts. Defaults to 0.95.
#' @param beta Beta parameter of tree structure prior. Called power in dbarts. Defaults to 2.
#' @param k Parameter for leaf node parameter prior. Same as in dbarts. Defaults to 2.
#' @param sigest Only for continuous response. Estimate of error variance
#' used to set inverse gamma prior on error variance. Same as in dbarts. Defaults to 0.667*var(Y)
#' @param sigdf Only for continuous response. Degrees of freedom for error variance prior
#' Same as in dbarts. Defaults to 10.
#' @param sigquant Only for continuous response. Quantile at which sigest is placed Same as in dbarts.
#' Defaults to 0.75.
#'
#' @return A list object with the estimated variable weights, i.e the probabilities that variables are selected
#' in the splitting rules. Additionaly, the final co-data model is returned. If EB is set to TRUE, estimates of k and alpha
#' are also returned.
#' The prior parameter estimates can then be used in your favorite BART R package that supports
#' manually setting the splitting variable probability vector (dbarts and BARTMachine).
#' @export
#'
#' @examples
#' ###################################
#' ### Continuous response example ###
#' ###################################
#'
#' # Simulate data from Friedman function and define Co-Data as grouping structure
#'
#' sigma <- 1.0
#' N <- 100
#' p <- 500
#' G <- 5   #number of groups
#' CoDat = rep(1:G, rep(p/G,G)) # specify grouping structure
#' CoDat = data.frame(factor(CoDat))
#' CoDat <- model.matrix(~., CoDat) # encode grouping structure by dummies yourself (include intercept)
#' colnames(CoDat)  = paste0("Group ",1:G)

#' g <- function(x) {
#'  10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,101] - 0.5)^2 + 10 * x[,102] + 10 * x[,3]
#' }

#' X <- matrix(runif(N * p), N, p)#'
#' Y <- g(X)+ rnorm(N, 0, sigma)
#'
#' Fit <- EBcoBART(Y=Y,X=X,CoData = CoDat, nIter = 15, model = "continuous",
#'                 EB = F, Info = T, Seed = T,
#'                 nchain = 5, nskip = 1000, ndpost = 1000,
#'                 Prob_Init = rep(1/ncol(X),ncol(X)),
#'                 k = 2, alpha = .95, beta = 2)
#' EstProbs <- Fit$SplittingProbs #estimated prior probabilities of variables getting selected in splitting rules
#'
#' # The prior parameter estimate EstProbs can then be used
#' # in your favorite BART fitting package
#' # We use dbarts:
#'
#' FinalFit <- dbarts::bart(x.train = X, y.train = Y, # training data
#'                         ndpost = 5000, # number of posterior samples
#'                         nskip = 5000,  # number of "warmup" samples to discard
#'                         nchain = 5,    # number of independent chains
#'                         ntree = 50,    # number of trees
#'                         k = 2, base = .95, power = beta, # prior parameters tree
#'                         sigest = .667*var(Y),sigdf = 10, sigquant = .75,  # prior parameters error variance
#'                         splitprobs = EstProbs,   # prior variable weights
#'                         combinechains = T, verbose = F)
#'
#' @examples
#' ###################################
#' ### Binary response example ######
#' ###################################
#'
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
#' CoDat <- model.matrix(~., CoDat) # encode grouping structure by dummies yourself (include intercept)
#' remove(dat)
#'
#' Fit <- EBcoBART(Y=Ytr,X=Xtr,CoData = CoDat, nIter = 15, model = "binary",
#'                 EB = T, Info = T, Seed = T,
#'                 nchain = 5, nskip = 1000, ndpost = 1000,
#'                 Prob_Init = rep(1/ncol(Xtr),ncol(Xtr)),
#'                 k = 2, alpha = .95, beta = 2)
#' EstProbs <- Fit$SplittingProbs #estimated prior probabilities of variables getting selected in splitting rules
#' alpha_EB <- Fit$alpha_est
#' k_EB <- Fit$k_est
#'
#' # The prior parameter estimates EstProbs, alpha_EB,
#' # and k_EB can then be used in your favorite BART fitting package
#' # We use dbarts:
#'
#' FinalFit <- dbarts::bart(x.train = Xtr, y.train = Ytr, # training data
#'                          x.test = Xte, # test X for predictions
#'                          ndpost = 5000,   # number of posterior samples
#'                          nskip = 5000, # number of "warmup" samples to discard
#'                          nchain = 5,   # number of independent, parallel chains
#'                          ntree = 50,    # number of trees
#'                          k = k_EB, base = alpha_EB, power = 2, # prior parameters tree
#'                          splitprobs = EstProbs, # prior variable weights
#'                          combinechains = T, verbose = F)
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

EBcoBART <- function(Y,X,CoData, model,
                     nIter = 10, EB = F,
                     Prob_Init = c(rep(1/ncol(X),ncol(X))),
                     Info = F, Seed = T,
                     ndpost = 5000, # number of posterior samples after burn-in period
                     nskip = 5000, # number of burn-in samples
                     nchain = 5, # number of independent mcmc chains
                     keepevery = 1, # optional thinning
                     ntree = 50, # number of trees in the BART model
                     alpha = .95, beta = 2, k = 2, #Tree hyperparameters
                     sigest = sd(Y)*0.667, sigdf = 10, sigquant = .75 #Error variance hyperparameters, only required for continuous outcome
){
  ## ---------------------------------------------------------------------
  ## Co-data guided Empirical Bayes estimates of splitting probabilities of
  ## BART model
  ## ---------------------------------------------------------------------

  ## Arguments ##
  # Y: response should be a numeric. For a binary model, specify Y as numeric with only zero's and ones

  # X: covariate model matrix, should be a matrix object. If X contains factor variables,
  # please consider Dat_EBcoBART function to set X and CoData matrix in the right format
  # dbarts uses dummy encoding for factors(with dummies for each category and no intercept).
  # This means that co-data is required for each dummy variable.
  # Dat_EBcoBART does that for you and also sets X in the right format for dbarts.

  # CoData: The codata model matrix with co-data information on the p covariates.
  # Should be a matrix object.
  # The number of rows of the co-data matrix should equal the number of columns of X.
  # Again, if X contains factor variables than co-data is required for each unique
  # category of the categorical variable.

  # model: either continuous or binary (0,1)
  # nIter: number of iterations of EM algorithm. Usually, 10 is sufficient (check WAIC trajectory)

  # EB: logical, asks whether hyperparameters k and alpha (base in dbarts package) should be updated.
  # defaults to FALSE (fixed k and alpha). True increases computational time.

  # Prob_Init: initial covariate splitting probabilities, i.e the
  # probability that covariate is considered in the splitting rules
  # This parameter will be estimated

  # Info: print information on fitting process
  # seed: fix seed for reproducibility
  # additional parameters required for fitting a BART model using dbarts. See dbarts documentation for
  # meaning of these parameters. Note that our alpha and beta correspond to base and beta in dbarts, respectively



  # load required R packages
  #if (!suppressMessages(require(dbarts, quietly = T))) {stop("Package dbarts not installed")}
  #if (!suppressMessages(require(loo, quietly = T))) {stop("Package loo not installed")}

  # control statements

  if(!(model == "continuous" | model == "binary")){stop("model should be specified as continuous or binary")}
  if (class(EB) != "logical"){stop("EB is not logical, specify as either TRUE or FALSE")}
  if (ncol(X) == 0 | nrow(X) == 0){stop("X not specified")}
  if (ncol(CoData) == 0 | nrow(CoData) == 0){stop("CoData  not specified")}
  if (length(Y) == 0){stop("Y vector is empty")}

  if(!(class(Y) == "numeric")) {stop("Y is not a numeric. If Y is binary please specify it as numeric vector coded with 0 and 1")}
  if (class(X)[1] != "matrix"){stop("X should be specified as matrix or a data.frame)")}
  if(class(CoData)[1] != "matrix" ){stop("CoData should be specified as a matrix. Please encode dummies yourself")}


  if (model=="continuous" & length(unique(Y)) < 3){stop("Y has less than 3 distinct values while model = continuous is specified")}
  if (model=="binary" & !all(Y==1 | Y== 0)){stop("Binary model, specify binary response as numeric coded with 0 and 1")}
  if(ncol(X) != nrow(CoData)){stop("number of columns of X should equal number of rows of CoData")}
  if(!all(Prob_Init > 0 & Prob_Init < 1)){stop("All prior splitting probabilities in Prob_Init should be between 0 and 1")}
  if(sum(Prob_Init) != 1){stop("Sum of Prob_Init should equal 1")}


  if(nchain<3){stop("Use at least 3 independent chains")}
  if(!all(c(alpha,beta,k,nchain,ndpost,nskip,nIter,keepevery,ntree)>0)){stop("Check if input for bart are all positive numerics")}

  # Initialization
  p <- ncol(X)
  probs <- Prob_Init # initial probabilities a covariate gets selected in the splitting rules
  CoData <- data.frame(CoData) #required for glm fit

  # storage containers
  EstimatedProbs <- matrix(NA, nrow = nIter+1, ncol = ncol(X))
  Codatamodels <- vector("list", length = nIter)
  EstimatedProbs[1,] <- probs
  row.names(EstimatedProbs) <- c(paste("Iter",0:(nIter), sep = " "))
  WAICVector <- c()
  WAIC_Old <- 10e8

  if (EB == T){
    k_Update <- c()
    alpha_Update <- c()
  }

  for (i in 1:nIter) {

    if (Info==T){
      print(paste("EM iteration",i,sep = " "))
    }

    ### step 1: Fit BART model ###
    ##############################

    if(model == "continuous"){

      if(Seed){
        set.seed(4*i^2+202+3*i)
      }

      fit <- dbarts::bart(x.train = X, y.train = Y, # training data
                          ndpost = ndpost,   # number of posterior samples
                          nskip = nskip, # number of "warmup" samples to discard
                          nchain = nchain,   # number of independent, parallel chains
                          keepevery = keepevery, # thinning
                          ntree = ntree,    # number of trees per chain
                          keeptrees = EB,
                          verbose = F,
                          k = k, base = alpha, power = beta, # hyperparameters tree
                          sigest = sigest,sigdf = sigdf, sigquant = sigquant,  # hyperparameters error variance
                          splitprobs = probs,
                          combinechains = T)# hyperparameter that will be updated using EB and co-data

      ## MCMC Convergence check ##
      if (i==1){

        if (Info==T){
          print("Check convergence of mcmc chains")
        }

        samps=fit$sigma
        samps<-matrix(samps, nrow=ndpost,ncol = nchain, byrow = T)
        Rhat <- .MCMC_convergence(samps)

        if (Info==T){
          print(paste("Rhat equals: ",Rhat))
        }

        remove(samps)
        if(Rhat<1.1){
          if (Info==T){
            print("convergence okay")
          }
        } else {
          stop("Not converged yet, please change mcmc sampling settings")
        }
      }


      ## Estimate WAIC ##
      Ypred = fit$yhat.train
      LogLikMatrix = .LikelihoodCont(Ypred = Ypred, Y = Y, sigma = fit$sigma)
      WAICVector[i] <- suppressWarnings(loo::waic(LogLikMatrix)$estimates[3,1])

      if (Info==T){
        print(paste("WAIC equals: ",WAICVector[i]))
      }
    }

    if(model=="binary"){

      if(Seed){
        set.seed(4*i^2+202+3*i)
      }

      fit <- dbarts::bart(x.train = X, y.train = Y,
                          #x.test = Xtest,
                          ndpost = ndpost,                   # number of posterior samples
                          nskip = nskip,                    # number of "warmup" samples to discard
                          nchain = nchain,                      # number of independent, parallel chains
                          keepevery = keepevery, # thinning
                          ntree = ntree,                       # number of trees per chain
                          verbose = F,
                          usequants = F,
                          k = k, base = alpha, power = beta, # hyperparameters tree
                          splitprobs = probs,                # prob that variable is chosen for split
                          keeptrees = EB,             #set to True if updating alpha and k and to False if not
                          combinechains = T)

      ## MCMC Convergence check
      if (i==1){

        if (Info==T){
          print("Check convergence of mcmc chains")
        }

        samps=fit$yhat.train[,sample(1:length(Y),1)]
        samps<-matrix(samps, nrow=ndpost,ncol = nchain, byrow = T)
        Rhat <- .MCMC_convergence(samps)

        if (Info==T){
          print(paste("Rhat equals: ",Rhat))
        }

        remove(samps)
        if(Rhat<1.1){
          if (Info==T){
            print("convergence okay")
          }
        } else {
          stop("Not converged yet, please change mcmc sampling settings")
        }
      }

      ## Estimate WAIC
      Ypred = pnorm(fit$yhat.train)
      Ypred[which(Ypred==0)] <- .0000000000000001
      Ypred[which(Ypred==1)] <- .9999999999999999
      LogLikMatrix = .LikelihoodBin(Ypred = Ypred, Y = Y)
      WAICVector[i] <- suppressWarnings(loo::waic(LogLikMatrix)$estimates[3,1])

      if (Info==T){
        print(paste("WAIC equals: ",WAICVector[i]))
      }
    }
    #### convergence check of EM algorithm


    if (WAICVector[i]>WAIC_Old) {
      EstProb = EstimatedProbs[i-1,]
      EstWAIC = WAIC_Old
      CodataModel <-  Codatamodels[[i-1]]
      if (EB == T) {
        Estk = k_Update[i-1]
        EstAlpha = alpha_Update[i-1]
      }
      break
    } else {
      WAIC_Old <- WAICVector[i]
    }




    # obtain average number of times each variable occurs in the splitting rules
    VarsUsed <- colSums(fit$varcount)  # count of each variable occuring in the splitting rules
    VarsUsed <- VarsUsed/sum(VarsUsed) # normalize count of each variable to probabilities = pure EB updates of hyperparameter S

    ### STEP 2: Fit co-data model ###
    coDataModel <- stats::glm(VarsUsed ~.-1,
                       data=CoData,family=quasibinomial) # the model

    Codatamodels[[i]] <- coDataModel
    probs <- predict(coDataModel, type="response", newdata = CoData) # estimating the co-data moderated estimates of hyperparameter S
    probs[is.na(probs)] <- 0
    probs <-unname(probs)

    EstimatedProbs[i+1,] <- probs

    ## Optional step: update other hyperparameters (alpha and k) of BART using Empirical Bayes ##
    if (EB == T) {
      trees <- dbarts::extract(fit,"trees",chainNum = c(1:nchain), sampleNum=c(base::sample(1:ndpost,0.25*ndpost,replace = F))) # tree structures, for computation, we only select a part of the tree
      #trees <- extract(fit,"trees") # tree structures, for computation, we only select a part of the tree

      # Update leaf node parameter k
      k <- .EstimateLeafNode(Trees = trees, ntree = ntree, model = model)[2]
      k_Update[i] <- k

      # Update tree structure parameter alpha, we keep beta fixed
      trees <- trees[c("tree", "sample", "chain", "n","var","value" )]
      trees$depth <- unname(unlist(by(trees, trees[,c("tree", "sample", "chain")], .getDepth)))
      alpha <- stats::optim(alpha,.LikelihoodTreeStructure, beta = beta, Trees=trees, method = 'Brent', lower = .00001, upper = .9999999)$par
      alpha_Update[i] <- alpha
      remove(trees)
    }



    if (i==nIter){
      print("EM algorithm not converged yet, consider increasing nIter")
      print("Return estimates at last iteration")
      EstProb = EstimatedProbs[i,]
      EstWAIC = WAICVector[i]
      CodataModel <-  Codatamodels[[i]]
      if (EB == T) {
        Estk = k_Update[1]
        EstAlpha = alpha_Update[i]
      }
    }

  }
  # collect results
  if (EB == T){
    res <- list(SplittingProbs = EstProb, Codatamodel= CodataModel, k_est = Estk,alpha_est = EstAlpha)
  } else {
    res <- list(SplittingProbs = EstProb, Codatamodels= CodataModel)
  }
  return(res)
}






###################################
####### Auxiliary Functions #######
###################################
.FiniteSum <- function(x) {
  base::sum(x[is.finite(x)])
}

.MCMC_convergence <- function(Samples){

  ## ---------------------------------------------------------------------
  ## Compute improved Rhat from Vehtari and Gelman to assess mcmc convergence
  ## for BART samples
  ## ---------------------------------------------------------------------

  if (class(Samples)[1] != "matrix") {stop("Samples should be specified as matrix with nsample rows and nchain columns")}
  if(nrow(Samples)<ncol(Samples)){print("Are you sure Samples is specified by nsample rows and nchain columns")}
  Rhat = posterior::rhat(Samples)
  return(Rhat)
}

.LikelihoodBin <- function(Ypred,Y){

  ## ---------------------------------------------------------------------
  ## Compute likelihood for mcmc samples for binary response
  ## ---------------------------------------------------------------------

  result <- apply(Ypred,1, function(x) Y*log(x)+(1-Y)*log(1-x))
  return(t(result))
}

.LikelihoodCont <- function(Ypred, Y,sigma){

  ## ---------------------------------------------------------------------
  ## Compute likelihood for mcmc samples for continuous response
  ## ---------------------------------------------------------------------

  loglik <- -(0.5*(1/sigma^2))*(base::sweep(Ypred,2,Y)^2)-.5*log(sigma^2)-.5*log(2*pi)
  return(loglik)
}

.getDepth <- function(tree) {

  ## ---------------------------------------------------------------------
  ## Compute detph for all nodes, required for ..LikelihoodTreeStructure
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

  LogLike <- ifelse(Trees$var==-1,log(1-alpha*(1+Trees$depth)^(-beta)),log(alpha*(1+Trees$depth)^(-beta)))
  S <- .FiniteSum(LogLike)
  return(-S)
}

.EstimateLeafNode <- function(Trees, ntree, model) {

  ## ---------------------------------------------------------------------
  ## Estimate leaf node prior parameter k from tree output
  ## ---------------------------------------------------------------------

  if(!(model == "continuous" || model =="binary")){stop("model should be specified as continuous or binary")}

  ids <- which(Trees$var==-1) #check which rows correspond to leaf nodes
  samples <- Trees$value[ids] #obtain samples of leaf nodes
  varhat <- (1/length(samples))*.FiniteSum(samples^2) #maximum likelihood estimate of variance for known mean (equals 0)

  if (model=="continuous"){cnst = 0.5}
  if (model=="binary"){cnst = 3}

  k_hat <- cnst/(sqrt(varhat)*sqrt(ntree))
  return(c(varhat=varhat,k_hat=k_hat))
}

