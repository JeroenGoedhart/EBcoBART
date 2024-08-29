# Short description

R package to estimate prior variable weights (the probabilities that
variables are selected in the splitting rules) for Bayesian Additive
Regression Trees (BART). These prior variable weights are estimated
using empirical Bayes and auxiliary information on the variables (termed
co-data). For all details on the method see:

<https://arxiv.org/abs/2311.09997>

In addition to prior variable weights, the package also provides the
option to estimate other hyperparameters of BART using empirical Bayes.

# Installation

    #install.packages("devtools")
    devtools::install_github("JeroenGoedhart/EBcoBART")
    #install.packages("EBcoBART") (if accepted on CRAN)

# Example

Simulate data from Friedman function (function g) and define grouped
Co-data, i.e. assign each covariate in X to a group. EBcoBART then
estimates group specific prior weights. These estimated weights may then
be used in a BART sampler (e.g. dbarts).

    sigma <- 1.0
    N <- 100
    p <- 500
    G <- 5   # number of groups
    CoDat = rep(1:G, rep(p/G,G)) #specify grouping structure
    CoDat = data.frame(factor(CoDat))
    CoDat <- stats::model.matrix(~., CoDat) #encode groups  by dummies yourself(include intercept)
    colnames(CoDat)  = paste0("Group ",1:G)
    g <- function(x) {
     10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,101] - 0.5)^2 + 10 * x[,102] +
     10 * x[,3]
    }
    X <- matrix(runif(N * p), N, p)
    Y <- g(X)+ rnorm(N, 0, sigma)
    set.seed(4) # for reproducible results
    Fit <- EBcoBART(Y=Y,X=X,CoData = CoDat, nIter = 15, model = "continuous",
                    EB_k = FALSE, EB_alpha = FALSE, EB_sigma = FALSE, #EB estimation
                    verbose = FALSE,
                    nchain = 5, nskip = 1000, ndpost = 1000,
                    Prob_Init = rep(1/ncol(X),ncol(X)), # initial prior covariate weights
                    k = 2, alpha = .95, beta = 2)
                    
    EstProbs <- Fit$SplitProbs # estimated prior weights of variables (group-specific)
    EstProbs[c(1,101,201,301,401)] # check weights for each group
    print(Fit)
    summary(Fit)

The prior parameter estimate EstProbs can then be used in your favorite
BART fitting package. We use dbarts:

    FinalFit <- dbarts::bart(x.train = X, y.train = Y,
                            ndpost = 5000,
                            nskip = 5000,
                            nchain = 5,
                            ntree = 50,
                            k = 2, base = .95, power = 2,
                            sigest = stats::sd(Y)*0.667,
                            sigdf = 10, sigquant = .75,
                            splitprobs = EstProbs,
                            combinechains = TRUE, verbose = FALSE)
