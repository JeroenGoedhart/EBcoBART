#' Exemplary Data Set
#'
#' Contains training data and test data to predict 2 year progression free
#' survival (yes/no) #' based on four types of variables: copy number variation,
#' point mutations, translocations, #' and clinical. For the variables,
#' auxiliary information (co-data) is available which may be used to give more
#' weight to certain variables in the prediction model. This data set is used
#' in the manuscript "Co-data Learning for Bayesian Additive Regression Trees"
#'
#'
#' @format A list object with five data sets:
#' \describe{
#'   \item{Xtrain}{Dataframe with 101 rows (samples) and 140 columns (variables.
#'   Explanatory variables used for fitting BART.
#'   Variable names are anonymized.}
#'   \item{Ytrain}{Numeric of length 101. Binary training response
#'   (0: 2 year progression free survival, 1: disease came back within 2 years)}
#'   \item{Xtest}{Dataframe with 83 rows (samples) and 140 columns (variables).
#'   Explanatory variables used for fitting BART.
#'   Variable names are anonymized.}
#'   \item{Ytest}{Numeric of length 83 Binary training response
#'   (0: 2 year progression free survival, 1: disease came back within 2 years)}
#'   \item{CoData}{Dataframe with 140 rows and 2 columns. Auxiliary information
#'   on the 140 variables. Contains a grouping structure indicating which type
#'   a variable is (copy number variation (CNV), mutation, translocation, or
#'   clinical), and p values (logit scale) for each variable obtained from a
#'   previous study}
#'}
#'
#' @usage data(dat)
#'
#' @references
#' Jeroen M. Goedhart, Thomas Klausch, Jurriaan Janssen, Mark A. van de Wiel.
#' "Co-data Learning for Bayesian Additive Regression Trees."
#' arXiv preprint arXiv:2311.09997. 2023 Nov 16.
#'
#' @author Jeroen M. Goedhart, \email{j.m.goedhart@@amsterdamumc.nl}
#'
#' Jurriaan Janssen
"dat"
