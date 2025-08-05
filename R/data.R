#' Lymphoma
#'
#' Contains training data and test data to predict 2 year progression free
#' survival (yes/no) based on four types of variables: copy number variation,
#' point mutations, translocations, and clinical. For the variables,
#' auxiliary information (co-data) is available which may be used to give more
#' weight to certain variables in the prediction model. This data set is used
#' in the manuscript "Co-data Learning for Bayesian Additive Regression Trees"
#'
#'
#' @format A list object with five objects:
#' \describe{
#'   \item{Xtrain}{Dataframe with 101 rows (samples) and 140 columns (variables).
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
#' @usage data(Lymphoma)
#'
#' @references
#' Jeroen M. Goedhart, Thomas Klausch, Jurriaan Janssen, Mark A. van de Wiel.
#' "Co-data Learning for Bayesian Additive Regression Trees."
#' arXiv preprint arXiv:2311.09997. 2023 Nov 16.
#'
#' @author Jeroen M. Goedhart, \email{j.m.goedhart@@amsterdamumc.nl}
#'
#' Jurriaan Janssen
"Lymphoma"

#' Bloodplatelet
#'
#' Contains not standardized messenger-RNA expression measurements, derived from
#' blood platelets, which are used to classify breast cancer versus non-small-
#' cell lung cancer patients. For the 500 m-RNA variables, co-data is available.
#' Co-data is defined by estimated p-values (- logit scale) of all the 500 m-RNA
#' for three different classification tasks: 1) colorectal cancer vs. control
#' patients, 2) pancreas cancer vs. control patients, and 3) pancreas cancer vs.
#' colorectal cancer. Co-data is therefore informative if different cancer
#' classification tasks have similar important m-RNA variables.
#' See Novianti and others (2017) \doi{doi:10.1093/bioinformatics/btw837} for details
#' on the complete data set, from which this data is derived.
#'
#'
#' @format A list object with five objects:
#' \describe{
#'   \item{Xtrain}{Data frame with 101 rows (samples) and 140 columns (variables).
#'   Explanatory variables used for fitting BART.
#'   Variable names are present.}
#'   \item{Y}{Numeric of length 100. Binary training response
#'   (0: Breast cancer, 1: non-small-cell lung cancer)}
#'   \item{CoData}{Matrix with 500 rows and 4 columns. Auxiliary information
#'   on the 500 variables. Contains, for each variable, estimated p-values
#'   from three different classification tasks. P-values are -logit transformed.
#'   An intercept is included to the co-data matrix.}
#'}
#'
#' @usage data(Bloodplatelet)
#'
#' @references
#' P. W. Novianti, B.C. Snoek, S. Wilting, and M. A. Van De Wiel,
#' Better diagnostic signatures from RNAseq data through use of auxiliary co-data
#' 2017 Bioinformatics, Vol. 33, No. 10, p. 1572-1574
#'
#' @author Jeroen M. Goedhart, \email{j.m.goedhart@@amsterdamumc.nl}
#'
#' Mark A van de Wiel
"Bloodplatelet"
