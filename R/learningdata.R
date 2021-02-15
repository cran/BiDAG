#' Asia dataset
#'
#' A synthetic dataset from Lauritzen and Spiegelhalter (1988) about lung
#' diseases (tuberculosis, lung cancer or bronchitis) and visits to Asia.
#'
#' @source \url{https://www.bnlearn.com/bnrepository/}
#' @format A data frame with 5000 rows and 8 binary variables:
#' \itemize{
#' \item D (dyspnoea), binary 1/0 corresponding to "yes" and "no"
#' \item T (tuberculosis), binary 1/0 corresponding to "yes" and "no"
#' \item L (lung cancer), binary 1/0 corresponding to "yes" and "no"
#' \item B (bronchitis), binary 1/0 corresponding to "yes" and "no"
#' \item A (visit to Asia), binary 1/0 corresponding to "yes" and "no"
#' \item S (smoking), binary 1/0 corresponding to "yes" and "no"
#' \item X (chest X-ray), binary 1/0 corresponding to "yes" and "no"
#' \item E (tuberculosis versus lung cancer/bronchitis), binary 1/0 corresponding to "yes" and "no"
#' }
#'@references Lauritzen S, Spiegelhalter D (1988). `Local Computation with Probabilities on Graphical Structures and their Application to Expert Systems (with discussion)'.
#'Journal of the Royal Statistical Society: Series B 50, 157-224.
#'
"Asia"

#' Boston housing data
#'
#' A dataset containing information collected by the U.S Census Service concerning housing
#' in the area of Boston, originally published by Harrison and Rubinfeld (1978).
#'
#' @source \url{http://lib.stat.cmu.edu/datasets/boston}
#' @format A data frame with 506 rows and 14 variables:
#' \itemize{
#' \item CRIM - per capita crime rate by town
#' \item ZN - proportion of residential land zoned for lots over 25,000 sq.ft.
#' \item INDUS - proportion of non-retail business acres per town.
#' \item CHAS - Charles River dummy variable (1 if tract bounds river; 0 otherwise)
#' \item NOX - nitric oxides concentration (parts per 10 million)
#' \item RM - average number of rooms per dwelling
#' \item AGE - proportion of owner-occupied units built prior to 1940
#' \item DIS - weighted distances to five Boston employment centres
#' \item TAX - full-value property-tax rate per $10,000
#' \item RAD - index of accessibility to radial highways
#' \item PTRATIO - pupil-teacher ratio by town
#' \item B - 1000(Bk - 0.63)^2 where Bk is the proportion of blacks by town
#' \item LSTAT - percentage lower status of the population
#' \item MEDV - Median value of owner-occupied homes in $1000's
#' }
#'
#'@references Harrison, D and Rubinfeld, DL (1978)
#' `Hedonic prices and the demand for clean air',
#' Journal of Environmental Economics and Management 5, 81-102.
#'
"Boston"

#' A simulated data set from a Gaussian continuous Bayesian network
#'
#' A synthetic dataset containing 1000 observations generated from a random DAG with 100 continuous nodes.
#'
#' @format A data frame with 1000 rows representing observations of 100 continuous variables: V1, ..., V100
#'
"gsim"

#' A simulated data set from a Gaussian continuous Bayesian network
#'
#' A synthetic dataset containing 100 observations generated from a random DAG with 100 continuous nodes.
#'
#' @format A data frame with 100 rows representing observations of 100 continuous variables: V1, ..., V100
#'
"gsim100"




