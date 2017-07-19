#'Initialising score object
#'
#'This function returns an object of class scoreparameters containing the data and parameters needed for calculation of the BDe/BGe score.
#' @param n number of nodes (variables) in the Bayesian network
#' @param data the data matrix with n columns (the number of variables) and a number of rows equal to the number of observations
#' @param scoretype the score to be used to assess the DAG structure: "bde" for binary data, "bge" for Gaussian data
#' @param weightvector (optional) a numerical vector of positive values representing the weight of each observation; should be NULL(default) for non-weighted data 
#' @param bgepar a list which contains parameters for BGe score:
#' \itemize{
#' \item am (optional) a positive numerical value, 1 by default
#' \item aw (optional) a positive numerical value should be more than \code{n+1}, \code{n+am+1} by default
#' }
#' 
#' @param bdepar a list which contains parameters for BDe score:
#' \itemize{
#' \item chi (optional) a positive number of prior pseudo counts used by the BDe score, 0.5 by default
#' \item edgepf (optional) a positive numerical value providing the edge penalization factor to be combined with the BDe score, 2 by default
#' }
#'@return an object of class \code{scoreparameters}, which includes all necessary information for calculating the BDe/BGe score
#'@references Geiger D and Heckerman D (2002). Parameter priors for directed acyclic graphical models and the characterization of several probability distributions. The Annals of Statistics 30, 1412-1440.
#'@references Kuipers J, Moffa G and Heckerman D (2014). Addendum on the scoring of Gaussian acyclic graphical models. The Annals of Statistics 42, 1689-1691.
#'@references Heckerman D and Geiger D (1995). Learning Bayesian networks: A unification for discrete and Gaussian domains. In Eleventh Conference on Uncertainty in Artificial Intelligence, pages 274-284.
#'@references Scutari M (2016). An Empirical-Bayes Score for Discrete Bayesian Networks. Journal of Machine Learning Research 52, 438-448 
#'@examples
#' myDAG<-pcalg::randomDAG(20, prob=0.15, lB = 0.4, uB = 2) 
#' myData<-pcalg::rmvDAG(200, myDAG) 
#' myScore<-scoreparameters(20, "bge", myData)
#'@export
# a constructor function for the "scoreparameters" class
scoreparameters<-function(n, scoretype=c("bge","bde"), data, weightvector=NULL, bgepar=list(am=1, aw=NULL),
                          bdepar=list(edgepf=2, chi=0.5)) {
  if (!(scoretype%in%c("bde","bge"))) {
    stop("Scoretype should be bge (for continious data) or bde (for binary data)")
  }
  if (ncol(data)!=n) {
    stop("n and number of columns in data do not match")
  }
  if (!is.null(weightvector)) {
  if (length(weightvector)!=nrow(data)) {
    stop("length of weightvector does not match number of columns (observations) in data")
  }
  }
  initparam<-list()
  initparam$type=scoretype
  initparam$weightvector<-weightvector
  if(scoretype=="bge") {
    initparam$data<-data
    if (is.null(weightvector)) {
        N<-nrow(data)
        covmat<-cov(data)
        means<-colMeans(data)
    } else {
        N=sum(weightvector)
        forcov<-cov.wt(data,wt=weightvector,cor=TRUE)
        covmat<-forcov$cov
        means<-forcov$center
      }
    if (is.null(bgepar$aw)) {
    bgepar$aw<-n+bgepar$am+1
    }
    mu0<-numeric(n)
    T0scale <- bgepar$am*(bgepar$aw-n-1)/(bgepar$am+1) # This follows from equations (19) and (20) of [GH2002]
    T0<-diag(T0scale,n,n)
    initparam$TN <- T0 + (N-1)* covmat + ((bgepar$am*N)/(bgepar$am+N))* (mu0 - means)%*%t(mu0 - means)
    initparam$awpN<-bgepar$aw+N
    constscorefact<- -(N/2)*log(pi) + (1/2)*log(bgepar$am/(bgepar$am+N))

    initparam$scoreconstvec<-numeric(n)
    for (j in 1:n) {# j represents the number of parents plus 1
      awp<-bgepar$aw-n+j
      initparam$scoreconstvec[j]<-constscorefact - lgamma(awp/2) + lgamma((awp+N)/2) + ((awp+j-1)/2)*log(T0scale)
    }
  }
  else if (scoretype=="bde") {
    if (is.null(weightvector)) {
      initparam$d1<-t(data)
      initparam$d0<-t((1-data))
    } else {
      initparam$d1<-t(data*weightvector)
      initparam$d0<-t((1-data)*weightvector)
    }
    initparam$data<-t(data)
    maxparents<-n-1
    initparam$scoreconstvec<-rep(0,maxparents+1)
    initparam$chi<-bdepar$chi #1
    initparam$pf<-bdepar$edgepf
    for(i in 0:maxparents){ # constant part of the score depending only on the hyperparameters
      noparams<-2^i
      initparam$scoreconstvec[i+1]<-noparams*lgamma(initparam$chi/noparams)-2*noparams*lgamma(initparam$chi/(2*noparams))-i*log(initparam$pf) #
    }
  }  
  
  attr(initparam, "class") <- "scoreparameters"
  initparam
}

