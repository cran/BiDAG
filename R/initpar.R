#started changing for background nodes


#'Initialising score object
#'
#'This function returns an object of class scoreparameters containing the data and parameters needed for calculation of the BDe/BGe score, or a user defined score.
#' @param n number of nodes (variables) in the Bayesian network (excluding background nodes)
#' @param data the data matrix with n columns (the number of variables) and a number of rows equal to the number of observations
#' @param scoretype the score to be used to assess the DAG structure:
#'  "bge" for Gaussian data, "bde" for binary data, 
#'  "bdecat" for categorical data, "dbn" for dynamic Bayesian networks,
#'  "usr" for a user defined score
#' @param weightvector (optional) a numerical vector of positive values representing the weight of each observation; should be NULL(default) for non-weighted data 
#' @param bgnodes (optional) a numerical vector which contains numbers of columns in the data defining background nodes, background nodes are nodes which have no parents but can be parents of other nodes in the network
#' @param bgepar a list which contains parameters for BGe score:
#' \itemize{
#' \item am (optional) a positive numerical value, 1 by default
#' \item aw (optional) a positive numerical value should be more than \code{n+1}, \code{n+am+1} by default
#' }
#' @param bdepar a list which contains parameters for BDe score for binary data:
#' \itemize{
#' \item chi (optional) a positive number of prior pseudo counts used by the BDe score, 0.5 by default
#' \item edgepf (optional) a positive numerical value providing the edge penalization factor to be combined with the BDe score, 2 by default
#' }
#' @param bdecatpar a list which contains parameters for BDe score for categorical data:
#' \itemize{
#' \item chi (optional) a positive number of prior pseudo counts used by the BDe score, 0.5 by default
#' \item edgepf (optional) a positive numerical value providing the edge penalization factor to be combined with the BDe score, 2 by default
#' }
#' @param dbnpar which type of score to use for the slices
#' \itemize{
#' \item dbnscoretype (optional) "bge" for continuous data, "bde" for binary data or "bdecat" for categorical
#' \item slices the number of time slices
#' }
#' @param usrpar a list which contains parameters for the user defined score
#' \itemize{
#' \item pctesttype (optional) "bde" to use the independence test for binary data, "bge" for continuous data
#' }
#'@param edgepmat (optional) a matrix of positive numerical values providing the per edge penalization factor to be added to the score, NULL by default
#'@param nodeslabels (optional) a vector of characters which denote the names of nodes in the Bayesian network
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
scoreparameters<-function(n, scoretype=c("bge","bde","bdecat","dbn","usr"), data, weightvector=NULL,
                          bgnodes=NULL,
                          bgepar=list(am=1, aw=NULL),
                          bdepar=list(chi=0.5, edgepf=2),
                          bdecatpar=list(chi=0.5, edgepf=2),
                          dbnpar=list(dbnscoretype=c("bge","bde","bdecat"),slices=2), 
                          usrpar=list(pctesttype=c("bge","bde","bdecat")), 
                          edgepmat=NULL, nodeslabels=NULL) {
  
  if (!(scoretype%in%c("bge", "bde", "bdecat", "dbn","usr"))) {
    stop("Scoretype should be bge (for continuous data), bde (for binary data) bdecat (for categorical data) or usr (for user defined)")
  }
  
  if (anyNA(data)) {
    stop("Dataset contains missing data")  
  }
  
  if (ncol(data)!=n) {
    stop("n and number of columns in data do not match")
  }
  
  if (!is.null(weightvector)) {
    if (length(weightvector)!=nrow(data)) {
      stop("Length of the weightvector does not match the number of columns (observations) in data")
    }
  }
  
  if (scoretype=="bde") {
    if (!all(sapply(data,function(x)x%in%c(0,1)))) {
      stop("Dataset contains non-binary values")  
    }
  }
  
  if (scoretype=="bdecat") { # convert factors to levels  
    indx <- sapply(data, is.factor)
    data[indx] <- lapply(data[indx], function(x) as.numeric(x)-1)
    if (!all(unlist(lapply(data, function(x) setequal(unique(x),c(0:max(x))))))) {
      stop("Some variable levels are not present in the data")  
    }    
  }
  
  if (is.null(nodeslabels)) {
    if(all(is.character(colnames(data)))){
      nodeslabels<-colnames(data)
    } else {
      nodeslabels<-sapply(c(1:n), function(x)paste("v",x,sep=""))
    }
  }
  
  colnames(data)<-nodeslabels
  
  initparam<-list()
  initparam$labels<-nodeslabels
  initparam$type<-scoretype
  initparam$weightvector<-weightvector
  initparam$data<-data
  initparam$bgnodes<-bgnodes
  initparam$bgn<-length(bgnodes)
  if(scoretype=="dbn") {
    initparam$nsmall<-n/dbnpar$slices  
  } else {
    initparam$nsmall<-n-initparam$bgn
  }
  
  if (is.null(bdepar$edgepmat)) {
    initparam$logedgepmat <- NULL
  } else {
    initparam$logedgepmat <- log(bdepar$edgepmat)
  }
  
  if(scoretype=="bge") {
    if (is.null(weightvector)) {
      N<-nrow(data)
      covmat<-cov(data)*(N-1)
      means<-colMeans(data)
    } else {
      N<-sum(weightvector)
      forcov<-cov.wt(data,wt=weightvector,cor=TRUE,method="ML")
      covmat<-forcov$cov*N
      means<-forcov$center
    }
    if (is.null(bgepar$aw)) {
      bgepar$aw<-n+bgepar$am+1
    }
    
    initparam$am <- bgepar$am # store parameters
    initparam$aw <- bgepar$aw
    
    initparam$N <- N # store effective sample size
    #initparam$covmat <- (N-1)*covmat
    initparam$means <- means # store means
    
    mu0<-numeric(n)
    T0scale <- bgepar$am*(bgepar$aw-n-1)/(bgepar$am+1) # This follows from equations (19) and (20) of [GH2002]
    T0<-diag(T0scale,n,n)
    initparam$TN <- T0 + covmat + ((bgepar$am*N)/(bgepar$am+N))* (mu0 - means)%*%t(mu0 - means)
    initparam$awpN<-bgepar$aw+N
    constscorefact<- -(N/2)*log(pi) + (1/2)*log(bgepar$am/(bgepar$am+N))
    
    initparam$muN <- (N*means + bgepar$am*mu0)/(N + bgepar$am) # posterior mean mean
    initparam$SigmaN <- initparam$TN/(initparam$awpN-n-1) # posterior mode covariance matrix
    
    initparam$scoreconstvec<-numeric(n)
    for (j in (1:n)) {# j represents the number of parents plus 1
      awp<-bgepar$aw-n+j
      initparam$scoreconstvec[j]<-constscorefact - lgamma(awp/2) + lgamma((awp+N)/2) + ((awp+j-1)/2)*log(T0scale)
    }
  } else if (scoretype=="bde") {
    if(is.null(bdepar$chi)) {bdepar$chi<-0.5}
    if(is.null(bdepar$edgepf)) {bdepar$edgepf <- 2}
    if (is.null(weightvector)) {
      initparam$N<-nrow(data)
      initparam$d1<-data
      initparam$d0<-(1-data)
    } else {
      initparam$N<-sum(weightvector)
      initparam$d1<-data*weightvector
      initparam$d0<-(1-data)*weightvector
    }
    maxparents<-n-1
    initparam$scoreconstvec<-rep(0,maxparents+1)
    initparam$chi<-bdepar$chi #1
    initparam$pf<-bdepar$edgepf
    for(i in 0:maxparents){ # constant part of the score depending only on the hyperparameters
      noparams<-2^i
      initparam$scoreconstvec[i+1]<-noparams*lgamma(initparam$chi/noparams)-2*noparams*lgamma(initparam$chi/(2*noparams))-i*log(initparam$pf) #
    }
  } else if (scoretype=="bdecat") {
    if(is.null(bdecatpar$chi)) {bdecatpar$chi<-0.5}
    if(is.null(bdecatpar$edgepf)) {bdecatpar$edgepf <- 2}
    maxparents<-n-1
    initparam$chi<-bdecatpar$chi 
    initparam$pf<-bdecatpar$edgepf
    initparam$scoreconstvec <- -c(0:maxparents)*log(initparam$pf) # just edge penalisation here
    initparam$Cvec <- apply(initparam$data,2,max)+1 # number of levels of each variable
  } else if (scoretype=="usr"){
    if(is.null(usrpar$pctesttype)){usrpar$pctesttype <- "usr"}
    initparam$pctesttype <- usrpar$pctesttype
    initparam <- usrscoreparameters(initparam, usrpar)
  }  else if (scoretype=="dbn"){
    initparam$bgnodes<-c(1:nsmall+nsmall)
    nsmall <- n/dbnpar$slices # number of nodes in each slice
    
    # other slices we layer the data, 
    datalocal <- data[,1:(2*nsmall)]
    if(dbnpar$slices > 2){ # layer on later time slices
      for(jj in 1:(dbnpar$slices-2)){
        datalocal <- rbind(datalocal,data[,nsmall*jj+1:(2*nsmall)])
      }
    }
    
    # and have earlier times on the right hand side!
    datalocal <- datalocal[,c(1:nsmall+nsmall,1:nsmall)]
    
    initparam$otherslices <- scoreparameters(n=2*nsmall, scoretype=dbnpar$dbnscoretype, 
                                             datalocal, weightvector=weightvector, 
                                             bgnodes=c(1:nsmall+nsmall),
                                             bgepar=bgepar, bdepar=bdepar, bdecatpar=bdecatpar, dbnpar=dbnpar, 
                                             edgepmat=edgepmat, nodeslabels=NULL)
    
    # first slice we just take the first block
    bdecatpar$edgepf <- 1 # we don't want any additional edge penalisation
    bdepar$edgepf <- 1
    datalocal <- data[,bgnodes]
    initparam$firstslice <- scoreparameters(n=nsmall, scoretype=dbnpar$dbnscoretype, 
                                            datalocal, weightvector=weightvector, 
                                            bgepar=bgepar, bdepar=bdepar, bdecatpar=bdecatpar, dbnpar=dbnpar, 
                                            edgepmat=edgepmat, nodeslabels=NULL)
  }  
  attr(initparam, "class") <- "scoreparameters"
  initparam
}



#'Initialising score object
#'
#'This function returns an object of class scoreparameters containing the data and parameters needed for calculation of the BDe/BGe score, or a user defined score.
#' @param n number of nodes (variables) in the Bayesian network (excluding background nodes)
#' @param data the data matrix with n columns (the number of variables) and a number of rows equal to the number of observations
#' @param scoretype the score to be used to assess the DAG structure: "bge" for Gaussian data, "bde" for binary data, "bdecat" for categorical data, "usr" for a user defined score
#' @param weightvector (optional) a numerical vector of positive values representing the weight of each observation; should be NULL(default) for non-weighted data 
#' @param bgnodes (optional) a numerical vector which contains numbers of columns in the data defining background nodes, background nodes are nodes which have no parents but can be parents of other nodes in the network
#' @param bgepar a list which contains parameters for BGe score:
#' \itemize{
#' \item am (optional) a positive numerical value, 1 by default
#' \item aw (optional) a positive numerical value should be more than \code{n+1}, \code{n+am+1} by default
#' }
#' @param bdepar a list which contains parameters for BDe score for binary data:
#' \itemize{
#' \item chi (optional) a positive number of prior pseudo counts used by the BDe score, 0.5 by default
#' \item edgepf (optional) a positive numerical value providing the edge penalization factor to be combined with the BDe score, 2 by default
#' }
#' @param bdecatpar a list which contains parameters for BDe score for categorical data:
#' \itemize{
#' \item chi (optional) a positive number of prior pseudo counts used by the BDe score, 0.5 by default
#' \item edgepf (optional) a positive numerical value providing the edge penalization factor to be combined with the BDe score, 2 by default
#' }
#' @param usrpar a list which contains parameters for the user defined score
#' \itemize{
#' \item pctesttype (optional) "bde" to use the independence test for binary data, "bge" for continuous data
#' }
#'@param edgepmat (optional) a matrix of positive numerical values providing the per edge penalization factor to be added to the score, NULL by default
#'@param nodeslabels (optional) a vector of characters which denote the names of nodes in the Bayesian network
#'@return an object of class \code{scoreparameters}, which includes all necessary information for calculating the BDe/BGe score
#'@references Geiger D and Heckerman D (2002). Parameter priors for directed acyclic graphical models and the characterization of several probability distributions. The Annals of Statistics 30, 1412-1440.
#'@references Kuipers J, Moffa G and Heckerman D (2014). Addendum on the scoring of Gaussian acyclic graphical models. The Annals of Statistics 42, 1689-1691.
#'@references Heckerman D and Geiger D (1995). Learning Bayesian networks: A unification for discrete and Gaussian domains. In Eleventh Conference on Uncertainty in Artificial Intelligence, pages 274-284.
#'@references Scutari M (2016). An Empirical-Bayes Score for Discrete Bayesian Networks. Journal of Machine Learning Research 52, 438-448 
#'@examples
#' myDAG<-pcalg::randomDAG(20, prob=0.15, lB = 0.4, uB = 2) 
#' myData<-pcalg::rmvDAG(200, myDAG) 
#' myScore<-scoreparameters(20, "bge", myData)
# a constructor function for the "scoreparameters" class
scoreparameters.tested<-function(n, scoretype=c("bge","bde","bdecat","usr"), data, weightvector=NULL,
                          bgnodes=NULL,
                          bgepar=list(am=1, aw=NULL),
                          bdepar=list(chi=0.5, edgepf=2),
                          bdecatpar=list(chi=0.5, edgepf=2),
                          usrpar=list(pctesttype=c("bge","bde","bdecat")), 
                          edgepmat=NULL, nodeslabels=NULL) {
  
  if (!(scoretype%in%c("bge", "bde", "bdecat", "usr"))) {
    stop("Scoretype should be bge (for continuous data), bde (for binary data) bdecat (for categorical data) or usr (for user defined)")
  }
  
  if (anyNA(data)) {
    stop("Dataset contains missing data")  
  }
  
  if (ncol(data)!=n) {
    stop("n and number of columns in data do not match")
  }
  
  if (!is.null(weightvector)) {
    if (length(weightvector)!=nrow(data)) {
      stop("Length of the weightvector does not match the number of columns (observations) in data")
    }
  }
  
  if (scoretype=="bde") {
    if (!all(sapply(data,function(x)x%in%c(0,1)))) {
      stop("Dataset contains non-binary values")  
    }
  }
  
  if (scoretype=="bdecat") { # convert factors to levels  
    indx <- sapply(data, is.factor)
    data[indx] <- lapply(data[indx], function(x) as.numeric(x)-1)
    if (!all(unlist(lapply(data, function(x) setequal(unique(x),c(0:max(x))))))) {
      stop("Some variable levels are not present in the data")  
    }    
  }
  
  if (is.null(nodeslabels)) {
    if(all(is.character(colnames(data)))){
      nodeslabels<-colnames(data)
    } else {
      nodeslabels<-sapply(c(1:n), function(x)paste("v",x,sep=""))
    }
  }
  
  colnames(data)<-nodeslabels
  
  initparam<-list()
  initparam$labels<-nodeslabels
  initparam$type<-scoretype
  initparam$weightvector<-weightvector
  initparam$data<-data
  initparam$bgnodes<-bgnodes
  initparam$bgn<-length(bgnodes)
  initparam$nsmall<-n-initparam$bgn
  
  if (is.null(bdepar$edgepmat)) {
    initparam$logedgepmat <- NULL
  } else {
    initparam$logedgepmat <- log(bdepar$edgepmat)
  }
  
  if(scoretype=="bge") {
    if (is.null(weightvector)) {
      N<-nrow(data)
      covmat<-cov(data)*(N-1)
      means<-colMeans(data)
    } else {
      N<-sum(weightvector)
      forcov<-cov.wt(data,wt=weightvector,cor=TRUE,method="ML")
      covmat<-forcov$cov*N
      means<-forcov$center
    }
    if (is.null(bgepar$aw)) {
      bgepar$aw<-n+bgepar$am+1
    }
    
    initparam$am <- bgepar$am # store parameters
    initparam$aw <- bgepar$aw
    
    initparam$N <- N # store effective sample size
    #initparam$covmat <- (N-1)*covmat
    initparam$means <- means # store means
    
    mu0<-numeric(n)
    T0scale <- bgepar$am*(bgepar$aw-n-1)/(bgepar$am+1) # This follows from equations (19) and (20) of [GH2002]
    T0<-diag(T0scale,n,n)
    initparam$TN <- T0 + covmat + ((bgepar$am*N)/(bgepar$am+N))* (mu0 - means)%*%t(mu0 - means)
    initparam$awpN<-bgepar$aw+N
    constscorefact<- -(N/2)*log(pi) + (1/2)*log(bgepar$am/(bgepar$am+N))
    
    initparam$muN <- (N*means + bgepar$am*mu0)/(N + bgepar$am) # posterior mean mean
    initparam$SigmaN <- initparam$TN/(initparam$awpN-n-1) # posterior mode covariance matrix
    
    initparam$scoreconstvec<-numeric(n)
    for (j in (1:n)) {# j represents the number of parents plus 1
      awp<-bgepar$aw-n+j
      initparam$scoreconstvec[j]<-constscorefact - lgamma(awp/2) + lgamma((awp+N)/2) + ((awp+j-1)/2)*log(T0scale)
    }
  } else if (scoretype=="bde") {
    if(is.null(bdepar$chi)) {bdepar$chi<-0.5}
    if(is.null(bdepar$edgepf)) {bdepar$edgepf <- 2}
    if (is.null(weightvector)) {
      initparam$N<-nrow(data)
      initparam$d1<-data
      initparam$d0<-(1-data)
    } else {
      initparam$N<-sum(weightvector)
      initparam$d1<-data*weightvector
      initparam$d0<-(1-data)*weightvector
    }
    maxparents<-n-1
    initparam$scoreconstvec<-rep(0,maxparents+1)
    initparam$chi<-bdepar$chi #1
    initparam$pf<-bdepar$edgepf
    for(i in 0:maxparents){ # constant part of the score depending only on the hyperparameters
      noparams<-2^i
      initparam$scoreconstvec[i+1]<-noparams*lgamma(initparam$chi/noparams)-2*noparams*lgamma(initparam$chi/(2*noparams))-i*log(initparam$pf) #
    }
  } else if (scoretype=="bdecat") {
    if(is.null(bdecatpar$chi)) {bdecatpar$chi<-0.5}
    if(is.null(bdecatpar$edgepf)) {bdecatpar$edgepf <- 2}
    maxparents<-n-1
    initparam$chi<-bdecatpar$chi 
    initparam$pf<-bdecatpar$edgepf
    initparam$scoreconstvec <- -c(0:maxparents)*log(initparam$pf) # just edge penalisation here
    initparam$Cvec <- apply(initparam$data,2,max)+1 # number of levels of each variable
  } else if (scoretype=="usr"){
    if(is.null(usrpar$pctesttype)){usrpar$pctesttype <- "usr"}
    initparam$pctesttype <- usrpar$pctesttype
    initparam <- usrscoreparameters(initparam, usrpar)
  }  
  
  attr(initparam, "class") <- "scoreparameters"
  initparam
}




