#'Calculating the score of a sample against a DAG
#'
#'This function calculates the score of a given sample against a DAG represented by its incidence matrix. 
#'@param n number of nodes in the Bayesian network
#'@param scoreparam an object of class \code{scoreparameters}; see constructor function \code{\link{scoreparameters}}
#'@param incidence a square matrix of dimensions equal to the number of variables with entries in \code{\{0,1\}}, representing the adjacency matrix of the DAG against which the score is calculated
#'@param datatoscore (optional) a matrix (vector) containing binary observations (or just one observation)  to be scored; the number of columns should be equal to the number of variables in the Bayesian network, the number of rows should be equal to the number of observations; by default all data from \code{scoreparam} parameter is used
#'@return the log of the BDe score of given observations against a DAG
#'@references Heckerman D and Geiger D, (1995). Learning Bayesian networks: A unification for discrete and Gaussian domains. In Eleventh Conference on Uncertainty in Artificial Intelligence, pages 274-284, 1995.
#'@examples
#'\dontrun{
#' myDAG<-pcalg::randomDAG(20, prob=0.15, lB = 0.4, uB = 2) 
#' myData<-pcalg::rmvDAG(200, myDAG) 
#' adjacency<-dag2adjacencymatrix(myDAG)
#' param<-scoreparameters()
#' scoreagainstDAG(20, myData, adjacency, "bge")
#' }
#'@export
scoreagainstDAG <- function(n, scoreparam, incidence, datatoscore=NULL){
  if (scoreparam$type=="bge") {
    stop("function implemented only for binary data and BDe score")}
  if (!is.null(datatoscore)) {
    if (scoreparam$type=="bge") {
      stop("function implemented only for binary data and BDe score")
      scoreparam$data<-datatoscore
    } else {
      if (is.null(scoreparam$weightvector)) {
        scoreparam$d1<-t(datatoscore)
        scoreparam$d0<-t((1-datatoscore))
      } else {
        scoreparam$d1<-t(datatoscore*scoreparam$weightvector)
        scoreparam$d0<-t((1-datatoscore)*scoreparam$weightvector)
      }
      scoreparam$data<-t(datatoscore)
    }
  }
  samplescores <- matrix(0,nrow=n,ncol=ncol(scoreparam$data))  
  for (j in 1:n)  {
    parentnodes <- which(incidence[,j]==1)
    samplescores[j,]<-scoreagainstDAGcore(j,parentnodes,n,scoreparam)
  }
  return(colSums(samplescores))
}


scoreagainstDAGcore<-function(j,parentnodes,n,param) {
  samplenodescores<-rep(0,ncol(param$data)) # store
  lp<-length(parentnodes) # number of parents
  noparams<-2^lp # number of binary states of the parents
  switch(as.character(lp),
         "0"={# no parents
           N1<-sum(param$d1[j,])
           N0<-sum(param$d0[j,])
           NT<-N0+N1
           theta<-(N1+param$chi/(2*noparams))/(NT+param$chi/noparams) # the probability of each state
           samplenodescores[which(param$data[j,]==1)]<-log(theta) # log scores of 1s
           samplenodescores[which(param$data[j,]==0)]<-log(1-theta) # log scores of 0s
         },
         "1"={# one parent
           corescore<-param$scoreconstvec[lp+1]  
           summys<-param$data[parentnodes,]
           summysfull<-param$data[parentnodes,]
           
           for(i in 1:noparams-1){
             totest<-which(summys==i)
             N1<-sum(param$d1[j,totest])
             N0<-sum(param$d0[j,totest])
             NT<-N0+N1
             theta<-(N1+param$chi/(2*noparams))/(NT+param$chi/noparams) # the probability of each state
             toscore<-which(summysfull==i)
             samplenodescores[toscore[which(param$data[j,toscore]==1)]]<-log(theta) # log scores of 1s
             samplenodescores[toscore[which(param$data[j,toscore]==0)]]<-log(1-theta) # log scores of 0s
           }
         },     
         { # more parents
           corescore<-param$scoreconstvec[lp+1]  
           summys<-colSums(2^(c(0:(lp-1)))*param$data[parentnodes,])
           summysfull<-colSums(2^(c(0:(lp-1)))*param$data[parentnodes,])
           
           for(i in 1:noparams-1){
             totest<-which(summys==i)
             N1<-sum(param$d1[j,totest])
             N0<-sum(param$d0[j,totest])
             NT<-N0+N1
             theta<-(N1+param$chi/(2*noparams))/(NT+param$chi/noparams) # the probability of each state
             toscore<-which(summysfull==i)
             samplenodescores[toscore[which(param$data[j,toscore]==1)]]<-log(theta) # log scores of 1s
             samplenodescores[toscore[which(param$data[j,toscore]==0)]]<-log(1-theta) # log scores of 0s
           }
         })
  
  return(samplenodescores)
}
