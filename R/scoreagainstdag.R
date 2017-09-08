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
    stop("function implemented only for binary data and BDe score")
  } else {
    if (is.null(datatoscore)) {
      datatoscore<-scoreparam$data
    }
  }
  samplescores <- matrix(0,nrow=nrow(datatoscore),ncol=n)  
  for (j in 1:n)  {
    parentnodes <- which(incidence[,j]==1)
    samplescores[,j]<-scoreagainstDAGcore(j,parentnodes,n,scoreparam,datatoscore)
  }
  
  return(rowSums(samplescores))
}

scoreagainstDAGcore<-function(j,parentnodes,n,param,datatoscore) {
  samplenodescores<-rep(0,nrow(datatoscore)) # store
  lp<-length(parentnodes) # number of parents
  noparams<-2^lp # number of binary states of the parents
  switch(as.character(lp),
         "0"={# no parents
           N1<-sum(param$d1[,j],na.rm=TRUE)
           N0<-sum(param$d0[,j],na.rm=TRUE)
           NT<-N0+N1
           theta<-(N1+param$chi/(2*noparams))/(NT+param$chi/noparams) # the probability of each state
           samplenodescores[which(datatoscore[,j]==1)]<-log(theta) # log scores of 1s
           samplenodescores[which(datatoscore[,j]==0)]<-log(1-theta) # log scores of 0s
         },
         "1"={# one parent
           corescore<-param$scoreconstvec[lp+1]  
           summys<-param$data[,parentnodes]
           summysfull<-datatoscore[,parentnodes]
           
           for(i in 1:noparams-1){
             totest<-which(summys==i)
             N1<-sum(param$d1[totest,j],na.rm=TRUE)
             N0<-sum(param$d0[totest,j],na.rm=TRUE)
             NT<-N0+N1
             theta<-(N1+param$chi/(2*noparams))/(NT+param$chi/noparams) # the probability of each state
             toscore<-which(summysfull==i)
             samplenodescores[toscore[which(datatoscore[toscore,j]==1)]]<-log(theta) # log scores of 1s
             samplenodescores[toscore[which(datatoscore[toscore,j]==0)]]<-log(1-theta) # log scores of 0s
           }
         },     
         { # more parents
           corescore<-param$scoreconstvec[lp+1]  
           summys<-colSums(2^(c(0:(lp-1)))*t(param$data[,parentnodes]))
           tokeep<-which(!is.na(summys+param$d1[,j])) # remove NAs either in the parents or the child
           if(length(tokeep)<length(summys)){
             N1s<-collectC(summys[tokeep],param$d1[tokeep,j],noparams)
             N0s<-collectC(summys[tokeep],param$d0[tokeep,j],noparams)
           } else {
             N1s<-collectC(summys,param$d1[,j],noparams)
             N0s<-collectC(summys,param$d0[,j],noparams)
           }
           NTs<-N0s+N1s
           thetas<-(N1s+param$chi/(2*noparams))/(NTs+param$chi/noparams) # the probability of each state
           
           summysfull<-colSums(2^(c(0:(lp-1)))*t(datatoscore[,parentnodes]))
           ones<-which(datatoscore[,j]==1)
           samplenodescores[ones]<-log(thetas[summysfull[ones]+1])
           zeros<-which(datatoscore[,j]==0)
           samplenodescores[zeros]<-log(1-thetas[summysfull[zeros]+1])
         })
  
  return(samplenodescores)
}
