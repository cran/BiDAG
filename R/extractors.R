#' Extracting adjacency matrix (DAG) from MCMC object
#' 
#' This function extracts an adjacency matrix of
#' a maximum scoring DAG from the result of the MCMC run.
#' 
#' @param x object of class 'orderMCMC','partitionMCMC' or 'iterativeMCMC'
#' @param amat logical, when TRUE adjacency matrix is returned and object of class 'graphNEL' otherwise
#' @param cp logical, when TRUE the CPDAG (equivalence class) is returned and DAG otherwise; FALSE by default
#' @return adjacency matrix of a maximum scoring DAG (or CPDAG) discovered/sampled in one MCMC run
#' @examples
#'myscore<-scoreparameters("bge", Boston)
#'\dontrun{
#'itfit<-learnBN(myscore,algorithm="orderIter")
#'maxEC<-getDAG(itfit,cp=TRUE)
#'}
#' @export
getDAG<-function(x,amat=TRUE,cp=FALSE) {
  if (cp) graphy<-x$CPDAG else graphy<-x$DAG
  if(amat)  return(graphy) else return(m2graph(graphy))
}

#' Extracting scorespace from MCMC object
#' 
#' This function extracts an object of class 'scorespace'
#' from the result of the MCMC run when the parameter 'scoreout' was set to TRUE; otherwise extracts
#' only adjacency matrix of the final search space without the score tables.
#' 
#' @param x object of class 'orderMCMC','partitionMCMC' or 'iterativeMCMC'
#' @return an object of class 'scorespace' or an adjacency binary matrix corresponding to a search space last used in MCMC
#' @examples
#'myscore<-scoreparameters("bge", Boston)
#'\dontrun{
#'itfit<-learnBN(myscore,algorithm="orderIter",scoreout=TRUE)
#'itspace<-getSpace(itfit)
#'}
#' @export
getSpace<-function(x) {
  if(is.null(x$scoretable)) {
    warning("object x does not contain score tables! set the parameter 'scoreout' to TRUE when running MCMC 
            only adjacency matrix is returned")
    return(x$endspace)
  } else {
    return(x$scoretable)
  }
}

#' Extracting score from MCMC object
#' 
#' This function extracts the score of a maximum DAG sampled in the MCMC run.
#' 
#' @param x object of class 'orderMCMC','partitionMCMC' or 'iterativeMCMC'
#' @return a score of a maximum-scoring DAG found/sampled in one MCMC run
#' @examples
#'myscore<-scoreparameters("bge", Boston)
#'\dontrun{
#'itfit<-learnBN(myscore,algorithm="orderIter")
#'getMCMCscore(itfit)
#'}
#' @export
getMCMCscore<-function(x) {
    return(x$score)
}

#' Extracting trace from MCMC object
#' 
#' This function extracts a trace of
#' \itemize{
#' \item DAG scores 
#' \item DAG adjacency matrices
#' \item orders
#' \item order scores 
#' }
#' from the result of the MCMC run. Note that the last three options
#' work only when the parameter 'scoreout' was set to TRUE.
#' 
#' @param x object of class 'orderMCMC','partitionMCMC' or 'iterativeMCMC'
#' @param which integer, indication which trace is returned: DAG scores (which = 0), DAGs (which = 1),
#' orders (which = 2), order scores (which = 3) 
#' @return a list or a vector of objects representing MCMC trace, depends on parameter 'which'; by default, the trace of DAG scores is returned
#' @examples
#'myscore<-scoreparameters("bge",Boston)
#'\dontrun{
#'orderfit<-sampleBN(myscore,algorithm="order")
#'DAGscores<-getTrace(orderfit,which=0)
#'DAGtrace<-getTrace(orderfit,which=1)
#'orderscores<-getTrace(orderfit,which=3)
#'}
#' @export
getTrace<-function(x,which=0) {
  if(is.null(x$traceadd) & which>0) {
    warning("the result does not contain all sampled DAGs (only maximum)! 
            set the parameter 'chainout' to TRUE when running MCMC 
            returning DAG scores")
    return(x$trace)
  } else {
    if(which==0) return(x$trace) else return(x$traceadd[[which]])
  }
}


#' Extracting runtime
#' 
#' This function extracts runtime of a particular step of order and partition MCMC.
#' 
#' @param x object of class 'orderMCMC'or 'partitionMCMC'
#' @param which integer, defines if the runtime is extracted for: computing score tables (which = 1), running MCMC chain (which = 2)
#' @return runtime of a particular step of MCMC scheme or total runtime
#' @examples
#'myscore<-scoreparameters("bge",Boston)
#'\dontrun{
#'orderfit<-sampleBN(myscore,algorithm="order")
#'(getRuntime(orderfit,1))
#'(getRuntime(orderfit,2))
#'}
#' @export
getRuntime<-function(x,which=0) {
    if(which==1) return(x$info$runtimes["scoretables"]) else return(x$info$runtimes["MCMCchain"])
}

