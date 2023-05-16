#'Bayesian network structure sampling from the posterior distribution
#'
#'This function can be used for structure sampling using three different MCMC schemes. Order MCMC scheme (\code{algorithm="order"}) is the most computationally
#'efficient however it imposes a non-uniform prior in the space of DAGs. Partition MCMC (\code{algorithm="partition"}) is less computationally efficient and requires more iterations
#'to reach convergence, however it implements sampling using a uniform prior in the space of DAGs.
#'Due to the superexponential size of the search space as the number of nodes increases, the 
#'MCMC search is performed on a reduced search space. By default the search space is limited to the skeleton found through the PC algorithm by means of conditional independence tests 
#'(using the functions \code{\link[pcalg]{skeleton}} and \code{\link[pcalg]{pc}} from the `pcalg' package [Kalisch et al, 2012]).
#'It is also possible to define an arbitrary search space by inputting an adjacency matrix, for example estimated by partial correlations or other network algorithms.
#'Also implemented is the possibility to expand the default or input search space, by allowing each node in the network to have one additional parent.  
#'This offers improvements in the learning and sampling of Bayesian networks. The iterative MCMC scheme (\code{algorithm="orderIter"}) allows for iterative expansions of the search space.
#'This is useful in cases when the initial search space is poor in a sense that it contains only a limited number of true positive edges. Iterative expansions of the search space
#'efficiently solve this issue. However this scheme requires longer runtimes due to the need of running multiple consecutive MCMC chains.  
#'This function is a wrapper for the three individual structure learning and sampling functions that implement each of the described algorithms; for details see \code{\link{orderMCMC}},
#'\code{\link{partitionMCMC}},\code{\link{iterativeMCMC}}.
#' @param scorepar an object of class \code{scoreparameters}, containing the data and score parameters, see constructor function \code{\link{scoreparameters}}
#' @param algorithm MCMC scheme to be used for sampling from posterior distribution; possible options are "order" (\code{\link{orderMCMC}}), "orderIter" (\code{\link{iterativeMCMC}}) or "partition" (\code{\link{partitionMCMC}})
#' @param chainout logical, if TRUE the saved MCMC steps are returned, TRUE by default
#' @param scoreout logical, if TRUE the search space and score tables are returned, FALSE by default
#' @param moveprobs a numerical vector of 4 (for "order" and "orderIter" algorithms) or 5 values (for "partition" algorithm) representing probabilities of the different moves in the space of
#' order and partitions accordingly. The moves are described in the corresponding algorithm specific functions \code{\link{orderMCMC}} and \code{\link{partitionMCMC}}
#' @param iterations integer, the number of MCMC steps, the default value is \eqn{6n^{2}\log{n}} orderMCMC, \eqn{20n^{2}\log{n}} for partitionMCMC and \eqn{3.5n^{2}\log{n}} for iterativeMCMC; where n is the number of nodes in the Bayesian network
#' @param stepsave integer, thinning interval for the MCMC chain, indicating the number of steps between two output iterations, the default is \code{iterations/1000}
#' @param alpha numerical significance value in \code{\{0,1\}} for the conditional independence tests at the PC algorithm stage
#' @param gamma tuning parameter which transforms the score by raising it to this power, 1 by default
#' @param cpdag logical, if TRUE the CPDAG returned by the PC algorithm will be used as the search
#'space, if FALSE (default) the full undirected skeleton will be used as the search space
#' @param hardlimit integer, limit on the size of parent sets in the search space; 
#' @param verbose logical, if TRUE messages about the algorithm's progress will be printed, FALSE by default
#' @param compress logical, if TRUE adjacency matrices representing sampled graphs will be stored as a sparse Matrix (recommended); TRUE by default
#' @param startspace (optional) a square sparse or ordinary matrix, of dimensions equal to the number of nodes, which defines the search space for the order MCMC in the form of an adjacency matrix. If NULL, the skeleton obtained from the PC-algorithm will be used. If \code{startspace[i,j]} equals to 1 (0) it means that the edge from node \code{i} to node \code{j} is included (excluded) from the search space. To include an edge in both directions, both \code{startspace[i,j]} and \code{startspace[j,i]} should be 1.
#' @param blacklist (optional) a square sparse or ordinary matrix, of dimensions equal to the number of nodes, which defines edges to exclude from the search space. If \code{blacklist[i,j]} equals to 1 it means that the edge from node \code{i} to node \code{j} is excluded from the search space.
#' @param scoretable (optional) object of class \code{scorespace} containing list of score tables calculated for example by the last iteration of the function \code{iterativeMCMC}. When not NULL, parameter \code{startspace} is ignored.
#' @param startpoint (optional) integer vector of length n (representing an order when \code{algorithm="order"} or \code{algorithm="orderIter"}) or an adjacency matrix or sparse adjacency matrix (representing a DAG when \code{algorithm="partition"}), which will be used as the starting point in the MCMC algorithm, the default starting point is random
#' @param plus1 logical, if TRUE (default) the search is performed on the extended search space; only changable for orderMCMC; for other algorithms is fixed to TRUE
#' @param iterpar addition list of parameters for the MCMC scheme implemeting iterative expansions of the search space; for more details see \code{\link{iterativeMCMC}}; list(posterior = 0.5, softlimit = 9, mergetype = "skeleton", accum = FALSE, 
#'plus1it = NULL, addspace = NULL, alphainit = NULL)
#' @return Depending on the value or the parameter \code{algorithm} returns an object of class \code{orderMCMC}, \code{partitionMCMC} or \code{iterativeMCMC} which contains log-score trace of sampled DAGs as well 
#' as adjacency matrix of the maximum scoring DAG(s), its score and the order or partition score. The output can optionally include DAGs sampled in MCMC iterations and the score tables. 
#' Optional output is regulated by the parameters \code{chainout} and \code{scoreout}. See \code{\link{orderMCMC class}}, \code{\link{partitionMCMC class}}, \code{\link{iterativeMCMC class}} for a detailed description of the classes' structures.
#' @note see also extractor functions \code{\link{getDAG}}, \code{\link{getTrace}}, \code{\link{getSpace}}, \code{\link{getMCMCscore}}.
#'@references P. Suter, J. Kuipers, G. Moffa, N.Beerenwinkel (2023) <doi:10.18637/jss.v105.i09>
#'@references Friedman N and Koller D (2003). A Bayesian approach to structure discovery in bayesian networks. Machine Learning 50, 95-125.
#'@references Kalisch M, Maechler M, Colombo D, Maathuis M and Buehlmann P (2012). Causal inference using graphical models with the R package pcalg. Journal of Statistical Software 47, 1-26.
#'@references Geiger D and Heckerman D (2002). Parameter priors for directed acyclic graphical models and the characterization of several probability distributions. The Annals of Statistics 30, 1412-1440.
#'@references Kuipers J, Moffa G and Heckerman D (2014). Addendum on the scoring of Gaussian acyclic graphical models. The Annals of Statistics 42, 1689-1691.
#'@references Spirtes P, Glymour C and Scheines R (2000). Causation, Prediction, and Search, 2nd edition. The MIT Press.
#'@examples
#'\dontrun{
#'Asiascore <- scoreparameters("bde", Asia)
#'iterativefit <- learnBN(Asiascore, algorithm = "orderIter")
#'orderfit <- sampleBN(Asiascore, scoretable = iterativefit)
#'
#'myScore<-scoreparameters("bge",Boston)
#'MCMCchains<-list()
#'MCMCchains[[1]]<-sampleBN(myScore,"partition")
#'MCMCchains[[2]]<-sampleBN(myScore,"partition")
#'edge_posterior<-lapply(MCMCchains,edgep,pdag=TRUE)
#'plotpcor(edge_posterior)
#'}
#'@author Polina Suter, Jack Kuipers, the code partly derived from the order MCMC implementation from Kuipers J, Moffa G (2017) <doi:10.1080/01621459.2015.1133426>
#'@export
sampleBN<-function(scorepar, algorithm = c("order", "orderIter", "partition"), chainout = TRUE,
                   scoreout = FALSE, alpha = 0.05, moveprobs = NULL, iterations = NULL, stepsave = NULL, 
                   gamma = 1, verbose = FALSE, compress = TRUE, startspace = NULL, blacklist = NULL, 
                   scoretable = NULL, startpoint = NULL, plus1 = TRUE, cpdag = FALSE, hardlimit = 12,
                   iterpar = list(posterior = 0.5, softlimit = 9, mergetype = "skeleton", accum = FALSE, 
                                  plus1it = NULL, addspace = NULL, alphainit = NULL)){
  
  if(length(algorithm)>1) algorithm<-"order"
  
  if(algorithm=="order") {
    
    MCMCresult<-orderMCMC(scorepar,MAP=FALSE,chainout=chainout,scoreout=scoreout,alpha=alpha,moveprobs=moveprobs,
                          iterations=iterations,stepsave=stepsave,gamma=gamma,verbose=verbose,compress=compress,
                          startspace=startspace,blacklist=blacklist,scoretable=scoretable,startorder=startpoint,
                          plus1=plus1,cpdag=cpdag,hardlimit=hardlimit)
    
  } else if (algorithm=="partition") {
    if(is(scoretable,"iterativeMCMC")){
      scoretable<-getSpace(scoretable)
    }
    
    MCMCresult<-partitionMCMC(scorepar,scoreout=scoreout,alpha=alpha,moveprobs=moveprobs,
                              iterations=iterations,stepsave=stepsave,gamma=gamma,verbose=verbose,compress=compress,
                              startspace=startspace,blacklist=blacklist,scoretable=scoretable,startDAG=startpoint)
    
  } else {
    
    iterpardef<-list(posterior = 0.5, softlimit = 9, mergetype = "skeleton", accum = FALSE, plus1it = NULL, addspace = NULL, alphainit = NULL)
    iterpardef[names(iterpar)]<-iterpar[names(iterpar)]
    iterpar<-iterpardef
    
    MCMCresult<-iterativeMCMC(scorepar,MAP=FALSE,chainout=chainout,scoreout=scoreout,alpha=alpha,moveprobs=moveprobs,
                              iterations=iterations,stepsave=stepsave,gamma=gamma,verbose=verbose,compress=compress,
                              startspace=startspace,blacklist=blacklist,scoretable=scoretable,startorder=startpoint,
                              cpdag=cpdag,hardlimit=hardlimit,posterior=iterpar$posterior,softlimit=iterpar$softlimit,
                              mergetype=iterpar$mergetype,accum=iterpar$accum,plus1it=iterpar$plus1it,
                              addspace=iterpar$addspace,alphainit=iterpar$alphainit)
  }
  return(MCMCresult)
}

