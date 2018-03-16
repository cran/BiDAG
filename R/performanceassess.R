#'Estimating posterior probabilities of single edges
#'
#'This function estimates the posterior probabilities of edges by averaging over a sample of DAGs
#'obtained via an MCMC scheme.
#'
#'@param MCMCchain list of square matrices with elements in \code{\{0,1\}} and representing adjacency matrices of a sample of DAGs obtained via an MCMC scheme
#'@param pdag logical, if TRUE (FALSE by default) all DAGs in the MCMCchain are first converted to equivalence class (CPDAG) before the averaging
#'@param burnin (optional) number between \code{0} and \code{1}, indicates the percentage of the samples which will be  the discarded as `burn-in' of the MCMC chain; the rest  of the samples will be used to calculate the posterior probabilities; 0.2 by default
#'@return a square matrix with dimensions equal to the number of variables; each entry \code{[i,j]} is an estimate of the posterior probability of the edge from node \code{i} to node \code{j}
#'@examples
#'Bostonscore<-scoreparameters(14, "bge", Boston)
#'\dontrun{
#'samplefit<-orderMCMC(14, Bostonscore, iterations=25000,chainout=TRUE)
#'MCMCchain<-samplefit$chain$incidence
#'edgesposterior<-edges.posterior(MCMCchain, burnin=0.2)
#'edgesposterior<-edges.posterior(MCMCchain, pdag=TRUE, burnin=0.2)
#'}
#'@export
edges.posterior<-function(MCMCchain,pdag=FALSE,burnin=0.2) {
  endstep<-length(MCMCchain)
  startstep<-as.integer(burnin*endstep)
  if (pdag) {
    cpdags<-lapply(MCMCchain[startstep:endstep],dagadj2cpadj)
    return(Reduce('+', cpdags)/(endstep-startstep+1))
  } else {
    return(Reduce('+', MCMCchain[startstep:endstep])/(endstep-startstep+1))
  }
}


#'Estimating a graph corresponding to a posterior probability threshold
#'
#'This function constructs a directed graph (not necessarily acyclic) including all edges with a posterior probability above a certain threshold.  The posterior probability is evaluated as the Monte Carlo estimate from a sample of DAGs obtained via an MCMC scheme.
#'
#'@param n number of nodes in the Bayesian network
#'@param MCMCchain list of adjacency matrices with dimensions equal to n and elements in \code{\{0,1\}}, representing a sample of DAGs from an MCMC scheme
#'@param pbarrier threshold such that only edges with a higher posterior probability will be retained in the directed graph summarising the sample of DAGs
#'@param pdag logical, if TRUE (FALSE by default) all DAGs in the MCMCchain are first converted to equivalence class (CPDAG) before the averaging
#'@param burnin (optional)  number between \code{0} and \code{1}, indicates the percentage of the samples which will be  the discarded as `burn-in' of the MCMC chain; the rest  of the samples will be used to calculate the posterior probabilities; 0.2 by default
#'@return a square matrix with dimensions equal to the number of variables representing the adjacency matrix of the directed graph summarising the sample of DAGs
#'@examples
#'Bostonscore<-scoreparameters(14, "bge", Boston)
#'\dontrun{
#'orderfit<-orderMCMC(14, Bostonscore, MAP=FALSE, iterations=25000, chainout=TRUE)
#'MCMCchain<-orderfit$chain$incidence
#'hdag<-dag.threshold(MCMCchain, pbarrier=0.9)
#'}
#'@export
dag.threshold<-function(n,MCMCchain, pbarrier, pdag=FALSE, burnin=0.2) {
  n<-nrow(MCMCchain[[1]])
  incidence<-matrix(rep(0, n*n), nrow=n, ncol=n)
  endstep<-length(MCMCchain)
  startstep<-as.integer(burnin*endstep)
  if (pdag) {
    cpdags<-lapply(MCMCchain[startstep:endstep],dagadj2cpadj)
    incidence[which(Reduce('+', cpdags)/(endstep-startstep+1)>pbarrier)]<-1
  } else {
  incidence[which(Reduce('+', MCMCchain[startstep:endstep])/(endstep-startstep+1)>pbarrier)]<-1
  }
  return(incidence)
}


#'Performance assessment of iterative MCMC scheme against a known Bayesian network
#'
#'This function calculates the number of true and false positives, the true positive rate, the structural Hamming distance and score for each iteration in the search procedure implemented in the \code{\link{iterativeMCMCsearch}} function
#'
#'@param MCMCmult an object which of class \code{MCMCmult}, contained in the output of the function \code{\link{iterativeMCMCsearch}}, when its \code{chainout} argument is set to TRUE; contains adjacency matrices sampled at each iteration of search space expansion; acessible by \code{MCMCout$chain}, where \code{MCMCout} is the output of function \code{\link{iterativeMCMCsearch}}  
#'@param truedag ground truth DAG which generated the data used in the search procedure; represented by an object of class \code{\link[graph]{graphNEL}}
#'@param sample logical (FALSE by default), indicates if \code{MCMCmult} contains sample or maximum score DAGs
#'@param cpdag logical, if TRUE (FALSE by default) all DAGs in the \code{MCMCmult} are first converted to their respective equivalence class (CPDAG) before the averaging if parameter \code{sample} set to TRUE
#'@param pbarrier threshold such that only edges with a higher posterior probability will be retained in the directed graph summarising the sample of DAGs at each iteration from \code{MCMCmult} if parameter \code{sample} set to TRUE
#'@return A matrix with the number of rows equal to the number of elements in \code{MCMCmult}, and 5 columns reporting for 
#'the maximally scoring DAG uncovered at each iteration (or for a summary over the sample of DAGs if \code{sample} parameter set to TRUE) 
#'the number of true positive edges (`TP'), the number of false positive edges (`FP'), 
#'the true positive rate (`TPR'), the structural Hamming distance (`SHD') and the score of the DAG (`SC'). 
#'Note that the maximum estimated DAG as well as the true DAG are first converted to 
#'the corresponding equivalence class (CPDAG) when calculating the SHD.
#' @examples
#' myDAG<-pcalg::randomDAG(20, prob=0.15, lB = 0.4, uB = 2) 
#' myData<-pcalg::rmvDAG(200, myDAG)
#' myScore<-scoreparameters(20, "bge", myData)
#'\dontrun{
#' MAPestimate<-iterativeMCMCsearch(20, myScore, chainout=TRUE, scoreout=TRUE)
#' iterations.check(MAPestimate, myDAG)
#' }
#'@export
iterations.check<-function(MCMCmult, truedag, sample=FALSE,cpdag=TRUE, pbarrier=0.5) {
  TP<-vector()
  FP<-vector()
  TPR<-vector()
  SHD<-vector()
  SC<-vector()
  trueskeleton<-dag2skeletonadjacency(truedag)
  numedges<-sum(trueskeleton)
  if (!sample) {
  for (j in  1:length(MCMCmult$chain)) {
    maxN<-which.max(unlist(MCMCmult$chain[[j]][[2]]))
    SC[j]<-MCMCmult$chain[[j]][[2]][[maxN]]
    maxadj<-MCMCmult$chain[[j]][[1]][[maxN]]
    estskelmcmc<-adjacency2skeleton(maxadj)
    diffmcmc<-estskelmcmc-trueskeleton
    mc.dag<-adjacency2dag(maxadj)
    cpmcmc<-pcalg::dag2cpdag(mc.dag)
    truecp<-pcalg::dag2cpdag(truedag)
    TP[j]<-numedges-sum(diffmcmc<0)
    FP[j]<-sum(diffmcmc>0)
    SHD[j]<-pcalg::shd(truecp, cpmcmc)
  }
  TPR<-TP/numedges
  result<-cbind(TP, FP, TPR, SHD, SC)
  colnames(result)<-c("TP", "FP", "TPR", "SHD", "SC")
  return(result) } else {
    n<-nrow(MCMCmult$max$DAG)
    for (j in  1:length(MCMCmult$chain)) {
      adj<-dag.threshold(n,MCMCmult$chain[[j]][[1]],pbarrier,pdag=cpdag)
      estskelmcmc<-adjacency2skeleton(adj)
      diffmcmc<-estskelmcmc-trueskeleton
      mc.dag<-dag2cpdag(adjacency2dag(adj))
      truecp<-pcalg::dag2cpdag(truedag)
      TP[j]<-numedges-sum(diffmcmc<0)
      FP[j]<-sum(diffmcmc>0)
      SHD[j]<-pcalg::shd(truecp, mc.dag)
    }
      TPR<-TP/numedges
      result<-cbind(TP, FP, TPR, SHD)
      colnames(result)<-c("TP", "FP", "TPR", "SHD")
      return(result)
  }
}

#'Performance assessment of sampling algorithms against a known Bayesian network
#'
#'This function calculates the number of true and false positives and the structural Hamming distance between a ground truth DAG and a directed graph summarising a sample of DAGs obtained from an MCMC scheme, as the posterior probability threshold is varied
#'
#'@param n number of nodes in the Bayesian network
#'@param MCMCchain list of adjacency matrices with dimensions equal to n and elements in \code{\{0,1\}}, representing a sample of DAGs from an MCMC scheme
#'@param truedag ground truth DAG which generated the data used in the search procedure; represented by an object of class  \code{\link[graph]{graphNEL}}
#'@param pbarrier (optional) a vector of numeric values between 0 and 1, defining posterior probabilities according to which the edges of assessed structures are drawn, please note very low barriers can lead to very dense structures; by default 
#'\eqn{pbarrier=c(0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2)}
#'@param pdag logical, if TRUE (default) all DAGs in the MCMCchain are first converted to equivalence class (CPDAG) before the averaging
#'@param burnin (optional)  number between \code{0} and \code{1}, indicates the percentage of the samples which will be  the discarded as `burn-in' of the MCMC chain; the rest  of the samples will be used to calculate the posterior probabilities; 0.2 by default
#'@return A matrix with the number of rows equal to the number of posterior thresholds tested, and 4 columns reporting for each  thresholded directed graphs the number of true positive edges (`TP'), the number of false positive edges (`FP'), the structural Hamming distance (`SHD') and the posterior threshold
#' @examples
#' myDAG<-pcalg::randomDAG(n=20, prob=0.1, lB = 0.4, uB = 2)
#' myData<-pcalg::rmvDAG(200,myDAG)
#' myScore<-scoreparameters(20, "bge", myData)
#' \dontrun{
#' ordersample<-orderMCMC(n=20, myScore, chainout=TRUE)
#' MCMCchain<-ordersample$chain$incidence
#' sample.check(MCMCchain, myDAG, pdag=TRUE, burnin=0.2)
#' }
#'@export
sample.check<-function(n,MCMCchain, truedag, pbarrier=c(0.99,0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2),
                                         pdag=TRUE, burnin=0.2) {
  n<-nrow(MCMCchain[[1]])
  trueskeleton<-dag2skeletonadjacency(truedag)
  results<-matrix(ncol=4,nrow=length(pbarrier))
  results[,4]<-pbarrier
  numedges<-sum(trueskeleton)
  endstep<-length(MCMCchain)
  startstep<-as.integer(burnin*endstep)
  if(pdag) {
  dags<-lapply(MCMCchain[startstep:endstep],dagadj2cpadj) #first convert every DAG in the sample to equivalence class
  } else {dags<-MCMCchain[startstep:endstep]}
  for (p in 1:length(pbarrier)) {
    sampledag<-matrix(0, nrow=n,ncol=n)
    sampledag[which(Reduce('+', dags)/(endstep-startstep+1)>pbarrier[p])]<-1 #average over CPDAGs
    sampleest<-adjacency2skeleton(sampledag)
    diffmcmc<-sampleest-trueskeleton
    mc.dag<-adjacency2dag(sampledag)
    truecp<-pcalg::dag2cpdag(truedag)
    results[p,1]<-numedges-sum(diffmcmc<0)
    results[p,2]<-sum(diffmcmc>0)
    results[p,3]<-pcalg::shd(truecp, mc.dag)
  }
  colnames(results)<-c("TP", "FP", "SHD", "post.thr.")
  return(results)
}
