#'Structure learning with the order MCMC algorithm
#'
#'This function implements the order MCMC algorithm for the structure learning of Bayesian networks. This function can be used
#'for MAP discovery and for sampling from the posterior distribution of DAGs given the data.
#'Due to the superexponential size of the search space as the number of nodes increases, the 
#'MCMC search is performed on a reduced search space.
#'By default the search space is limited to the skeleton found through the PC algorithm by means of conditional independence tests 
#'(using the functions \code{\link[pcalg]{skeleton}} and \code{\link[pcalg]{pc}} from the `pcalg' package [Kalisch et al, 2012]).
#'It is also possible to define an arbitrary search space by inputting an adjacency matrix, for example estimated by partial correlations or other network algorithms.
#'Also implemented is the possibility to expand the default or input search space, by allowing each node in the network to have one additional parent.  This offers improvements in the learning and sampling of Bayesian networks. 
#' @param n number of nodes in the Bayesian network
#' @param scoreparam an object of class \code{scoreparameters}, containing the data and score parameters, see constructor function \code{\link{scoreparameters}}
#' @param MAP logical, if TRUE (default) the search targets the MAP DAG (a DAG with maximum score),
#' if FALSE at each MCMC step a DAG is sampled from the order proportionally to its score
#' @param plus1 logical, if TRUE (default) the search is performed on the extended search space
#' @param moveprobs (optional) a numerical vector of 3 values in \code{\{0,1\}} corresponding to the probabilities of the following MCMC moves in the order space
#' \itemize{
#' \item exchanging 2 random nodes in the order
#' \item exchanging 2 adjacent nodes in the order
#' \item placing a single node elsewhere in the order
#' \item staying still
#' }
#' @param iterations (optional) integer, the number of MCMC steps, the default value is \eqn{5n^{2}\log{n}}
#' @param stepsave (optional) integer, thinning interval for the MCMC chain, indicating the number of steps between two output iterations, the default is \code{iterations/1000}
#' @param alpha (optional) numerical significance value in \code{\{0,1\}} for the conditional independence tests at the PC algorithm stage (by default \eqn{0.4} for \eqn{n<50}, \eqn{20/n} for \eqn{n>50})
#' @param gamma (optional) tuning parameter which transforms the score by raising it to this power, 1 by default
#' @param startspace (optional) a square matrix, of dimensions equal to the number of nodes, which defines the search space for the order MCMC in the form of an adjacency matrix. If NULL, the skeleton obtained from the PC-algorithm will be used. If \code{startspace[i,j]} equals to 1 (0) it means that the edge from node \code{i} to node \code{j} is included (excluded) from the search space. To include an edge in both directions, both \code{startspace[i,j]} and \code{startspace[j,i]} should be 1.
#' @param blacklist (optional) a square matrix, of dimensions equal to the number of nodes, which defines edges to exclude from the search space. If \code{blacklist[i,j]} equals to 1 it means that the edge from node \code{i} to node \code{j} is excluded from the search space.
#' @param scoretable (optional) list of score tables calculated for example by the last iteration of the iterativeMCMCsearch function, to avoid their recomputation  The score tables must match the permissible parents in the search space defined by the startspace parameter.
#' @param startorder (optional) integer vector of length n, which will be used as the starting order in the MCMC algorithm, the default order is \code{c(1:n)}
#' @param cpdag (optional) logical, if TRUE the CPDAG returned by the PC algorithm will be used as the search
#'space, if FALSE (default) the full undirected skeleton will be used as the search space
#' @param chainout logical, if TRUE the saved MCMC steps are returned, FALSE by default
#' @param scoreout logical, if TRUE the search space and score tables are returned, FALSE by default
#' @param verbose logical, if TRUE messages about the algorithm's progress will be printed, FALSE by default
#' @return Depends on the logical parameters \code{chainout} and \code{scoreout}.
#'  If both are FALSE (default), an object of class \code{MCMCmax}, containing a list of 3 elements:
#' \itemize{
#' \item {DAG -} {the adjacency matrix of the DAG with maximal score}
#' \item {order -} {an order it belongs to}
#' \item {score -} {the score of the reported DAG}
#' }
#' If \code{chainout} is TRUE an object of class \code{MCMCtrace} is additionally returned, contains 4 lists  (each of the 4 lists has length
#' \code{iterations/stepsave}, i.e. the number of saved MCMC steps):
#' \itemize{
#' \item incidence - contains a list of adjacency matrices of DAGs sampled at each step of MCMC
#' \item DAGscores - contains a list of scores of DAGs sampled at each step of MCMC
#' \item orderscores - contains a list of scores of orders of DAGs sampled at each step of MCMC
#' \item order - contains a list of permutations of the nodes of DAGs sampled at each step of MCMC
#' }
#' If \code{scoreout} is TRUE an object of class \code{MCMCspace} is additionally returned, contains a list of 2 elements:
#' \itemize{
#'  \item adjacency - the adjacency matrix representing the search space
#'  \item scoretable - the list of score tables corresponding to this search space
#' }
#'@references Friedman N and Koller D (2003). A Bayesian approach to structure discovery in bayesian networks. Machine Learning 50, 95-125.
#'@references Kalisch M, Maechler M, Colombo D, Maathuis M and Buehlmann P (2012). Causal inference using graphical models with the R package pcalg. Journal of Statistical Software 47, 1-26.
#'@references Geiger D and Heckerman D (2002). Parameter priors for directed acyclic graphical models and the characterization of several probability distributions. The Annals of Statistics 30, 1412-1440.
#'@references Kuipers J, Moffa G and Heckerman D (2014). Addendum on the scoring of Gaussian acyclic graphical models. The Annals of Statistics 42, 1689-1691.
#'@references Spirtes P, Glymour C and Scheines R (2000). Causation, Prediction, and Search, 2nd edition. The MIT Press.
#'@examples
#'\dontrun{
#'#find a MAP DAG with search space defined by PC and plus1 neighbourhood
#'Bostonscore<-scoreparameters(14,"bge",Boston)
#'orderMAPfit<-orderMCMC(14,Bostonscore)
#'orderMAPfit$max$score
#'#sample DAGs with order MCMC
#'ordersamplefit<-orderMCMC(14,Bostonscore,MAP=FALSE,chainout=TRUE)
#'}
#'
#'@export

orderMCMC<-function(n, scoreparam, MAP=TRUE, plus1=TRUE,
                    startspace=NULL, blacklist=NULL,startorder=c(1:n), scoretable=NULL,  
                    moveprobs=NULL, iterations=NULL, stepsave=NULL, alpha=NULL, cpdag=FALSE, gamma=1,
                    chainout=FALSE, scoreout=FALSE, verbose=FALSE) {
  if (is.null(moveprobs)) { 
    prob1<-99
    if(n>3){ prob1<-round(6*99*n/(n^2+10*n-24)) }
    prob1<-prob1/100
    moveprobs<-c(prob1,0.99-prob1,0.01)
    moveprobs<-moveprobs/sum(moveprobs)
    moveprobs<-c(moveprobs[c(1,2)],0,moveprobs[3])
  }
  if(is.null(iterations)){
    if(n<26){
      iterations<-30000
    } else {
      iterations<-(5*n*n*log(n))-(5*n*n*log(n)) %% 1000
    }
  }
  if(is.null(stepsave)){
    stepsave<-floor(iterations/1000)
  }
  
  result<-iterativeMCMCplus1(n,param=scoreparam,iterations,stepsave,plus1it=1,startorder=startorder,
                             moveprobs=moveprobs,alpha=alpha,cpdag=cpdag,scoretable=scoretable,
                             plus1=plus1,MAP=MAP,chainout=chainout, scoreout=scoreout,
                             startspace=startspace,blacklist=blacklist,gamma=gamma,verbose=verbose)
  return(result)
  
}

#'DAG structure sampling with partition MCMC
#'
#'This function implements the partition MCMC algorithm for the structure learning of Bayesian networks.  This procedure provides an unbiased sample from the posterior distribution of DAGs given the data. The search space can be defined either by a preliminary run of the iterativeMCMCsearch function or by a given adjacency matrix (which can be the full matrix with zero on the diagonal, to consider the entire space of DAGs, feasible only for a limited number of nodes). 
#'
#' @param n number of nodes in the Bayesian network
#' @param scoreparam an object of class \code{scoreparameters}, containing the data and scoring parameters;  see constructor function \code{\link{scoreparameters}}.
#' @param startspace (optional) a square matrix, of dimensions equal to the number of nodes, which defines the search space for the order MCMC in the form of an adjacency matrix; if NULL, the skeleton obtained from the PC-algorithm will be used. If \code{startspace[i,j]} equals to 1 (0) it means that the edge from node \code{i} to node \code{j} is included (excluded) from the search space. To include an edge in both directions, both \code{startspace[i,j]} and \code{startspace[j,i]} should be 1.
#' @param blacklist (optional) a square matrix, of dimensions equal to the number of nodes, which defines edges to exclude from the search space; if \code{blacklist[i,j]=1} it means that the edge from node \code{i} to node \code{j} is excluded from the search space
#' @param scoretable (optional) list of score tables calculated for example by the last iteration of the iterativeMCMCsearch function, to avoid their recomputation; the score tables must match the permissible parents in the search space defined by the startspace parameter
#' @param startDAG (optional) an adjacency matrix of dimensions equal to the number of nodes, representing a DAG in the search space defined by startspace.  If startspace is defined but \code{startDAG} is not, an empty DAG will be used by default
#' @param moveprobs (optional) a numerical vector of 5 values in \code{\{0,1\}} corresponding to the following MCMC move probabilities in the space of partitions:
#' \itemize{
#' \item swap any two elements from different partition elements
#' \item swap any two elements in adjacent partition elements
#' \item split a partition element or join one
#' \item move a single node into another partition element or into a new one
#' \item stay still
#' }
#' @param iterations (optional) integer, the number of MCMC steps, the default value is \eqn{8n^{2}\log{n}}
#' @param stepsave (optional) integer, thinning interval for the MCMC chain, indicating the number of steps between two output iterations, the default is \code{iterations/1000}
#' @param gamma (optional) tuning parameter which transforms the score by raising it to this power, 1 by default
#' @param verbose logical, if set to TRUE (default) messages about progress will be printed
#' @return an object of class \code{MCMCtrace}, which contains a list of 5 elements (each list contains \code{iterations/stepsave} elements):
#' \itemize{
#' \item incidence - contains a list of adjacency matrices of DAGs sampled at each step of MCMC
#' \item DAGscores - contains a list of scores of DAGs sampled at each step of MCMC
#' \item partitionscores - contains a list of scores of partitions of DAGs sampled at each step of MCMC
#' \item order - contains a list of permutations of the nodes in partitions of DAGs sampled at each step of MCMC
#' \item partition - contains a list of partitions of DAGs sampled at each step of MCMC
#' }
#'
#'@references Kuipers J and Moffa G (2017). Partition MCMC for inference on acyclic digraphs. Journal of the American Statistical Association 112, 282-299.
#'@references Geiger D and Heckerman D (2002). Parameter priors for directed acyclic graphical models and the characterization of several probability distributions. The Annals of Statistics 30, 1412-1440.
#'@references Heckerman D and Geiger D (1995). Learning Bayesian networks: A unification for discrete and Gaussian domains. In Eleventh Conference on Uncertainty in Artificial Intelligence, pages 274-284.
#'@references Kalisch M, Maechler M, Colombo D, Maathuis M and Buehlmann P (2012). Causal inference using graphical models with the R package pcalg. Journal of Statistical Software 47, 1-26.
#'@references Kuipers J, Moffa G and Heckerman D (2014). Addendum on the scoring of Gaussian directed acyclic graphical models. The Annals of Statistics 42, 1689-1691.
#'@examples
#'\dontrun{
#'myDAG<-pcalg::randomDAG(20, prob=0.15, lB = 0.4, uB = 2) 
#'myData<-pcalg::rmvDAG(200, myDAG) 
#'myScore<-scoreparameters(20, "bge", myData)
#'partfit<-partitionMCMC(20,myScore)
#'sample.check(20,partfit$chain$incidence,myDAG)
#'}
#'@import pcalg
#'@export
partitionMCMC<-function(n, scoreparam, startspace=NULL, blacklist=NULL,scoretable=NULL, startDAG=NULL, 
                        moveprobs=NULL, iterations=NULL,  stepsave=NULL,gamma=1,verbose=TRUE) {
  if (is.null(moveprobs)) {
    prob1start<-40/100
    prob1<-prob1start*100
    if(n>3){ prob1<-round(6*prob1*n/(n^2+10*n-24)) }
    prob1<-prob1/100
    prob2start<-99/100-prob1start
    prob2<-prob2start*100
    if(n>3){ prob2<-round(6*prob2*n/(n^2+10*n-24)) }
    prob2<-prob2/100
    moveprobs.partition<-c(prob1,prob1start-prob1,prob2start-prob2,prob2,0.01)
    moveprobs.partition<-moveprobs.partition/sum(moveprobs.partition) # normalisation
  }
  if(is.null(iterations)){
    if(n<26){
      iterations<-30000
    } else {
      iterations<-(8*n*n*log(n))-(8*n*n*log(n)) %% 1000
    }
  }
  if(is.null(stepsave)){
    stepsave<-floor(iterations/1000)
  }
  partres<-partitionMCMCplus1sample(n,param=scoreparam,startspace=startspace,blacklist=blacklist,
                                    moveprobs=moveprobs.partition,numit=iterations,DAG=startDAG,
                                    stepsave=stepsave,startorder=NULL,scoretable=scoretable,verbose=verbose,
                                    gamma=gamma)
  return(partres)
  
}


#'Structure learning with an iterative order MCMC algorithm on an expanded search space
#'
#'This function implements an iterative search for the maximum a posteriori (MAP) DAG, 
#'by means of order MCMC.  At each iteration, the current search space is expanded by 
#'allowing each node to have up to one additional parent not already included in the search space. 
#'By default the initial search space is obtained through the PC-algorithm (using the functions \code{\link[pcalg]{skeleton}} and \code{\link[pcalg]{pc}} from the `pcalg' package [Kalisch et al, 2012]).  
#'At each iteration order MCMC is employed to search for the MAP DAG.  
#'The edges in the MAP DAG are added to the initial search space to provide 
#'the search space for the next iteration.  The algorithm iterates until no 
#'further score improvements can be achieved by expanding the search space.  
#'The final search space may be used for the sampling versions of \code{\link{orderMCMC}} and \code{\link{partitionMCMC}}.
#'
#' @param n number of nodes in the Bayesian network
#' @param scoreparam an object of class \code{scoreparameters}, containing the data and scoring parameters; see constructor function \code{\link{scoreparameters}}
#' @param plus1it (optional) integer, a number of iterations of search space expansion; by default the algorithm iterates until no score improvement can be achieved by further expanding the search space
#' @param moveprobs (optional) a numerical vector of 4 values in \code{\{0,1\}} corresponding to the probabilities of the following MCMC moves in the order space:
#' \itemize{
#' \item exchanging 2 random nodes in the order
#' \item exchanging 2 adjacent nodes in the order
#' \item placing a single node elsewhere in the order
#' \item staying still
#' }
#' @param iterations (optional) integer, the number of MCMC steps, the default value is \eqn{3.5n^{2}\log{n}}
#' @param stepsave (optional) integer, thinning interval for the MCMC chain, indicating the number of steps between two output iterations, the default is \code{iterations}/1000
#' @param MAP logical, if TRUE (default) the search targets the MAP DAG (a DAG with maximum score),
#' if FALSE at each MCMC step a DAG is sampled from the order proportionally to its score; when expanding a search space when MAP=TRUE all edges from the maximum scoring DAG are added
#'  to the new space, when MAP=FALSE only edges with posterior probability higher than defined by parameter \code{posterior} are added to the search space
#' @param posterior logical, when \code{MAP} set to FALSE defines posterior probability threshold for adding the edges to the search space 
#' @param alpha (optional) numerical significance value in \code{\{0,1\}} for the conditional independence tests in the PC-stage (by default \eqn{0.4} for \eqn{n<50}, \eqn{20/n} for \eqn{n>50})
#' @param gamma (optional) tuning parameter which transforms the score by raising it to this power, 1 by default
#' @param startorder (optional) integer vector of length n, which will be used as the starting order in the MCMC algorithm, the default order is \code{c(1:n)}
#' @param startspace (optional) a square matrix, of dimensions equal to the number of nodes, which defines the search space for the order MCMC in the form of an adjacency matrix; if NULL, the skeleton obtained from the PC-algorithm will be used; if \code{startspace[i,j]} equals to 1 (0) it means that the edge from node \code{i} to node \code{j} is included (excluded) from the search space; to include an edge in both directions, both \code{startspace[i,j]} and \code{startspace[j,i]} should be 1
#' @param scoretable (optional) list of score tables which has to match startspace and addspace
#' @param addspace (optional) a square matrix, of dimensions equal to the number of nodes, which defines the edges, which are added at to the search space only at the first iteration of iterative seach and do not necessarily stay afterwards; defined in the form of an adjacency matrix;  if \code{addspace[i,j]} equals to 1 (0) it means that the edge from node \code{i} to node \code{j} is included (excluded) from the search space; to include an edge in both directions, both \code{addspace[i,j]} and \code{addspace[j,i]} should be 1
#' @param blacklist (optional) a square matrix, of dimensions equal to the number of nodes, which defines edges to exclude from the search space; if \code{blacklist[i,j]} equals to 1 it means that the edge from node \code{i} to node \code{j} is excluded from the search space
#' @param softlimit (optional) integer, limit on the size of parent sets beyond which adding undirected edges is restricted; below this
#' limit edges are added to expand the parent sets based on the undirected skeleton of the MAP DAG (or from its CPDAG, depending
#' on the parameter \code{mergecp}), above the limit only the directed edges are added from the MAP DAG;  the limit is 9 by default
#' @param hardlimit (optional) integer, limit on the size of parent sets beyond which the search space is not further expanded to prevent long runtimes; the limit is 12 by default
#' @param cpdag logical, if set to TRUE the equivalence class (CPDAG) found by the PC algorithm is used as a search
#'  space, when FALSE (default) the undirected skeleton used as a search space
#' @param mergetype defines which edges are added to the search space at each expansion iteration; if set to 
#' \itemize{
#' \item "dag", then edges from maximum scoring DAG are added;
#' \item "cpdag", then the maximum scoring DAG is first converted to the CPDAG, from which all edges are added to the search space;
#' \item "skeleton", then the maximum scoring DAG is first converted to the skeleton, from which all edges are added to the search space
#' }
#' @param verbose logical, if TRUE (default) prints messages on the progress of execution
#' @param chainout logical, if TRUE the saved MCMC steps are returned, FALSE by default
#' @param scoreout logical, if TRUE the search space from the last plus1 iterations and the corresponding score tables are returned, FALSE by default
#' @return Depends on the logical parameters \code{chainout} and \code{scoreout}. If both are FALSE (default), an oject of class \code{MCMCmax}, containing a list of 4 elements:
#' \itemize{
#' \item {DAG -} {the adjacency matrix of the DAG with maximal score}
#' \item {order -} {an order it belongs to}
#' \item {score -} {the score of the reported DAG}
#' \item {it -} {the iteration at which maximum was reached}
#' }
#' If \code{chainout} is TRUE an object of class \code{MCMCtrace} is additionally returned, contains 4 lists  (each of the 4 lists has length
#' \code{iterations/stepsave}, i.e. the number of saved MCMC steps):
#' \itemize{
#' \item incidence - contains a list of adjacency matrices of DAGs sampled at each step of MCMC
#' \item DAGscores - contains a list of scores of DAGs sampled at each step of MCMC
#' \item orderscores - contains a list of scores of orders of DAGs sampled at each step of MCMC
#' \item order - contains a list of permutations of the nodes of DAGs sampled at each step of MCMC
#' }
#' If \code{scoreout} is TRUE an object of class \code{MCMCspace} is additionally returned, contains a list of 2 elements:
#' \itemize{
#'  \item adjacency - the adjacency matrix representing the search space
#'  \item scoretable - the list of score tables corresponding to this search space
#' }
#'@references Friedman N and Koller D (2003). A Bayesian approach to structure discovery in bayesian networks. Machine Learning 50, 95-125.
#'@references Kalisch M, Maechler M, Colombo D, Maathuis M and Buehlmann P (2012). Causal inference using graphical models with the R package pcalg. Journal of Statistical Software 47, 1-26.
#'@references Geiger D and Heckerman D (2002). Parameter priors for directed acyclic graphical models and the characterization of several probability distributions. The Annals of Statistics 30, 1412-1440.
#'@references Kuipers J, Moffa G and Heckerman D (2014). Addendum on the scoring of Gaussian directed acyclic graphical models. The Annals of Statistics 42, 1689-1691.
#'@references Spirtes P, Glymour C and Scheines R (2000). Causation, Prediction, and Search, 2nd edition. The MIT Press.
#'@examples
#'\dontrun{
#'Bostonpar<-scoreparameters(14,"bge",Boston)
#'itfit<-iterativeMCMCsearch(14,Bostonpar, scoreout=TRUE)
#'sp<-itfit$space$adjacency 
#'scores<-itfit$space$scoretable
#'ordersample<-orderMCMC(14, Bostonpar, MAP=FALSE, startspace=sp, scoretable=scores)
#'}
#'@import pcalg
#'@importFrom methods new
#'@importFrom graphics lines
#'@importFrom graphics par
#'@importFrom stats cor 
#'@importFrom stats cov 
#'@importFrom stats cov.wt
#'@importFrom stats pchisq 
#'@importFrom stats runif 
#'@importFrom stats rnorm 
#'@importFrom utils data
#'@importFrom utils flush.console
#'@importFrom Rcpp evalCpp
#'@useDynLib BiDAG, .registration=TRUE
# pcalg pc
# pcalg skeleton
# pcalg shd
# pcalg dag2cpdag
#'@export
iterativeMCMCsearch<-function(n, scoreparam, plus1it=NULL, moveprobs=NULL, MAP=TRUE, posterior=0.5,
                              iterations=NULL, stepsave=NULL, softlimit=9, hardlimit=12, alpha=NULL, gamma=1, 
                              startspace=NULL, blacklist=NULL,verbose=TRUE, chainout=FALSE, scoreout=FALSE, cpdag=FALSE, 
                              mergetype="skeleton",addspace=NULL,scoretable=NULL,startorder=c(1:n)) {
  
  if (is.null(moveprobs)) {
    prob1<-99
    if(n>3){ prob1<-round(6*99*n/(n^2+10*n-24)) }
    prob1<-prob1/100
    moveprobs<-c(prob1,0.99-prob1,0.01)
    moveprobs<-moveprobs/sum(moveprobs) # normalisation
    moveprobs<-c(moveprobs[c(1,2)],0,moveprobs[3])
  }
  if(is.null(iterations)) {
    if(n<26){
      iterations<-25000
    } else {
      iterations<-(3.5*n*n*log(n))-(3.5*n*n*log(n)) %% 1000
    }
  }
  if(is.null(stepsave)) {
    stepsave<-floor(iterations/1000)
  }
  result<-iterativeMCMCplus1(n,param=scoreparam,iterations,stepsave,plus1it=plus1it, MAP=MAP, posterior=posterior,alpha=alpha,cpdag=cpdag,
                             moveprobs=moveprobs,softlimit=softlimit,hardlimit=hardlimit,
                             plus1=TRUE,startspace=startspace,blacklist=blacklist,gamma=gamma,
                             verbose=verbose, chainout=chainout,scoreout=scoreout,mergecp=mergetype,
                             addspace=addspace,scoretable=scoretable,startorder=startorder)
  return(result)
  
}

#'Calculating the BGe/BDe score of a single DAG
#'
#'This function calculates the score of a DAG defined by its adjacency matrix. 
#'Acceptable data matrices are homogeneous with all variables of the same type: 
#'continuous, binary or categorical.  The BGe score is evaluated in the case of 
#'continuous data and the BDe score is evaluated for binary and categorical variables.
#'
#' @param n number of nodes in the Bayesian network
#' @param scoreparam an object of class \code{scoreparameters}, containing the data and
#'  scoring parameters; see constructor function \code{\link{scoreparameters}}
#' @param incidence a square matrix of dimensions equal to the number of nodes, representing the adjacency matrix of a DAG;  the matrix entries are in \code{\{0,1\}} such that \code{incidence[i,j]} equals 1 if there is a directed edge from node \code{i} to node \code{j} in the DAG and 
#' \code{incidence[i,j]} equals 0 otherwise
#' @return the log of the BGe or BDe score of the DAG
#' @references Geiger D and Heckerman D (2002). Parameter priors for directed acyclic graphical models and the characterization of several probability distributions. The Annals of Statistics 30, 1412-1440.
#' @references Heckerman D and Geiger D (1995). Learning Bayesian networks: A unification for discrete and Gaussian domains. In Eleventh Conference on Uncertainty in Artificial Intelligence, pages 274-284.
#' @references Kuipers J, Moffa G and Heckerman D (2014). Addendum on the scoring of Gaussian directed acyclic graphical models. The Annals of Statistics 42, 1689-1691.
#' @examples
#' myDAG<-pcalg::randomDAG(20, prob=0.15, lB = 0.4, uB = 2) 
#' myData<-pcalg::rmvDAG(200, myDAG) 
#' adjacency<-dag2adjacencymatrix(myDAG)
#' myScore<-scoreparameters(20,"bge",myData)
#' DAGscore(20,myScore, adjacency)
#' @import pcalg
#' @export
#' 
DAGscore <- function(n,scoreparam, incidence){ 
  n<-ncol(scoreparam$data)
  P_local <- numeric(n)
  for (j in 1:n)  { #j is a node at which scoring is done
    parentnodes <- which(incidence[,j]==1)
    P_local[j]<-DAGcorescore(j,parentnodes,n,scoreparam)
  }
  #print(P_local)
  return(sum(P_local))
}
