#Computing score tables
#
#This function computes the score tables which can be further used by structure learning functions
#' @param scorepar an object of class \code{scoreparameters}, containing the data and score scorepareters, see constructor function \code{\link{scoreparameters}}
#' @param alpha numerical significance value in \code{\{0,1\}} for the conditional independence tests at the PC algorithm stage (by default \eqn{0.4} for \eqn{n<50}, \eqn{20/n} for \eqn{n>50})
#' @param hardlimit integer, limit on the size of parent sets in the search space; by default 14 when MAP=TRUE and 20 when MAP=FALSE
#' @param plus1 logical, if TRUE (default) the search is performed on the extended search space
#' @param cpdag logical, if TRUE the CPDAG returned by the PC algorithm will be used as the search
#'space, if FALSE (default) the full undirected skeleton will be used as the search space
#' @param startspace (optional) a square matrix, of dimensions equal to the number of nodes, which defines the search space for the order MCMC in the form of an adjacency matrix. If NULL, the skeleton obtained from the PC-algorithm will be used. If \code{startspace[i,j]} equals to 1 (0) it means that the edge from node \code{i} to node \code{j} is included (excluded) from the search space. To include an edge in both directions, both \code{startspace[i,j]} and \code{startspace[j,i]} should be 1.
#' @param blacklist (optional) a square matrix, of dimensions equal to the number of nodes, which defines edges to exclude from the search space. If \code{blacklist[i,j]} equals to 1 it means that the edge from node \code{i} to node \code{j} is excluded from the search space.
#' @param verbose logical, if TRUE messages about the algorithm's progress will be printed, FALSE by default
#' @return Object of class \code{scorespace}, a list of three objects: 'adjacency' matrix representiong the search space, 'blacklist' used to exclude edges from the search space and 'tables' containing score quantities for each node
#' needed to run MCMC schemes 
#'@references Friedman N and Koller D (2003). A Bayesian approach to structure discovery in bayesian networks. Machine Learning 50, 95-125.
#'@examples
#'#' #find a MAP DAG with search space defined by PC and plus1 neighbourhood
#'Bostonscore<-scoreparameters("bge",Boston)
#'Bostonspace<-scorespace(Bostonscore, 0.05, 14)
#'\dontrun{
#'orderfit<-orderMCMC(Bostonscore, scoretable=Bostonspace)
#'partitionfit<-orderMCMC(Bostonscore, scoretable=Bostonspace)
#'}
#'@author Polina Suter, Jack Kuipers
#'@export
scorespace<-function(scorepar, alpha=0.05, hardlimit=14, plus1=TRUE, cpdag=TRUE, 
                     startspace=NULL, blacklist=NULL, verbose=FALSE) {
  result<-list()
 
  n<-scorepar$n
  nsmall<-scorepar$nsmall
  matsize<-ifelse(scorepar$DBN,n+nsmall,n)
  
  #defining startorder and updatenodes 
  if(!scorepar$DBN) { 
    
    if(scorepar$bgn!=0) {
      updatenodes<-c(1:n)[-scorepar$bgnodes]
    } else { 
      updatenodes<-c(1:n)
    }
    
  } else { #for DBNs startorder is defined in main.R
    updatenodes<-c(1:nsmall)
  }

  #creating blacklist objects
  if (is.null(blacklist)) {
    blacklist<-matrix(0,nrow=matsize,ncol=matsize)
  }
  diag(blacklist)<-1
  if(!is.null(scorepar$bgnodes)) {
    for(i in scorepar$bgnodes) {
      blacklist[,i]<-1
    }
  }
  
  #defining startskel
  if (is.null(startspace)){
    startspace<-definestartspace(alpha,scorepar,cpdag=cpdag,algo="pc")
  }
  startskel<-1*(startspace&!blacklist)

  
  blacklistparents<-list()
  for  (i in 1:matsize) {
    blacklistparents[[i]]<-which(blacklist[,i]==1)
  }
  
  if(verbose) {
    cat(paste("maximum parent set size is", max(apply(startskel,2,sum))),"\n")
  }
  if(max(apply(startskel,2,sum))>hardlimit) {
    stop("the size of maximal parent set is higher that the hardlimit; redifine the search space or increase the hardlimit!")
  }
  
  ptab<-listpossibleparents.PC.aliases(startskel,isgraphNEL=FALSE,n,updatenodes)
  
  if (verbose) {
    cat("skeleton ready \n")
    flush.console()
  }
  
  if (plus1==FALSE) {
    
    parenttable<-ptab$parenttable # basic parenttable without plus1 lists
    aliases<-ptab$aliases #aliases for each node since all nodes in parent tables are named as 1,2,3,4. not real parent names
    numberofparentsvec<-ptab$numberofparentsvec
    numparents<-ptab$numparents
    rowmaps<-parentsmapping(parenttable,numberofparentsvec,n,updatenodes)
    scoretable<-scorepossibleparents.alias(parenttable,aliases,n,scorepar,updatenodes,rowmaps,
                                             numparents,numberofparentsvec)
  } else {
    parenttable<-ptab$parenttable # basic parenttable without plus1 lists
    aliases<-ptab$aliases #aliases for each node since all nodes in parent tables are done as 1,2,3,4... not real parent names
    numberofparentsvec<-ptab$numberofparentsvec
    numparents<-ptab$numparents
    plus1lists<-PLUS1(matsize,aliases,updatenodes,blacklistparents)
    rowmaps<-parentsmapping(parenttable,numberofparentsvec,n,updatenodes)
    scoretable<-scorepossibleparents.PLUS1(parenttable,plus1lists,n,scorepar,updatenodes,
                                             rowmaps,numparents,numberofparentsvec) 
  }
  
  colnames(blacklist)<-rownames(blacklist)<-colnames(startskel)
  result$tables<-scoretable
  result$blacklist<-blacklist
  result$adjacency<-startskel
  
  class(result)<-"scorespace"
  
  result
  
}