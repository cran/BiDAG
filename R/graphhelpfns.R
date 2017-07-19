#'Deriving an adjacency matrix of a graph
#'
#'This function derives the adjacency matrix corresponding to a graph object
#'
#'@param g graph, object of class \code{\link[graph]{graphNEL}} (package `graph')
#'@return a square matrix whose dimensions are the number of nodes in the graph g, where element
#' \code{[i,j]} equals \code{1} if there is a directed edge from node \code{i} to node \code{j} in the graph \code{g},
#'  and \code{0} otherwise
#'@export
dag2adjacencymatrix<-function(g) {
  l<-length(g@edgeL)
  adj<-matrix(rep(0,l*l),nrow=l,ncol=l)
  for (i in 1:l) {
    adj[g@edgeL[[i]]$edges,i]<-1}
  return(t(adj))
}

#'Deriving an adjacency matrix of the skeleton of a graph
#'
#'This function derives the skeleton matrix corresponding to a graph object
#'
#'@param g graph, object of class \code{\link[graph]{graphNEL}} (package `graph')
#'@return a symmetric square matrix whose dimensions are the number of nodes in the graph \code{g},
#' where element \code{[i,j]} equals \code{1} if there is a directed edge from node \code{i} to node \code{j},
#'  or from node \code{j} to node \code{i}, in the graph \code{g}, and \code{0} otherwise
#'@examples 
#'myDAG<-pcalg::randomDAG(20, prob=0.15, lB = 0.4, uB = 2)
#'dag2skeletonadjacency(myDAG)
#'@export
dag2skeletonadjacency <-function(g) {
  skeletonedges<-g@edgeL
  l1<-length(skeletonedges)
  edges<-matrix(rep(0,l1^2), ncol=l1,nrow=l1)
  for (j in 1:l1)
  {
    edges[j,skeletonedges[[j]]$edges]<-1
  }
  edges<-edges|t(edges)
  edges<-ifelse(upper.tri(edges)==TRUE,edges,0)
  return(edges)
}

#'Deriving a graph from an adjacancy matrix
#'
#'This function derives a graph object corresponding to an adjacency matrix
#'
#'@param adj square adjacency matrix with elements in \code{\{0,1\}}, representing a graph
#'@param nodes (optional) labels of the nodes, \code{c(1:n)} are used by default
#'@return object of class \code{\link[graph]{graphNEL}} (package `graph'); if element \code{adj[i,j]} equals \code{1}, then there is a directed edge from node \code{i} to node \code{j} in the graph, and no edge otherwise
#'@examples 
#'adj<-matrix(rep(0,16),nrow=4)
#'adj[2,1]<-1
#'adj[1,4]<-1
#'adjacency2dag(adj)
#'@export
adjacency2dag<-function(adj,nodes=NULL) {
  l<-ncol(adj)
  if (length(nodes)==0) {
    V <- c(1:l)
    edL <- vector("list", length=l)
    names(edL) <- sapply(V,toString)
    for(i in 1:l)
      edL[[i]] <- list(edges=which(adj[i,]==1))}
  else {
    V <- nodes
    edL <- vector("list", length=l)
    names(edL) <- V
    for(i in 1:l)
      edL[[i]] <- list(edges=which(adj[i,]==1))
  }
  gR <- new("graphNEL", nodes=sapply(V,toString), edgeL=edL,edgemode="directed")
  
  return(gR)
  
}

#'Comparing two DAGs
#'
#'This function compares one (estimated) DAG to another DAG (true DAG), returning a vector of 3 values: structural Hamming distance,
#'number of true positive edges and number of false positive edges.
#'
#'@param eDAG an object of class \code{\link[graph]{graphNEL}} (package `graph'), representing the DAG which should be compared to a ground truth DAG
#'@param trueDAG an object of class \code{\link[graph]{graphNEL}} (package `graph'), representing the ground truth DAG
#'@return a vector of 3: SHD, number of true positive edges and number of false positive edges
#'@examples
#'myDAG<-pcalg::randomDAG(20, prob=0.15, lB = 0.4, uB = 2)
#'myData<-pcalg::rmvDAG(200, myDAG) 
#'myScore<-scoreparameters(20,"bge",myData)
#'\dontrun{
#'eDAG<-orderMCMC(20,myScore)
#'compareDAGs(adjacency2dag(eDAG$max$DAG),myDAG)
#'}
#'@export
compareDAGs<-function(eDAG,trueDAG) {
  skeleton1<-dag2skeletonadjacency(eDAG)
  skeleton2<-dag2skeletonadjacency(trueDAG)
  numedges1<-sum(skeleton1)
  numedges2<-sum(skeleton2)
  diff2<-skeleton2-skeleton1
  res<-vector()
  res[1]<-pcalg::shd(eDAG, trueDAG)
  res[2]<-numedges1-sum(diff2<0)
  res[3]<-sum(diff2<0)
  attr(res, "class")<-"compDAGs"
  return(res)
}

#returns the adjacency matrix corresponding to a graph object such that children are in columns!!
dagadjacencymatrix<-function(g) {
  l<-length(g@edgeL)
  adj<-matrix(rep(0,l*l),nrow=l,ncol=l)
  for (i in 1:l) {
    adj[g@edgeL[[i]]$edges,i]<-1}
  return(adj)
}

#returns a graph object which represents the skeleton of a given adjacency matrix
adjacency2skelDAG<-function(adj,nodes=NULL) {
  l<-ncol(adj)
  if (length(nodes)==0) {
    V <- c(1:l)
    edL <- vector("list", length=l)
    names(edL) <- sapply(V,toString)
    for(i in 1:l)
      edL[[i]] <- list(edges=which(adj[i,]==1))}
  else {
    V <- nodes
    edL <- vector("list", length=l)
    names(edL) <- V
    for(i in 1:l)
      edL[[i]] <- list(edges=which(adj[i,]==1))
  }
  gR <- new("graphNEL", nodes=sapply(V,toString), edgeL=edL,edgemode="undirected")
  
  return(gR)
  
}

#returns a matrix of a CPDAG corresponding to a given DAG
dagadj2cpadj<-function(adj) {
  g<-adjacency2dag(adj)
  cpg<-pcalg::dag2cpdag(g)
  return(t(dagadjacencymatrix(cpg)))
}

#returns a symmetric matrix of a skeleton corresponding to a given CPDAG
adjacency2skeleton<-function(adj) {
  skel<-1*(adj|t(adj))
  skel<-ifelse(upper.tri(skel)==TRUE,skel,0)
  return(skel)
}

