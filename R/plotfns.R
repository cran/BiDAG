setMethod("plot", "MCMCmult",
          function(x, ...) {
            x<-x$chain
            if(is.null(x)) {
              stop("no saved MCMC steps found! set chainout=TRUE")
            }
            nchains<-length(x$DAGscores)
            scorevecmin<-unlist(x$DAGscores[[1]])
            scorevecprevmax<-unlist(x$DAGscores[[nchains-1]])
            minprev<-min(scorevecprevmax)
            scorevecmax<-unlist(x$DAGscores[[nchains]])
            vecl<-length(scorevecmin)
            scoremax<-max(scorevecmax)
            scoremin<-min(scorevecmin)
            scoremaxmin<-min(scorevecmax)
            par(mfrow=c(1,2))
            par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
            plot(scorevecmin,type="l",col=1,xlab="iteration",ylab="logscore",
                 ylim=c(scoremin,scoremax),main="DAG scores: all plus iterations",cex.main=1)
            for(i in 2:nchains) {
              lines(unlist(x$DAGscores[[i]]),col=i)
            }
            par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
            plot(c(scorevecprevmax,scorevecmax),type="l",col="blue",xlab="iteration",ylab="logscore",
                 ylim=c(minprev,scoremax),main="last 2 plus1 iterations",cex.main=1)
          })


setMethod("plot", "MCMCres",
          function(x, ...,burnin=0.2) {
  x<-x$chain
  if(is.null(x)) {
    stop("no saved MCMC steps found! set chainout=TRUE")
  }
  scorevec<-unlist(x$DAGscores)
  vecl<-length(scorevec)
  burnin<-ceiling(vecl*burnin)
  score20<-min(scorevec[burnin:vecl])
  scoremax<-max(scorevec)
  scoremin<-min(scorevec)
  
  par(mfrow=c(1,2))
  par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
  plot(scorevec,type="l",col="blue",xlab="iteration",ylab="logscore",
       ylim=c(scoremin,scoremax+(scoremax-scoremin)*0.02),main="DAG logscores",cex.main=1)
  par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
  plot(x=c(burnin:vecl),y=scorevec[burnin:vecl],type="l",col="blue",xlab="iteration",ylab="logscore",
       ylim=c(score20,scoremax),main="DAG logscores excluding burn-in",cex.main=1)
  
})

#'Plotting posterior probabilities of single edges
#'
#'This function plots posterior probabilities of all possible edges in the graph as a function of MCMC iterations. It can be used for convergence diagnostics of MCMC
#'sampling algorithms order MCMC and partition MCMC.
#' @param MCMCtrace an object of class MCMCres
#' @param cutoff number representing a threshold of posterior probability below which lines will not be plotted
#' @param pdag logical, when true DAGs in a sample will be first coverted to CPDAGs
#' @param onlyedges (optional) binary matrix, only edges corresponding to entries whuch equal 1 will be plotted
#'@return plots the graph which includes edges from edag and truedag, however edges which are different in edag compared to truedag are coloured according to the type of a mistake: false positive with red, false negative with dashed grey, error in direction with magenta
#'@examples
#'score100<-scoreparameters(8, "bde", Asia[1:100,])
#'orderfit100<-orderMCMC(score100,plus1=FALSE)
#'score5000<-scoreparameters(8, "bde", Asia)
#'orderfit5000<-orderMCMC(score5000,plus1=FALSE)
#'plotpedges(orderfit100)
#'plotpedges(orderfit5000)
#'@export
plotpedges<-function(MCMCtrace,cutoff=0.2,pdag=FALSE,onlyedges=NULL) {
  MCMCtrace<-MCMCtrace$chain$incidence
  if(is.null(MCMCtrace)) {
    stop("no saved MCMC steps found! try chainout=TRUE when running sampling")
  }
  cols5<-c("#cccccc","#d7b5d8","#df65b0","#dd1c77","#980043")
  cols5<-c("#cccccc","#bae4b3","#74c476","#31a354","#006d2c")
  cols5<-c("#f2f0f7","#b3cde3","#8c96c6","#8856a7","#810f7c")
  lchain<-length(MCMCtrace)
  if(pdag==TRUE) {
    MCMCtrace<-lapply(MCMCtrace,dagadj2cpadj)
  }
  countmatrix<-MCMCtrace[[1]]
  posteriors<-list()
  counter<-1
  posteriors[[1]]<-countmatrix/counter
  for (i in 2:lchain) {
    countmatrix<-countmatrix+MCMCtrace[[i]]
    counter<-counter+1
    posteriors[[i]]<-countmatrix/counter
  }
  if(!is.null(onlyedges)) {
    cutoffelems<-which(posteriors[[lchain]]>cutoff & onlyedges==1)
  } else {
    cutoffelems<-which(posteriors[[lchain]]>cutoff)
  }
  numelem<-length(cutoffelems)
  postvec<-list()
  colvec<-vector()
  k<-1
  for(i in cutoffelems){
    postvec[[k]]<-lapply(posteriors,function(x)x[i])
    colvec[k]<-defcolrange(posteriors[[lchain]][i])
    k<-k+1
  }
  plot(x=c(1:lchain),y=postvec[[1]],type="l",
       col=cols5[colvec[1]],xlab="iteration",ylab="posterior",
       main="posterior probabilities of all edges",ylim=c(0,1),cex.main=1)
  for(i in 2:numelem){
    lines(x=c(1:lchain),postvec[[i]],type="l",col=cols5[colvec[i]])
  }
}


#'Comparing posterior probabilitites of single edges based on two samples
#'
#'This function can be used to compare posterior probabilities of edges in a graph 
#'based on two samples of graphs
#'@param edgepmat1 binary matrix, representing posterior probabilities of single edges in a Bayesian network
#'@param edgepmat2 binary matrix, representing posterior probabilities of single edges in a Bayesian network
#'@param highlight numeric, defines maximum acceptable difference between posterior probabilities of an edge in two samples; points corresponding to higher differences are highlighted 
#'@param cut numeric value corresponding to a minimum posterior probabilitity which is included into calculation of squared correlation and MSE
#'@param main character string, a title for the plot
#'@param xlab character string, a title for the x-axis
#'@param ylab character string, a title for the y-axis
#'@return squared correlation and MSE of posterior probabilities higher than the value defined by the argument cut; also plots the posterior probabilitites from two samples against each other
#'@examples
#'Asiascore<-scoreparameters(8, "bde", Asia)
#'orderfit1<-orderMCMC(Asiascore,plus1=FALSE,iterations=10000)
#'orderfit2<-orderMCMC(Asiascore,plus1=FALSE,iterations=30000)
#'pedges1<-edges.posterior(orderfit1)
#'pedges2<-edges.posterior(orderfit2)
#'plotpcor(pedges1,pedges2)
#'@export
plotpcor<-function(edgepmat1,edgepmat2,highlight=0.3,cut=0.05,main="",
                   xlab="sample 1",
                   ylab="sample 2") {
  
  vec1<-as.vector(edgepmat1)
  vec2<-as.vector(edgepmat2)
  
  if(highlight<1){
    diffmat<-abs(edgepmat1-edgepmat2)
    pointstohighlight<-which(diffmat>highlight)
    if(length(pointstohighlight)>0) {
      vec1high<-vec1[pointstohighlight]
      vec2high<-vec2[pointstohighlight]
      vec1<-vec1[-pointstohighlight]
      vec2<-vec2[-pointstohighlight]
    } else {
      highlight<-1
    }
  }
  
  
  plot(x=seq(0,1,by=0.1),y=seq(0,1,by=0.1),type="l",col="blue",lty=2,
       xlab=xlab,ylab=ylab,
       xlim=c(0,1),ylim=c(0,1),main=main)
  lines(vec1,vec2,type="p",col="grey")
  if(highlight<1) lines(vec1high,vec2high,type="p",col="red")
  
  rsqindex<-intersect(which(vec1<cut),which(vec2<cut))
  Rsq<-cor(vec1[-rsqindex],vec2[-rsqindex])^2
  res<-list()
  res$MSE<-sum((vec1[-rsqindex]-vec2[-rsqindex])^2)/length(vec1[-rsqindex])
  res$R2<-Rsq
  return(res)
}


#' Plotting a DBN
#' 
#' This function can be used for plotting initial and transitional structures of a dynamic Bayesian network.
#' 
#'
#'@param DBN binary matrix (or a graph object) representing a 2-step DBN (compact or unrolled)
#'@param struct option used to determine if the initial or the transitional structure should be plotted; accaptable values are init or trans
#'@param n.dynamic number of dynamic variables in one time slice of a DBN
#'@param n.static number of static variables in one time slice of a DBN; note that for function to work correctly all static variables have to be in the first n.static columns of the matrix
#'@examples
#'plotDBN(DBNmat, "trans", n.dynamic=12,n.static=3)
#'plotDBN(DBNmat, "init", n.dynamic=12,n.static=3)
#' @export
plotDBN<-function(DBN,struct=c("init","trans"),n.dynamic,n.static){
  
  if(!is.matrix(DBN)) {
    DBN<-graph2m(DBN)
  }
  
  nodelabs<-colnames(DBN)
  statcol<-"lightgrey"
  dyn1col<-"#f7f4f9"
  dyn2col<-"#d4b9da"
  
  
  if(struct=="init") {
    
    nodelabs<-nodelabs[1:(n.dynamic+n.static)]
    
    if(n.static>0){
      legadj<-matrix(0,nrow=2,ncol=2)
      colnames(legadj)<-c("static","t=1")
      legadj[1,2]<-1
      legendG<-m2graph(legadj)
      staticnames<-nodelabs[1:n.static]
    }
    
    dynamicnames<-nodelabs[1:n.dynamic+n.static]
    
    adj<-DBN[1:(n.dynamic+n.static),1:(n.dynamic+n.static)]
    arcslist<-adjacency2edgel(adj,nodes=nodelabs)
    graph.obj = new("graphNEL", nodes = nodelabs, edgeL = arcslist,
                    edgemode = 'directed')
    
    subGList<-list()
    sg<-list()
    
    
    if(n.static!=0) {  
      twocolors<-c(statcol,dyn1col)
      nodeType <- 1 + (nodes(graph.obj) %in% dynamicnames)
      nA = makeNodeAttrs(graph.obj, fillcolor=twocolors[nodeType],shape="circle")
      sg1 = subGraph(dynamicnames, graph.obj)
      sg2 = subGraph(staticnames, graph.obj)
      
      sgL = list(list(graph=sg1, cluster = TRUE),
                 list(graph=sg2, cluster = TRUE))
      att = list(graph = list(rankdir = "TB"))
    } else {
      nA = makeNodeAttrs(graph.obj, fillcolor=dyn1col,shape="circle")
    }
    
    if(n.static>0) {
      layout(matrix(c(1,1,1,1,3,
                      1,1,1,1,2,
                      1,1,1,1,3), nrow = 3, ncol = 5, byrow = TRUE))
      plot(graph.obj, attrs = att, nodeAttrs=nA, subGList = sgL,main="DBN initial structure")
      nA.legend = makeNodeAttrs(legendG, fillcolor=c(statcol,
                                                     dyn1col))
      plot(legendG,attrs=list(graph=list(rankdir="TB"),
                              edge=list(color="white")),nodeAttrs=nA.legend,
           main="node codes")
    } else {
      plot(graph.obj, attrs = att, nodeAttrs=nA)
    }
  }
  
  if(struct=="trans") {
    if(n.static>0) {
      legadj<-matrix(0,nrow=3,ncol=3)
      legadj[1,2]<-1
      legadj[2,3]<-1
      colnames(legadj)<-c("static","t=i","t=i+1")
      legendG<-m2graph(legadj)
    } else {
      legadj<-matrix(0,nrow=2,ncol=2)
      legadj[1,2]<-1
      colnames(legadj)<-c("t=i","t=i+1")
      legendG<-m2graph(legadj)
    }
    
    adjt<-DBNcut(DBN[1:(n.static+2*n.dynamic),1:(n.static+2*n.dynamic)],n.dynamic,n.static)
    graph.obj<-m2graph(adjt)
    staticnames<-nodelabs[1:n.static]
    dyn1names<-nodelabs[1:n.dynamic+n.static]
    dyn2names<-nodelabs[1:n.dynamic+n.static+n.dynamic]
    
    
    nA = makeNodeAttrs(graph.obj, fillcolor=c(rep(statcol,n.static),
                                              rep(dyn1col,n.dynamic),
                                              rep(dyn2col,n.dynamic)),
                       shape="circle")
    sgStat = subGraph(staticnames, graph.obj)
    sgDyn1 = subGraph(dyn1names, graph.obj)
    sgDyn2 = subGraph(dyn2names, graph.obj)
    
    
    sgL = list(list(graph=sgStat, cluster = TRUE, attrs = c(rankdir="TB",rank="sink")),
               list(graph=sgDyn1, cluster = TRUE, attrs = c(rank="same")),
               list(graph=sgDyn2, cluster = TRUE, attrs = c(rank="same")))
    att = list(graph = list(rankdir = "TB", rank = ""))
    layout(matrix(c(1,1,1,1,3,
                    1,1,1,1,2,
                    1,1,1,1,3), nrow = 3, ncol = 5, byrow = TRUE))
    
    plot(graph.obj, attrs = att, nodeAttrs=nA, subGList = sgL,main="DBN transitional structure")
    
    nA.legend = makeNodeAttrs(legendG, fillcolor=c(statcol,
                                                   dyn1col,
                                                   dyn2col))
    plot(legendG,attrs=list(graph=list(rankdir="TB"),edge=list(color="white")),nodeAttrs=nA.legend,
         main="node codes")
  }
  
}

#' Plotting difference between two graphs
#' 
#' This function plots an estimated graph such that the edges which are different to the ground truth graph are highlighted. 
#' 
#'@param edag object of class graphNEL (or its adjacency matrix), representing estimated structure (not necessarily acyclic) to be compared to the ground truth graph
#'@param truedag object of class graphNEL (or its adjacency matrix), representing the ground truth structure (not necessarily acyclic)
#'@param clusters (optional) a list of nodes to be represented on the graph as clusters 
#'@return plots the graph which includes edges from edag and truedag, however edges which are different in edag compared to truedag are coloured according to the type of a mistake: false positive with red, false negative with dashed grey, error in direction with magenta
#'@examples
#'Asiascore<-scoreparameters(8,"bde",Asia)
#'Asiamap<-orderMCMC(Asiascore)
#'plotdiffs(Asiamap$max$DAG,Asiamat)
#'Asiacp<-pcalg::dag2cpdag(m2graph(Asiamat))
#'mapcp<-pcalg::dag2cpdag(m2graph(Asiamap$max$DAG))
#'plotdiffs(mapcp,Asiacp)
#'@export
plotdiffs<-function(edag,truedag,clusters=NULL) {
  
  if(!is.matrix(edag)) {
    adj<-graph2m(edag)
  } else {
    adj<-edag
  }
  if(!is.matrix(truedag)) {
    adjt<-graph2m(truedag)
  } else {
    adjt<-truedag
  }
  
  nodelabs<-colnames(adj)
  if(is.null(nodelabs)) {
    nodelabs<-paste("v",1:ncol(adj),sep="")
  }
  jointmat<-1*(adj|adjt)
  n<-nrow(adj)
  
  FPlist<-NULL
  FNlist<-NULL
  EDlist<-NULL
  MissDlist<-NULL #missing directions
  ExtrDlist<-NULL #extra directions
  
  BiFP<-NULL
  BiFN<-NULL
  
  
  diffmat<-matrix(0,nrow=nrow(jointmat),ncol=ncol(jointmat))
  for(i in 1:n) {
    for(j in 1:n) {
      bi1<-adj[i,j]+adj[j,i]
      bi2<-adjt[i,j]+adjt[j,i]
      if(bi1==2 | bi2==2) {
        if(bi1!=bi2) {
          if(bi2==2 & bi1==0) { #FN
            diffmat[i,j]<-3
            diffmat[j,i]<-3
            FNlist<-rbind(FNlist,c(nodelabs[i],nodelabs[j]))
            FNlist<-rbind(FNlist,c(nodelabs[j],nodelabs[i]))
          } else if (bi2==2 & bi1==1) { #ED FN
            diffmat[i,j]<-4
            diffmat[j,i]<-4
            MissDlist<-rbind(MissDlist,c(nodelabs[i],nodelabs[j]))
            MissDlist<-rbind(MissDlist,c(nodelabs[j],nodelabs[j]))
          } else if (bi2==1 & bi1==2) { #ED FP
            diffmat[i,j]<-4
            diffmat[j,i]<-4
            ExtrDlist<-rbind(ExtrDlist,c(nodelabs[i],nodelabs[j]))
            ExtrDlist<-rbind(ExtrDlist,c(nodelabs[j],nodelabs[j]))
            
          } else { #FP
            diffmat[i,j]<-2
            diffmat[j,i]<-2
            FPlist<-rbind(FPlist,c(nodelabs[i],nodelabs[j]))
            FPlist<-rbind(FPlist,c(nodelabs[j],nodelabs[i]))
          }
        }
      } else {
        if(adj[i,j]!=adjt[i,j]){
          if(adj[j,i]!=adjt[j,i]) {#ED
            if(adj[i,j]==1) {
              diffmat[i,j]<-4 
              jointmat[j,i]<-0
              EDlist<-rbind(EDlist,c(nodelabs[i],nodelabs[j]))
            } else {
              diffmat[j,i]<-4
              jointmat[i,j]<-0
              EDlist<-rbind(EDlist,c(nodelabs[j],nodelabs[i]))
            }
          } else if (adj[i,j]==1) { #FP
            diffmat[i,j]<-2
            FPlist<-rbind(FPlist,c(nodelabs[i],nodelabs[j]))
          } else {#FN
            diffmat[i,j]<-3
            FNlist<-rbind(FNlist,c(nodelabs[i],nodelabs[j]))
            
          }
        }
      }
    }
  }
  
  
  jointarcs<-adjacency2edgel(jointmat,nodes=nodelabs)
  graph.obj = new("graphNEL", nodes = nodelabs, edgeL = jointarcs,
                  edgemode = 'directed')
  if(is.null(clusters)){
    subGList = NULL } else {
      numsubg<-length(clusters)
      sg<-list()
      subGList<-list()
      for(i in 1:numsubg) {
        sg[[i]] <-subGraph(clusters[[i]],  graph.obj)
        subGList[[i]]<-list(graph = sg[[i]], cluster = TRUE)
      }
    }
  graph.plot = Rgraphviz::layoutGraph(graph.obj,subGList = subGList)
  if(!is.null(FPlist)) {
    FP<-apply(FPlist, 1, paste, collapse = "~")
    graph::edgeRenderInfo(graph.plot)[["col"]][FP] = "red"
  }
  if(!is.null(FNlist)) {
    FN<-apply(FNlist, 1, paste, collapse = "~")
    graph::edgeRenderInfo(graph.plot)[["col"]][FN] = "grey"
    graph::edgeRenderInfo(graph.plot)[["lty"]][FN] = "dashed"
  }
  if(!is.null(EDlist)) {
    ED<-apply(EDlist, 1, paste, collapse = "~")
    graph::edgeRenderInfo(graph.plot)[["col"]][ED] = "blue"
  }
  if(!is.null(BiFP)) {
    BiFP<-apply(BiFP, 1, paste, collapse = "~")
    graph::edgeRenderInfo(graph.plot)[["col"]][BiFP] = "red"
  }
  if(!is.null(BiFN)) {
    BiFN<-apply(BiFN, 1, paste, collapse = "~")
    graph::edgeRenderInfo(graph.plot)[["col"]][BiFN] = "grey"
    graph::edgeRenderInfo(graph.plot)[["lty"]][BiFN] = "dashed"
    
  }
  if(!is.null(MissDlist)) {
    MissD<-apply(MissDlist, 1, paste, collapse = "~")
    graph::edgeRenderInfo(graph.plot)[["col"]][MissD] = "blue"
    graph::edgeRenderInfo(graph.plot)[["lty"]][MissD] = "solid"
    
  }
  if(!is.null(ExtrDlist)) {
    ExtrD<-apply(ExtrDlist, 1, paste, collapse = "~")
    graph::edgeRenderInfo(graph.plot)[["col"]][ExtrD] = "blue"
    graph::edgeRenderInfo(graph.plot)[["lty"]][ExtrD] = "solid"
    
  }
  
  #par(mfrow=c(1,2))
  layout(matrix(c(1,1,2), nrow = 1, ncol = 3, byrow = TRUE))
  Rgraphviz::renderGraph(graph.plot)
  #par(xpd=TRUE)
  par(mar = c(0,0,0,0))
  plot(1:10,1:10,type="n", axes = FALSE, xlab = "", ylab = "")
  op <- par(cex = 1.7)
  legend("center",legend=c("true\npositive",
                           "false\npositive",
                           "false\nnegative",
                           "error in\ndirection"),lty=c(1,1,2,1),
         col=c("black","red","grey","blue"),bty="n")
}


#' Plotting difference between two DBNs
#' 
#' This function plots an estimated DBN such that the edges which are different to the ground truth DBN are highlighted. 
#' 
#'@param eDBN object of class graphNEL (or its adjacency matrix), representing estimated structure (not necessarily acyclic) to be compared to the ground truth graph
#'@param trueDBN object of class graphNEL (or its adjacency matrix), representing the ground truth structure (not necessarily acyclic)
#'@param struct option used to determine if the initial or the transitional structure should be plotted; accaptable values are init or trans
#'@param n.dynamic number of dynamic variables in one time slice of a DBN
#'@param n.static number of static variables in one time slice of a DBN; note that for function to work correctly all static variables have to be in the first n.static columns of the matrix
#'@return plots the graph which includes edges from edag and truedag, however edges which are different in edag compared to truedag are coloured according to the type of a mistake: false positive with red, false negative with dashed grey, error in direction with magenta
#'@examples
#'dbnscore<-scoreparameters(15,"bge",DBNdata,dbnpar = list(samestruct=TRUE, slices=5),
#'DBN=TRUE,bgnodes=c(1,2,3))
#'\dontrun{
#'orderDBNfit<-iterativeMCMC(dbnscore,chainout = TRUE, mergetype = "skeleton",scoreout=TRUE,alpha=0.4)
#'plotdiffs.DBN(orderDBNfit$max$DAG,DBNmat,struct="trans",n.dynamic=12,n.static=3)
#'plotdiffs.DBN(orderDBNfit$max$DAG,DBNmat,struct="init",n.dynamic=12,n.static=3)
#'}
#'@export
plotdiffs.DBN<-function(eDBN,trueDBN,struct=c("init","trans"),n.dynamic,n.static=0) {

  if(!is.matrix(eDBN)) {
    adj<-graph2m(eDBN)
  } else {
    adj<-eDBN
  }
  if(!is.matrix(trueDBN)) {
    adjt<-graph2m(trueDBN)
  } else {
    adjt<-trueDBN
  }
  
  nodelabs<-colnames(adj)
  statcol<-"lightgrey"
  dyn1col<-"#f7f4f9"
  dyn2col<-"#d4b9da"
  
  if(struct=="init") {
    n<-n.static+n.dynamic
    nodelabs<-nodelabs[1:(n.dynamic+n.static)]
    
    if(n.static>0){
      legadj<-matrix(0,nrow=2,ncol=2)
      colnames(legadj)<-c("static","t=1")
      legadj[1,2]<-1
      legendG<-m2graph(legadj)
      staticnames<-nodelabs[1:n.static]
      legcol<-c(statcol,dyn1col)
      names(legcol)<-c("static","t=1")
    }
    
    dynamicnames<-nodelabs[1:n.dynamic+n.static]
    
    adj<-adj[1:(n.dynamic+n.static),1:(n.dynamic+n.static)]
    adjt<-adjt[1:(n.dynamic+n.static),1:(n.dynamic+n.static)]
    jointmat<-1*(adj|adjt)
    
    #define edges with wrong directions
    FPlist<-NULL
    FNlist<-NULL
    EDlist<-NULL
    MissDlist<-NULL #missing directions
    ExtrDlist<-NULL #extra directions
    
    BiFP<-NULL
    BiFN<-NULL
    
    diffmat<-matrix(0,nrow=nrow(jointmat),ncol=ncol(jointmat))
    for(i in 1:n) {
      for(j in 1:n) {
        bi1<-adj[i,j]+adj[j,i]
        bi2<-adjt[i,j]+adjt[j,i]
        if(bi1==2 | bi2==2) {
          if(bi1!=bi2) {
            if(bi2==2 & bi1==0) { #FN
              diffmat[i,j]<-3
              diffmat[j,i]<-3
              FNlist<-rbind(FNlist,c(nodelabs[i],nodelabs[j]))
              FNlist<-rbind(FNlist,c(nodelabs[j],nodelabs[i]))
            } else if (bi2==2 & bi1==1) { #ED FN
              diffmat[i,j]<-4
              diffmat[j,i]<-4
              MissDlist<-rbind(MissDlist,c(nodelabs[i],nodelabs[j]))
              MissDlist<-rbind(MissDlist,c(nodelabs[j],nodelabs[j]))
            } else if (bi2==1 & bi1==2) { #ED FP
              diffmat[i,j]<-4
              diffmat[j,i]<-4
              ExtrDlist<-rbind(ExtrDlist,c(nodelabs[i],nodelabs[j]))
              ExtrDlist<-rbind(ExtrDlist,c(nodelabs[j],nodelabs[j]))
              
            } else { #FP
              diffmat[i,j]<-2
              diffmat[j,i]<-2
              FPlist<-rbind(FPlist,c(nodelabs[i],nodelabs[j]))
              FPlist<-rbind(FPlist,c(nodelabs[j],nodelabs[i]))
            }
          }
        } else {
          if(adj[i,j]!=adjt[i,j]){
            if(adj[j,i]!=adjt[j,i]) {#ED
              if(adj[i,j]==1) {
                diffmat[i,j]<-4 
                jointmat[j,i]<-0
                EDlist<-rbind(EDlist,c(nodelabs[i],nodelabs[j]))
              } else {
                diffmat[j,i]<-4
                jointmat[i,j]<-0
                EDlist<-rbind(EDlist,c(nodelabs[j],nodelabs[i]))
              }
            } else if (adj[i,j]==1) { #FP
              diffmat[i,j]<-2
              FPlist<-rbind(FPlist,c(nodelabs[i],nodelabs[j]))
            } else {#FN
              diffmat[i,j]<-3
              FNlist<-rbind(FNlist,c(nodelabs[i],nodelabs[j]))
              
            }
          }
        }
      }
    }
    
    
    jointarcs<-adjacency2edgel(jointmat,nodes=nodelabs)
    graph.obj = new("graphNEL", nodes = nodelabs, edgeL = jointarcs,
                    edgemode = 'directed')
  
    if(n.static!=0) {  
      sg1 = subGraph(dynamicnames, graph.obj)
      sg2 = subGraph(staticnames, graph.obj)
      sgL = list(list(graph=sg1, cluster = TRUE),
                 list(graph=sg2, cluster = TRUE))
      colvector<-c(rep(statcol,n.static),
                   rep(dyn1col,n.dynamic))
      names(colvector)<-nodelabs
    } else {
      colvector<-c(rep(dyn1col,n.dynamic))
      names(colvector)<-nodelabs
    }
    
    graph.plot = Rgraphviz::layoutGraph(graph.obj, subGList = sgL)
    graph::nodeRenderInfo(graph.plot)<-list(fill=colvector,shape="circle")
    if(!is.null(FPlist)) {
      FP<-apply(FPlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][FP] = "red"
    }
    if(!is.null(FNlist)) {
      FN<-apply(FNlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][FN] = "grey"
      graph::edgeRenderInfo(graph.plot)[["lty"]][FN] = "dashed"
    }
    if(!is.null(EDlist)) {
      ED<-apply(EDlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][ED] = "blue"
    }
    if(!is.null(BiFP)) {
      BiFP<-apply(BiFP, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][BiFP] = "red"
    }
    if(!is.null(BiFN)) {
      BiFN<-apply(BiFN, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][BiFN] = "grey"
      graph::edgeRenderInfo(graph.plot)[["lty"]][BiFN] = "dashed"
      
    }
    if(!is.null(MissDlist)) {
      MissD<-apply(MissDlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][MissD] = "blue"
      graph::edgeRenderInfo(graph.plot)[["lty"]][MissD] = "solid"
      
    }
    if(!is.null(ExtrDlist)) {
      ExtrD<-apply(ExtrDlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][ExtrD] = "blue"
      graph::edgeRenderInfo(graph.plot)[["lty"]][ExtrD] = "solid"
      
    }
    
    graph.par(list(graph=list(main="Comparison of DBNs initial structures")))
    
    if(n.static>0) {
      layout(matrix(c(1,1,1,1,3,
                      1,1,1,1,2,
                      1,1,1,1,3), nrow = 3, ncol = 5, byrow = TRUE))
      Rgraphviz::renderGraph(graph.plot)
      plot(legendG,attrs=list(graph=list(rankdir="TB"),
                              edge=list(color="white")),nodeAttrs=list(fillcolor=legcol),
           main="node codes")
    } else {
      Rgraphviz::renderGraph(graph.plot)
    }
  }
  
  if(struct=="trans") {
    
    n<-n.static+2*n.dynamic
    
    if(n.static>0) {
      legadj<-matrix(0,nrow=3,ncol=3)
      legadj[1,2]<-1
      legadj[2,3]<-1
      colnames(legadj)<-c("static","t=i","t=i+1")
      legendG<-m2graph(legadj)
      legcol<-c(statcol,dyn1col,dyn2col)
      names(legcol)<-c("static","t=i","t=i+1")

      colvector = c(rep(statcol,n.static), rep(dyn1col,n.dynamic),rep(dyn2col,n.dynamic))
      names(colvector)<-nodelabs
    } else {
      legadj<-matrix(0,nrow=2,ncol=2)
      legadj[1,2]<-1
      colnames(legadj)<-c("t=i","t=i+1")
      legendG<-m2graph(legadj)
      legcol<-c(dyn1col,dyn2col)
      names(legcol)<-c("t=i","t=i+1")
      colvector = c(rep(dyn1col,n.dynamic),rep(dyn2col,n.dynamic))
      names(colvector)<-nodelabs
    }
    
    adj<-DBNcut(adj[1:(n.static+2*n.dynamic),1:(n.static+2*n.dynamic)],n.dynamic,n.static)
    adjt<-DBNcut(adjt[1:(n.static+2*n.dynamic),1:(n.static+2*n.dynamic)],n.dynamic,n.static)
    jointmat<-1*(adj|adjt)
    
    arcslist<-adjacency2edgel(jointmat,nodes=nodelabs)
    
    graph.obj = new("graphNEL", nodes = nodelabs, edgeL = arcslist,
                    edgemode = 'directed')
    
    staticnames<-nodelabs[1:n.static]
    dyn1names<-nodelabs[1:n.dynamic+n.static]
    dyn2names<-nodelabs[1:n.dynamic+n.static+n.dynamic]
    
    sgStat = subGraph(staticnames, graph.obj)
    sgDyn1 = subGraph(dyn1names, graph.obj)
    sgDyn2 = subGraph(dyn2names, graph.obj)
    
    sgL = list(list(graph=sgStat, cluster = TRUE, attrs = c(rankdir="TB",rank="sink")),
               list(graph=sgDyn1, cluster = TRUE, attrs = c(rank="same")),
               list(graph=sgDyn2, cluster = TRUE, attrs = c(rank="same")))
    
    
    graph.plot = Rgraphviz::layoutGraph(graph.obj, subGList = sgL)
    graph::nodeRenderInfo(graph.plot)<-list(fill=colvector,shape="circle")

    
    #define edges with wrong directions
    FPlist<-NULL
    FNlist<-NULL
    EDlist<-NULL
    MissDlist<-NULL #missing directions
    ExtrDlist<-NULL #extra directions
    
    BiFP<-NULL
    BiFN<-NULL


    diffmat<-matrix(0,nrow=nrow(jointmat),ncol=ncol(jointmat))
    for(i in 1:n) {
      for(j in 1:n) {
        bi1<-adj[i,j]+adj[j,i]
        bi2<-adjt[i,j]+adjt[j,i]
        if(bi1==2 | bi2==2) {
          if(bi1!=bi2) {
            if(bi2==2 & bi1==0) { #FN
              diffmat[i,j]<-3
              diffmat[j,i]<-3
              FNlist<-rbind(FNlist,c(nodelabs[i],nodelabs[j]))
              FNlist<-rbind(FNlist,c(nodelabs[j],nodelabs[i]))
            } else if (bi2==2 & bi1==1) { #ED FN
              diffmat[i,j]<-4
              diffmat[j,i]<-4
              MissDlist<-rbind(MissDlist,c(nodelabs[i],nodelabs[j]))
              MissDlist<-rbind(MissDlist,c(nodelabs[j],nodelabs[j]))
            } else if (bi2==1 & bi1==2) { #ED FP
              diffmat[i,j]<-4
              diffmat[j,i]<-4
              ExtrDlist<-rbind(ExtrDlist,c(nodelabs[i],nodelabs[j]))
              ExtrDlist<-rbind(ExtrDlist,c(nodelabs[j],nodelabs[j]))
              
            } else { #FP
              diffmat[i,j]<-2
              diffmat[j,i]<-2
              FPlist<-rbind(FPlist,c(nodelabs[i],nodelabs[j]))
              FPlist<-rbind(FPlist,c(nodelabs[j],nodelabs[i]))
            }
          }
        } else {
          if(adj[i,j]!=adjt[i,j]){
            if(adj[j,i]!=adjt[j,i]) {#ED
              if(adj[i,j]==1) {
                diffmat[i,j]<-4 
                jointmat[j,i]<-0
                EDlist<-rbind(EDlist,c(nodelabs[i],nodelabs[j]))
              } else {
                diffmat[j,i]<-4
                jointmat[i,j]<-0
                EDlist<-rbind(EDlist,c(nodelabs[j],nodelabs[i]))
              }
            } else if (adj[i,j]==1) { #FP
              diffmat[i,j]<-2
              FPlist<-rbind(FPlist,c(nodelabs[i],nodelabs[j]))
            } else {#FN
              diffmat[i,j]<-3
              FNlist<-rbind(FNlist,c(nodelabs[i],nodelabs[j]))
              
            }
          }
        }
      }
    }
    

    if(!is.null(FPlist)) {
      FP<-apply(FPlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][FP] = "red"
    }
    if(!is.null(FNlist)) {
      FN<-apply(FNlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][FN] = "grey"
      graph::edgeRenderInfo(graph.plot)[["lty"]][FN] = "dashed"
    }
    if(!is.null(EDlist)) {
      ED<-apply(EDlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][ED] = "blue"
    }
    if(!is.null(BiFP)) {
      BiFP<-apply(BiFP, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][BiFP] = "red"
    }
    if(!is.null(BiFN)) {
      BiFN<-apply(BiFN, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][BiFN] = "grey"
      graph::edgeRenderInfo(graph.plot)[["lty"]][BiFN] = "dashed"
      
    }
    if(!is.null(MissDlist)) {
      MissD<-apply(MissDlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][MissD] = "blue"
      graph::edgeRenderInfo(graph.plot)[["lty"]][MissD] = "solid"
      
    }
    if(!is.null(ExtrDlist)) {
      ExtrD<-apply(ExtrDlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][ExtrD] = "blue"
      graph::edgeRenderInfo(graph.plot)[["lty"]][ExtrD] = "solid"
      
    }
    
    graph.par(list(graph=list(main="Comparison of DBNs transition structures")))
    layout(matrix(c(1,1,1,1,3,
                    1,1,1,1,2,
                    1,1,1,1,3), nrow = 3, ncol = 5, byrow = TRUE))
    
    Rgraphviz::renderGraph(graph.plot)
    plot(legendG,attrs=list(graph=list(rankdir="TB"),edge=list(color="white")),
         nodeAttrs=list(fillcolor=legcol),main="node codes")
  }

}

