#' Plotting a DBN
#' 
#' This function can be used for plotting initial and transition structures of a dynamic Bayesian network.
#'
#'@param DBN binary matrix (or a graph object) representing a 2-step DBN (compact or unrolled)
#'@param struct option used to determine if the initial or the transition structure should be plotted; acceptable values are init or trans
#'@param b number of static variables in the DBN, 0 by default; note that for function to work correctly all static variables have to be in the first b columns of the matrix
#'@param shape string, defining the shape of the box around each node; possible values are circle, ellipse, box
#'@param ... optional parameters passed to \code{Rgraphviz} plotting functions e.g. \code{main}, \code{fontsize}
#'@return plots the DBN defined by the adjacency matrix 'DBN' and number of static and dynamic variables. When 'struct' equals "trans" the transition structure is plotted,
#'otherwise initial structure is plotted
#'@examples
#'plotDBN(DBNmat, "init", b=3)
#'plotDBN(DBNmat, "trans", b=3)
#'
#' @author Polina Suter
#' @export
plotDBN<-function(DBN,struct=c("init","trans"),b=0,shape="circle",...){
  dyn<-(ncol(DBN)-b)/2
  
  old.par<-par(no.readonly = TRUE)
  oldgraphpar<-graph.par()
  on.exit(par(old.par))
  on.exit(graph.par(oldgraphpar),add=TRUE)
  
  a<-d<-1.2
  c<-12
  
  if(is.matrix(DBN)) {
    DBN<-DBN
  } else if (is(DBN,"graphNEL")) {
    DBN<-graph2m(DBN)
  } else {
    DBN<-as.matrix(DBN)
  }
  
  nodelabs<-colnames(DBN)
  statcol<-"lightgrey"
  dyn1col<-"#f7f4f9"
  dyn2col<-"#d4b9da"
  
  
  if(struct=="init") {
    
    shapevec <- rep(shape, dyn+b) 
    nodelabs<-nodelabs[1:(dyn+b)]
    
    if(b>0){
      legadj<-matrix(0,nrow=2,ncol=2)
      colnames(legadj)<-c("stat","1")
      legadj[1,2]<-1
      legendG<-m2graph(legadj)
      staticnames<-nodelabs[1:b]
      legcol<-c(statcol,dyn1col)
    } else {
      staticnames<-c()
    }
    
    dynamicnames<-nodelabs[1:dyn+b]
    
    adj<-DBN[1:(dyn+b),1:(dyn+b)]
    arcslist<-adjacency2edgel(adj,nodes=nodelabs)
    graph.obj = new("graphNEL", nodes = nodelabs, edgeL = arcslist,
                    edgemode = 'directed')
    subGList<-list()
    sg<-list()
    
 
    
    graph.par(list(nodes=list(col=dyn1col, lty="solid", lwd=1, ...),graph=list(...,cex.main=1.5)))
    if(b!=0) {  
      sg1 = subGraph(dynamicnames, graph.obj)
      sg2 = subGraph(staticnames, graph.obj)
      sgL = list(list(graph=sg1, cluster = TRUE),
                 list(graph=sg2, cluster = TRUE))
      graph.obj <- Rgraphviz::layoutGraph(graph.obj, subGList= sgL,nodeAttrs = list(shape = shapevec))
      graph::nodeRenderInfo(graph.obj)[["fill"]][staticnames] = statcol
      graph::nodeRenderInfo(graph.obj)[["fill"]][dynamicnames] = dyn1col
      graph::nodeRenderInfo(graph.obj)[["shape"]][c(staticnames,dynamicnames)] = shape
    } else {
      graph.obj <- Rgraphviz::layoutGraph(graph.obj,nodeAttrs = list(shape = shapevec))
    }
    
    if(b>0) {
      layout(matrix(c(1,1,1,1,1,3,3,
                      1,1,1,1,1,2,2,
                      1,1,1,1,1,3,3), nrow = 3, ncol = 7, byrow = TRUE))
      
      legendG <- Rgraphviz::layoutGraph(legendG,nodeAttrs = list(shape = shapevec))
      graph::nodeRenderInfo(legendG)[["fill"]]["stat"] = statcol
      graph::edgeRenderInfo(legendG)[["lwd"]]["stat~1"] = 0
      graph::nodeRenderInfo(legendG)[["shape"]][c("stat","1")] = shape
      Rgraphviz::renderGraph(graph.obj,nodeAttrs = list(shape = shapevec))
      graph.par(list(graph=list(main="nodes(t):")))
      Rgraphviz::renderGraph(legendG,nodeAttrs = list(shape = shapevec))
      
    } else {
      graph.par(list(nodes=list(col=dyn1col, lty="solid", lwd=1, ...)))
      Rgraphviz::renderGraph(graph.obj,nodeAttrs = list(shape = shape))
    }
  }
  
  if(struct=="trans") {
    shapevec = rep(shape, 2*dyn+b) #added
    if(b>0) {
      legadj<-matrix(0,nrow=3,ncol=3)
      legadj[1,2]<-1
      legadj[2,3]<-1
      colnames(legadj)<-c("stat","i","i+1")
      legendG<-m2graph(legadj)
      legcol<-c(statcol,dyn1col,dyn2col)
      names(legcol)<-c("stat","i","i+1")
    } else {
      legadj<-matrix(0,nrow=2,ncol=2)
      legadj[1,2]<-1
      colnames(legadj)<-c("i","i+1")
      legendG<-m2graph(legadj)
      staticnames<-c()
    }
    
    adjt<-DBNcut(DBN[1:(b+2*dyn),1:(b+2*dyn)],dyn,b)
    graph.obj<-m2graph(adjt)
    
    dyn1names<-nodelabs[1:dyn+b]
    dyn2names<-nodelabs[1:dyn+b+dyn]
    sgDyn1 = subGraph(dyn1names, graph.obj)
    sgDyn2 = subGraph(dyn2names, graph.obj)
    
    if(b>0) {
      staticnames<-nodelabs[1:b]
      sgStat = subGraph(staticnames, graph.obj)
      sgL = list(list(graph=sgStat, cluster = TRUE),
                 list(graph=sgDyn1, cluster = TRUE),
                 list(graph=sgDyn2, cluster = TRUE))
    } else {
      sgL = list(list(graph=sgDyn1, cluster = TRUE),
                 list(graph=sgDyn2, cluster = TRUE))
    }
    
    graph.par(list(nodes=list(col=dyn1col, lty="solid", lwd=1, ...),graph=list(...)))
    graph.obj <- Rgraphviz::layoutGraph(graph.obj, subGList= sgL,nodeAttrs = list(shape = shape))
    
    if(b>0) graph::nodeRenderInfo(graph.obj)[["fill"]][staticnames] = statcol
    graph::nodeRenderInfo(graph.obj)[["fill"]][dyn1names] = dyn1col
    graph::nodeRenderInfo(graph.obj)[["fill"]][dyn2names] = dyn2col
    graph::nodeRenderInfo(graph.obj)[["shape"]][c(staticnames,dyn1names,dyn2names)]<-shape
    
    layout(matrix(c(1,1,1,1,1,3,3,
                    1,1,1,1,1,2,2,
                    1,1,1,1,1,3,3), nrow = 3, ncol = 7, byrow = TRUE))
    Rgraphviz::renderGraph(graph.obj,nodeAttrs = list(shape = shape))
    
    legendG <- Rgraphviz::layoutGraph(legendG)
    if(b>0) {
      graph::nodeRenderInfo(legendG)[["shape"]][c("stat","i","i+1")] <- shape
      graph::nodeRenderInfo(legendG)[["fill"]]["stat"] = statcol
      graph::nodeRenderInfo(legendG)[["fill"]]["i"] = dyn1col
      graph::nodeRenderInfo(legendG)[["fill"]]["i+1"] = dyn2col
      graph::edgeRenderInfo(legendG)[["lwd"]]["stat~i"] = 0
      graph::edgeRenderInfo(legendG)[["lwd"]]["i~i+1"] = 0
    } else {
      graph::nodeRenderInfo(legendG)[["shape"]][c("i","i+1")] <- shape
      graph::nodeRenderInfo(legendG)[["fill"]]["i"] = dyn1col
      graph::nodeRenderInfo(legendG)[["fill"]]["i+1"] = dyn2col
      graph::edgeRenderInfo(legendG)[["lwd"]]["i~i+1"] = 0
    }
    #plot graph
    graph.par(list(graph=list(main="nodes(t):",cex.main=1.8),nodes=list(fontsize=16)))
    #plot legend
    Rgraphviz::renderGraph(legendG,nodeAttrs = list(shape = shapevec))
  }
}
assigncolor<-function(nit,ncol) {
  colind<-1:(ncol-1)
  ncolsmall<-ncol-1
  nitsmall<-nit-2
  colvec<-vector()
  
  if (nit<=ncol+1) {
    if(nit<3){
      return(rep(ncol,nit))
    } else {
      return(c(tail(colind,nitsmall),ncol,ncol))
    }
  } else {
    rmndr <- nitsmall %%  ncolsmall
    div<-nitsmall %/% ncolsmall
    if(rmndr==0) firstadd<-ncol else firstadd<-tail(colind,rmndr)[1]
    for(i in 1:ncolsmall) {
      if(i>=firstadd) {
        colvec<-c(colvec,rep(colind[i],div+1))
      } else {
        colvec<-c(colvec,rep(colind[i],div))
      }
    }
  }
  return(c(colvec,ncol,ncol))
}

#' Plotting difference between two DBNs
#' 
#' This function plots an estimated DBN such that the edges which are different to the ground truth DBN are highlighted. 
#' 
#'@param eDBN object of class graphNEL (or its adjacency matrix), representing estimated structure (not necessarily acyclic) to be compared to the ground truth graph
#'@param trueDBN object of class graphNEL (or its adjacency matrix), representing the ground truth structure (not necessarily acyclic)
#'@param struct option used to determine if the initial or the transition structure should be plotted; accaptable values are init or trans
#'@param b number of static variables in one time slice of a DBN; note that for function to work correctly all static variables have to be in the first b columns of the matrix
#'@param showcl logical, when TRUE (default) nodes are shown in clusters according to the time slice the belong to
#'@param orientation orientation of the graph layout, possible options are 'TB' (top-bottom) and 'LR' (left-right)
#'@param ... optional parameters passed to \code{Rgraphviz} plotting functions e.g. \code{main}, \code{fontsize}
#'@return plots the graph highlights differences between 'eDBN' (estimated DBN) and 'trueDBN' (ground truth); edges which are different in 'eDBN' compared to 'trueDBN' are coloured according to the type of a difference: false-positive, false-negative and error in direction.
#'@examples
#'dbnscore<-scoreparameters("bge",DBNdata,
#'dbnpar = list(samestruct=TRUE, slices=5, b=3),
#'DBN=TRUE)
#'\dontrun{
#'orderDBNfit<-learnBN(dbnscore,algorithm="order")
#'iterDBNfit<-learnBN(dbnscore,algorithm="orderIter")
#'plotdiffsDBN(getDAG(orderDBNfit),DBNmat,struct="trans",b=3)
#'plotdiffsDBN(getDAG(iterDBNfit),DBNmat,struct="trans",b=3)
#'}
#'@export
#'@author Polina Suter
plotdiffsDBN<-function(eDBN,trueDBN,struct=c("init","trans"),b=0, showcl=TRUE, orientation="TB",...) {
  
  old.par<-par(no.readonly = TRUE)
  oldgraphpar<-graph.par()
  on.exit(par(old.par))
  on.exit(graph.par(oldgraphpar),add=TRUE)
  shape<-"circle"
  a<-d<-1.2
  c<-12
  
  if(is.matrix(eDBN)) {
    adj<-eDBN
  } else if (is(eDBN,"graphNEL")) {
    adj<-graph2m(eDBN)
  } else {
    adj<-as.matrix(eDBN)
  }
  
  dyn<-(ncol(adj)-b)/2
  
  if(!is.matrix(trueDBN)) {
    if (is(trueDBN,"graphNEL")) {
      adjt<-graph2m(trueDBN)
    } else {
      adjt<-as.matrix(trueDBN)
    }
  } else {
    adjt<-trueDBN
  }
  
  nodelabs<-colnames(adj)
  statcol<-"lightgrey"
  dyn1col<-"#f7f4f9"
  dyn2col<-"#d4b9da"
  
  if(struct=="init") {
    n<-b+dyn
    nodelabs<-nodelabs[1:(dyn+b)]
    
    if(b>0){
      legadj<-matrix(0,nrow=2,ncol=2)
      colnames(legadj)<-c("stat","1")
      legadj[1,2]<-1
      legendG<-m2graph(legadj)
      staticnames<-nodelabs[1:b]
      legcol<-c(statcol,dyn1col)
      colvector = c(rep(statcol,b),rep(dyn1col,dyn))
      names(colvector)<-nodelabs
    } else {
      legcol<-c(dyn1col)
      colvector = c(rep(dyn1col,dyn))
      names(colvector)<-nodelabs
    }
    
    dynamicnames<-nodelabs[1:dyn+b]
    
    adj<-adj[1:(dyn+b),1:(dyn+b)]
    adjt<-adjt[1:(dyn+b),1:(dyn+b)]
    jointmat<-1*(adj|adjt)
    
    #define edges with wrong directions
    FPlist<-NULL
    FNlist<-NULL
    EDlist<-NULL
    MissDlist<-NULL #missing directions
    ExtrDlist<-NULL #extra directions
    
    BiFP<-NULL
    BiFN<-NULL
    comedges<-0
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
          } else comedges<-1
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
          } else if(adj[i,j]==1) comedges<-1
        }
      }
    }
    
    
    jointarcs<-adjacency2edgel(jointmat,nodes=nodelabs)
    graph.par(list(nodes=list(lty="solid", lwd=1, ...), graph=list(...)))
    graph.obj = new("graphNEL", nodes = nodelabs, edgeL = jointarcs,
                    edgemode = 'directed')
    
    if(b!=0) {  
      sg1 = subGraph(dynamicnames, graph.obj)
      sg2 = subGraph(staticnames, graph.obj)
      sgL = list(list(graph=sg1, cluster = TRUE),
                 list(graph=sg2, cluster = TRUE))
      colvector<-c(rep(statcol,b),
                   rep(dyn1col,dyn))
      names(colvector)<-nodelabs
    } else {
      sg1 = subGraph(dynamicnames, graph.obj)
      sgL = list(list(graph=sg1, cluster = TRUE))
      colvector<-c(rep(dyn1col,dyn))
      names(colvector)<-nodelabs
    }
    
    if(showcl) {
      graph.plot = Rgraphviz::layoutGraph(graph.obj, subGList = sgL) 
      } else {
        graph.plot = Rgraphviz::layoutGraph(graph.obj)
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
  
    tpname<-"true positive"
    fpname<-"false positive"
    fnname<-"false negative"
    edname<-"difference in direction"
    tpcol<-"black"
    edcol<-"blue"
    fpcol<-"red"
    fncol<-"grey"
    
    tplty<-1
    fplty<-1
    fnlty<-2
    edlty<-1
    
    if(is.null(FNlist) & is.null(BiFN)) {
      fncol<-NULL
      fnlwd<-NULL
      fnlty<-NULL
      fnname<-NULL
    }
    
    if(is.null(FPlist) & is.null(BiFP)) {
      fpcol<-NULL
      fplwd<-NULL
      fplty<-NULL
      fpname<-NULL
    }
    
    if(is.null(MissDlist) & is.null(ExtrDlist) & is.null(EDlist)) {
      edcol<-NULL
      edlwd<-NULL
      edlty<-NULL
      edname<-NULL
    }
    
    if(comedges==0) {
      tpcol<-NULL
      tplwd<-NULL
      tplty<-NULL
      tpname<-NULL
    }
    
    graph::nodeRenderInfo(graph.plot)<-list(fill=colvector,shape="circle",...)
    
    layout(matrix(c(1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,
                    2,2,2,3,3,3,4,4,4,
                    2,2,2,3,3,3,4,4,4), nrow = 9, ncol = 9, byrow = TRUE))
    
    if(b>0) {
      #graph
      Rgraphviz::renderGraph(graph.plot)
      #graph legend
      graph.par(list(graph=list(main="nodes(t):",cex.main=1.8),...))
      legendG <- Rgraphviz::layoutGraph(legendG,attrs = list(graph = list(rankdir = "TB")))
      graph::nodeRenderInfo(legendG)[["shape"]][c("stat","1")] <- shape
      graph::nodeRenderInfo(legendG)[["fill"]]["stat"] = statcol
      graph::nodeRenderInfo(legendG)[["fill"]]["1"] = dyn1col
      graph::edgeRenderInfo(legendG)[["lwd"]]["stat~1"] = 0
     
      Rgraphviz::renderGraph(legendG)

    } else {
      Rgraphviz::renderGraph(graph.plot)
    }
    par(mar = c(0,0,0,0))
    plot(1:10,1:10,type="n", axes = FALSE, xlab = "", ylab = "")
    plot(1:10,1:10,type="n", axes = FALSE, xlab = "", ylab = "")
    op <- par(cex = 1.7)
    legend("topleft",legend=c(tpname,
                             fpname,
                             fnname,
                             edname),lty=c(tplty,fplty,fnlty,edlty),
           col=c(tpcol,fpcol,fncol,edcol),bty="n",title="edges:",cex=0.7) 
  }
  
  if(struct=="trans") {
    
    n<-b+2*dyn
    
    if(b>0) {
      legadj<-matrix(0,nrow=3,ncol=3)
      legadj[1,2]<-1
      legadj[2,3]<-1
      colnames(legadj)<-c("stat","i","i+1")
      legendG<-m2graph(legadj)
      legcol<-c(statcol,dyn1col,dyn2col)
      colvector = c(rep(statcol,b),rep(dyn1col,dyn),rep(dyn2col,dyn))
      names(colvector)<-nodelabs
      names(colvector)<-nodelabs
    } else {
      legadj<-matrix(0,nrow=2,ncol=2)
      legadj[1,2]<-1
      colnames(legadj)<-c("i","i+1")
      legendG<-m2graph(legadj)
      legcol<-c(dyn1col,dyn2col)
      colvector = c(rep(dyn1col,dyn),rep(dyn2col,dyn))
      names(colvector)<-nodelabs
    }
    
    adj<-DBNcut(adj[1:(b+2*dyn),1:(b+2*dyn)],dyn,b)
    adjt<-DBNcut(adjt[1:(b+2*dyn),1:(b+2*dyn)],dyn,b)
    jointmat<-1*(adj|adjt)
    
    arcslist<-adjacency2edgel(jointmat,nodes=nodelabs)
    graph.obj = new("graphNEL", nodes = nodelabs, edgeL = arcslist,
                    edgemode = 'directed')
    
    staticnames<-nodelabs[1:b]
    dyn1names<-nodelabs[1:dyn+b]
    dyn2names<-nodelabs[1:dyn+b+dyn]
    
    sgStat = subGraph(staticnames, graph.obj)
    sgDyn1 = subGraph(dyn1names, graph.obj)
    sgDyn2 = subGraph(dyn2names, graph.obj)
    
    sgL = list(list(graph=sgStat, cluster = TRUE, attrs = c(rankdir="LR",rank="sink")),
               list(graph=sgDyn1, cluster = TRUE, attrs = c(rank="same")),
               list(graph=sgDyn2, cluster = TRUE, attrs = c(rank="same")))
    
    
    if(showcl) {
      graph.plot = Rgraphviz::layoutGraph(graph.obj, subGList = sgL,attrs = list(graph = list(rankdir = orientation))) 
    } else {
      graph.plot = Rgraphviz::layoutGraph(graph.obj,attrs = list(graph = list(rankdir = orientation)))
    }    
    
    graph::nodeRenderInfo(graph.plot)<-list(fill=colvector,shape="circle",...)
    
    
    #define edges with wrong directions
    FPlist<-NULL
    FNlist<-NULL
    EDlist<-NULL
    MissDlist<-NULL #missing directions
    ExtrDlist<-NULL #extra directions
    
    BiFP<-NULL
    BiFN<-NULL
    
    comedges<-0
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
          } else if(adj[i,j]==1) comedges<-1
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
          } else if(adj[i,j]==1) comedges<-1
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
    
    tpname<-"true positive"
    fpname<-"false positive"
    fnname<-"false negative"
    edname<-"difference in direction"
    tpcol<-"black"
    edcol<-"blue"
    fpcol<-"red"
    fncol<-"grey"
    
    tplty<-1
    fplty<-1
    fnlty<-2
    edlty<-1
    
    if(is.null(FNlist) & is.null(BiFN)) {
      fncol<-NULL
      fnlwd<-NULL
      fnlty<-NULL
      fnname<-NULL
    }
    
    if(is.null(FPlist) & is.null(BiFP)) {
      fpcol<-NULL
      fplwd<-NULL
      fplty<-NULL
      fpname<-NULL
    }
    
    if(is.null(MissDlist) & is.null(ExtrDlist) & is.null(EDlist)) {
      edcol<-NULL
      edlwd<-NULL
      edlty<-NULL
      edname<-NULL
    }
    
    if(comedges==0) {
      tpcol<-NULL
      tplwd<-NULL
      tplty<-NULL
      tpname<-NULL
    }
    
    graph.par(list(graph=list(...,cex.main=1.5),
                   nodes=list(lty="solid", lwd=1, fixedsize=FALSE,...)))

    
    layout(matrix(c(1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,
                    2,2,2,3,3,4,4,4,4,
                    2,2,2,3,3,4,4,4,4), nrow = 9, ncol = 9, byrow = TRUE))
    

   #plot graph
   Rgraphviz::renderGraph(graph.plot)
   #plot legend
   legendG <- Rgraphviz::layoutGraph(legendG,attrs = list(graph = list(rankdir = "TB")))
   if(b>0) {
     graph::nodeRenderInfo(legendG)[["shape"]][c("stat","i","i+1")] <- shape
     graph::nodeRenderInfo(legendG)[["fill"]]["stat"] = statcol
     graph::nodeRenderInfo(legendG)[["fill"]]["i"] = dyn1col
     graph::nodeRenderInfo(legendG)[["fill"]]["i+1"] = dyn2col
     graph::edgeRenderInfo(legendG)[["lwd"]]["stat~i"] = 0
     graph::edgeRenderInfo(legendG)[["lwd"]]["i~i+1"] = 0
   } else {
     graph::nodeRenderInfo(legendG)[["shape"]][c("i","i+1")] <- shape
     graph::nodeRenderInfo(legendG)[["fill"]]["i"] = dyn1col
     graph::nodeRenderInfo(legendG)[["fill"]]["i+1"] = dyn2col
     graph::edgeRenderInfo(legendG)[["lwd"]]["i~i+1"] = 0
   }
   
   graph.par(list(graph=list(main="nodes(t):",cex.main=1.8),...))
   Rgraphviz::renderGraph(legendG)
    
    par(mar = c(0,0,0,0))
    plot(1:10,1:10,type="n", axes = FALSE, xlab = "", ylab = "")
    plot(1:10,1:10,type="n", axes = FALSE, xlab = "", ylab = "")
    op <- par(cex = 1.7)
    legend(1,10,legend=c(tpname,
                             fpname,
                             fnname,
                             edname),lty=c(tplty,fplty,fnlty,edlty),
           col=c(tpcol,fpcol,fncol,edcol),bty="n",title="edges:",cex=0.7) 
  }
}

#' Plotting difference between two graphs
#' 
#' This function plots edges from two graphs in one and indicates similarities and differences between these graphs.
#' It is also possible to use this function for plotting mistakes in estimated graph when the ground truth graph is known.
#' 
#'@param graph1 object of class graphNEL or its adjacency matrix
#'@param graph2 object of class graphNEL or its adjacency matrix
#'@param estimated logical, indicates if graph1 is estimated graph and graph2 is ground truth DAG, TRUE by default; this affects the legend and colouring of the edges
#'@param name1 character, custom name for 'graph1'
#'@param name2 character, custom name for 'graph2'
#'@param clusters (optional) a list of nodes to be represented on the graph as clusters 
#'@param ... optional parameters passed to \code{Rgraphviz} plotting functions e.g. \code{main}, \code{fontsize}
#'@return plots the graph which includes edges from graph1 and graph2; edges which are different in graph1 compared to graph2 are coloured according to the type of a difference
#'@examples
#'Asiascore<-scoreparameters("bde",Asia)
#'Asiamap<-orderMCMC(Asiascore)
#'plotdiffs(Asiamap$DAG,Asiamat)
#'Asiacp<-pcalg::dag2cpdag(m2graph(Asiamat))
#'mapcp<-pcalg::dag2cpdag(m2graph(Asiamap$DAG))
#'plotdiffs(mapcp,Asiacp)
#'@author Polina Suter
#'@export
plotdiffs<-function(graph1,graph2,estimated=TRUE,name1="graph1",
                    name2="graph2",clusters=NULL, ...) {
  
  old.par<-par(no.readonly = TRUE)
  on.exit(par(old.par))
  
  if(is(graph1,"graphNEL")) {
    adj<-graph2m(graph1)
  } else {
    adj<-as.matrix(graph1)
  }
  if(is(graph2,"graphNEL")) {
    adjt<-graph2m(graph2)
  } else {
    adjt<-as.matrix(graph2)
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
  
  graph.par(list(nodes=list(lty="solid", lwd=1, ...),
                 graph=list(...)))
  comedges<-0
  if(estimated) {

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
        } else comedges<-1
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
        } else if(adj[i,j]==1) comedges<-1
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
  
  tpname<-"common edge"
  fpname<-"false positive"
  fnname<-"false negative"
  edname<-"difference in direction"
  tpcol<-"black"
  edcol<-"blue"
  fpcol<-"red"
  fncol<-"grey"
  
  tplty<-1
  fplty<-1
  fnlty<-2
  edlty<-1
  
  if(is.null(FNlist) & is.null(BiFN)) {
    fncol<-NULL
    fnlwd<-NULL
    fnlty<-NULL
    fnname<-NULL
  }
  
  if(is.null(FPlist) & is.null(BiFP)) {
    fpcol<-NULL
    fplwd<-NULL
    fplty<-NULL
    fpname<-NULL
  }
  
  if(is.null(MissDlist) & is.null(ExtrDlist) & is.null(EDlist)) {
    edcol<-NULL
    edlwd<-NULL
    edlty<-NULL
    edname<-NULL
  }
  
  if(comedges==0) {
    tpcol<-NULL
    tplwd<-NULL
    tplty<-NULL
    tpname<-NULL
  }
  graph::edgeRenderInfo(graph.plot)[["lwd"]] = 2
  #graph::edgeRenderInfo(graph.plot)[["col"]] = "black"
  #graph::edgeRenderInfo(graph.plot)[["lty"]] = 1
  
  layout(matrix(c(1,1,2), nrow = 1, ncol = 3, byrow = TRUE))
  Rgraphviz::renderGraph(graph.plot)

  par(mar = c(0,0,0,0))
  plot(1:10,1:10,type="n", axes = FALSE, xlab = "", ylab = "")
  op <- par(cex = 1.7)
  legend("center",legend=c(tpname,
                           fpname,
                           fnname,
                           edname),lty=c(tplty,fplty,fnlty,edlty),
         col=c(tpcol,fpcol,fncol,edcol),bty="n",cex=0.6) 
  } else {
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
          } else {
            if(adj[i,j]==1) comedges<-1
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
          } else {
            if(adj[i,j]==1)  comedges<-1
          }
        }
      }
    }
    
    fpcol<-"#74c476"
    fncol="#df65b0"
    
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
    graph::edgeRenderInfo(graph.plot)[["lwd"]] = 2
    
    if(!is.null(FPlist)) {
      FP<-apply(FPlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][FP] = fpcol
      graph::edgeRenderInfo(graph.plot)[["lty"]][FP] = "solid"
      graph::edgeRenderInfo(graph.plot)[["lwd"]][FP] = 1
    }
    if(!is.null(FNlist)) {
      FN<-apply(FNlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][FN] = fncol
      graph::edgeRenderInfo(graph.plot)[["lty"]][FN] = "solid"
      graph::edgeRenderInfo(graph.plot)[["lwd"]][FN] = 1
    } 
    if(!is.null(EDlist)) {
      ED<-apply(EDlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][ED] = "blue"
    }
    if(!is.null(BiFP)) {
      BiFP<-apply(BiFP, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][BiFP] = fpcol
      graph::edgeRenderInfo(graph.plot)[["lty"]][BiFP] = "solid"
      graph::edgeRenderInfo(graph.plot)[["lwd"]][BiFP] = 1
    }
    if(!is.null(BiFN)) {
      BiFN<-apply(BiFN, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][BiFN] = fncol
      graph::edgeRenderInfo(graph.plot)[["lty"]][BiFN] = "solid"
      graph::edgeRenderInfo(graph.plot)[["lwd"]][BiFN] = 1
      
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
    
    tpname<-"common edge"
    fpname<-paste("present only in",name1)
    fnname<-paste("present only in",name2)
    edname<-"difference in direction"
    tpcol<-"black"
    edcol<-"blue"
    
    tplwd<-2
    fplwd<-1
    fnlwd<-1
    edlwd<-2
    
    if(is.null(FNlist) & is.null(BiFN)) {
      fncol<-NULL
      fnlwd<-NULL
      fnlty<-NULL
      fnname<-NULL
    }
    
    if(is.null(FPlist) & is.null(BiFP)) {
      fpcol<-NULL
      fplwd<-NULL
      fplty<-NULL
      fpname<-NULL
    }
    
    if(is.null(MissDlist) & is.null(ExtrDlist)) {
      edcol<-NULL
      edlwd<-NULL
      edlty<-NULL
      edname<-NULL
    }
    
    if(comedges==0) {
      tpcol<-NULL
      tplwd<-NULL
      tplty<-NULL
      tpname<-NULL
    }
    
  
    layout(matrix(c(1,1,1,1,
                    1,1,1,1,
                    1,1,1,1,
                    1,1,1,1,
                    1,1,1,1,
                    2,2,3,3), nrow = 6, ncol = 4, byrow = TRUE))
    Rgraphviz::renderGraph(graph.plot)
    
    par(mar = c(0,0,0,0))
    plot(1:10,1:10,type="n", axes = FALSE, xlab = "", ylab = "")
    op <- par(cex = 1.7)
    legend("top",legend=c(tpname,fpname,fnname,edname),lwd=c(tplwd,fplwd,fnlwd,edlwd),
           col=c(tpcol,fpcol,fncol,edcol),bty="n",cex=0.6) 
  }
}

#' Highlighting similarities between two graphs
#' 
#' This function plots nodes and edges from two graphs in one and indicates similarities between these graphs.
#' 
#'@param graph1 binary adjacency matrix of a graph
#'@param graph2 binary adjacency matrix of a graph, column names should coincide with column names of 'graph1'
#'@param name1 character, custom name for 'graph1'; when NULL no legend will be plotted
#'@param name2 character, custom name for 'graph2'
#'@param bidir logical, defines if arrows of bidirected edges are drawn; FALSE by defauls.
#'@param ... optional parameters passed to \pkg{Rgraphviz} plotting functions e.g. \code{main}, \code{fontsize}
#'@return plots the graph which includes nodes and edges two graphs; nodes which are connected to at least one other node in both graphs are plotted only once and coloured orange, edges which are shared by two graphs
#'are coloured orange; all other nodes and edges a plotted once for each 'graph1' and 'graph2' and coloured blue and green accordingly.
#'@author Polina Suter
#'@export
plot2in1<-function(graph1, graph2, name1=NULL,
                   name2=NULL,bidir=FALSE, ...) {
  
  if(is(graph1,"graphNEL")) {
    graph1<-graph2m(graph1)
  } else if (!is.matrix(graph1)){
    graph1<-as.matrix(graph1)
  }
  
  if(is(graph2,"graphNEL")) {
    graph2<-graph2m(graph2)
  } else if (!is.matrix(graph2)){
    graph2<-as.matrix(graph2)
  }
  
  if(!all(colnames(graph1)==colnames(graph2))) stop("adjacency matrices 'graph1' and 'graph2' have different column names!")
  
  old.par<-par(no.readonly = TRUE)
  on.exit(par(old.par))
  
  mycolors <- c("#cbd5e8", "#ccebc5","#fdcdac","#E8CED4FF")
  glist<-list()
  glist[[1]]<- connectedSubGraph(graph1)
  glist[[2]]<- connectedSubGraph(graph2)
  
  clustm<-list()
  
  
  clustm[[1]]<-setdiff(colnames(glist[[1]]),colnames(glist[[2]]))
  clustm[[2]]<-setdiff(colnames(glist[[2]]),colnames(glist[[1]]))
  clustm[[3]]<-intersect(colnames(glist[[1]]),colnames(glist[[2]]))
  numsubg<-length(clustm)
  
  allmuts<-unique(c(colnames(glist[[1]]),colnames(glist[[2]])))
  
  glist[[1]]<-getSubGraph(graph1,allmuts)
  glist[[2]]<-getSubGraph(graph2,allmuts)

  nodelabs<-colnames(glist[[2]])
  n<-nrow(glist[[2]])
  jointgraph<-1*Reduce("|",glist)
  jointarcs<-adjacency2edgel(jointgraph,nodes=nodelabs)
  graph.obj = new("graphNEL", nodes = nodelabs, edgeL = jointarcs,
                  edgemode = 'directed')

  sg<-list()
  subGList<-list()
  for(i in 1:numsubg) {
    if(!is.null(clustm[[i]])) {
      sg[[i]] <-subGraph(clustm[[i]],  graph.obj)
      subGList[[i]]<-list(graph = sg[[i]], cluster = TRUE)
    }
  }
  
  graph.par(list(nodes=list(lty="solid", lwd=1, ...), graph=list(...)))
  
  graph.plot = Rgraphviz::layoutGraph(graph.obj,subGList = subGList)
  graph::nodeRenderInfo(graph.plot)[["fill"]][clustm[[1]]] = mycolors[1] #M
  graph::nodeRenderInfo(graph.plot)[["fill"]][clustm[[2]]] = mycolors[2] #T
  graph::nodeRenderInfo(graph.plot)[["fill"]][clustm[[3]]] = mycolors[3] #CNA

  sumgraph<-Reduce("+",glist)
  arcl<-list()
  
  for(i in 1:2) {
    edgy<-which(glist[[i]]>0,arr.ind = TRUE)
    arcl[[i]]<-matrix(ncol=2,nrow=nrow(edgy))
    for(j in 1:nrow(edgy)) {
      arcl[[i]][j,]<-c(nodelabs[edgy[j,1]],nodelabs[edgy[j,2]])
    }
  }
  
  commedges<-which(sumgraph>1,arr.ind = TRUE)
  if(nrow(commedges>0)) {
    arcl[[3]]<-matrix(ncol=2,nrow=nrow(commedges))
    for(j in 1:nrow(commedges)) {
      arcl[[3]][j,]<-c(nodelabs[commedges[j,1]], nodelabs[commedges[j,2]])
    }
  }
  
  for(i in 1:length(arcl)) {
    if(!is.null(arcl[[i]])) {
      arcl[[i]]<-apply(arcl[[i]], 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][arcl[[i]]] = mycolors[i]
    }
  }
  
  
  u <- names(which(graph::edgeRenderInfo(graph.plot)[["direction"]] == "both"))
  if(!bidir) {
    graph::edgeRenderInfo(graph.plot)[["arrowhead"]][u] = "none"
    graph::edgeRenderInfo(graph.plot)[["arrowtail"]][u] = "none"
  }
  graph::edgeRenderInfo(graph.plot)[["lwd"]]<-2
  
  layout(matrix(c(1,1,1,1,1,
                  1,1,1,1,1,
                  1,1,1,1,1,
                  1,1,1,1,1,
                  1,1,1,1,1,
                  1,1,1,2,1), nrow = 6, ncol = 5, byrow = TRUE))
  Rgraphviz::renderGraph(graph.plot)
  if(!is.null(name1)){
      par(mar = c(0,0,0,0))
      plot(1:10,1:10,type="n", axes = FALSE, xlab = "", ylab = "")
      op <- par(cex = 1.7)
      legend("topright",
             col=c("#cbd5e8", "#ccebc5","#fdcdac"),
             lwd = 2,
             lty = 1,
             legend = c(name1,name2,"both"),
             cex = 0.6,
             bg = NA,bty = "n")
      legend("topright",
             col = "black",
             pt.bg = c("#cbd5e8", "#ccebc5","#fdcdac"),
             pch = 21,
             lwd = 1,
             legend = c(name1,name2,"both"),
             cex = 0.6,
             lty = 0,
             bty = "n")
  }
  
}

