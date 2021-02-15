#' Plotting a DBN
#' 
#' This function can be used for plotting initial and transition structures of a dynamic Bayesian network.
#' 
#'
#'@param DBN binary matrix (or a graph object) representing a 2-step DBN (compact or unrolled)
#'@param struct option used to determine if the initial or the transition structure should be plotted; accaptable values are init or trans
#'@param n.dynamic number of dynamic variables in one time slice of a DBN
#'@param n.static number of static variables in one time slice of a DBN; note that for function to work correctly all static variables have to be in the first n.static columns of the matrix
#'@examples
#'plotDBN(DBNmat, "init", n.dynamic=12,n.static=3)
#'plotDBN(DBNmat, "trans", n.dynamic=12,n.static=3)
#'
#' @author Polina Suter
#' @export
plotDBN<-function(DBN,struct=c("init","trans"),n.dynamic,n.static){
  
  old.par <- par(no.readonly = TRUE)
  
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
      colnames(legadj)<-c("stat","1")
      legadj[1,2]<-1
      legendG<-m2graph(legadj)
      staticnames<-nodelabs[1:n.static]
      legendG<-Rgraphviz::layoutGraph(legendG)
    }
    
    dynamicnames<-nodelabs[1:n.dynamic+n.static]
    
    adj<-DBN[1:(n.dynamic+n.static),1:(n.dynamic+n.static)]
    arcslist<-adjacency2edgel(adj,nodes=nodelabs)
    graph.obj = new("graphNEL", nodes = nodelabs, edgeL = arcslist,
                    edgemode = 'directed')
    subGList<-list()
    sg<-list()
    
    graph.par(list(nodes=list(col=dyn1col, lty="solid", lwd=2, fontsize=14,cex=1.1),graph=list(main="Initial structure",cex.main=1.5)))
    if(n.static!=0) {  
      sg1 = subGraph(dynamicnames, graph.obj)
      sg2 = subGraph(staticnames, graph.obj)
      sgL = list(list(graph=sg1, cluster = TRUE),
                 list(graph=sg2, cluster = TRUE))
      graph.obj <- Rgraphviz::layoutGraph(graph.obj, subGList= sgL)
      graph::nodeRenderInfo(graph.obj)[["fill"]][staticnames] = statcol
      graph::nodeRenderInfo(graph.obj)[["fill"]][dynamicnames] = dyn1col
      
    } else {
      graph.obj <- Rgraphviz::layoutGraph(graph.obj)
    }
    
    if(n.static>0) {
      layout(matrix(c(1,1,1,1,3,
                      1,1,1,1,2,
                      1,1,1,1,3), nrow = 3, ncol = 5, byrow = TRUE))
      
      graph::nodeRenderInfo(legendG)[["fill"]]["stat"] = statcol
      graph::edgeRenderInfo(legendG)[["lwd"]]["stat~1"] = 0
      Rgraphviz::renderGraph(graph.obj )
      graph.par(list(nodes=list(col=dyn1col, lty="solid", lwd=2, fontsize=14,cex=1.1),graph=list(main="nodes (t):")))
      Rgraphviz::renderGraph(legendG)
    } else {
      graph.par(list(nodes=list(col=dyn1col, lty="solid", lwd=2, fontsize=14,cex=1.1),graph=list(main="nodes (t):")))
      Rgraphviz::renderGraph(graph.obj )
    }
  }
  
  if(struct=="trans") {
    if(n.static>0) {
      legadj<-matrix(0,nrow=3,ncol=3)
      legadj[1,2]<-1
      legadj[2,3]<-1
      colnames(legadj)<-c("stat","i","i+1")
      legendG<-m2graph(legadj)
      
    } else {
      legadj<-matrix(0,nrow=2,ncol=2)
      legadj[1,2]<-1
      colnames(legadj)<-c("i","i+1")
      legendG<-m2graph(legadj)
    }
    legendG<-Rgraphviz::layoutGraph(legendG)
    
    adjt<-DBNcut(DBN[1:(n.static+2*n.dynamic),1:(n.static+2*n.dynamic)],n.dynamic,n.static)
    graph.obj<-m2graph(adjt)
    
    dyn1names<-nodelabs[1:n.dynamic+n.static]
    dyn2names<-nodelabs[1:n.dynamic+n.static+n.dynamic]
    sgDyn1 = subGraph(dyn1names, graph.obj)
    sgDyn2 = subGraph(dyn2names, graph.obj)
    
    if(n.static>0) {
      staticnames<-nodelabs[1:n.static]
      sgStat = subGraph(staticnames, graph.obj)
      sgL = list(list(graph=sgStat, cluster = TRUE, attrs = c(rankdir="TB",rank="sink")),
                 list(graph=sgDyn1, cluster = TRUE, attrs = c(rank="same")),
                 list(graph=sgDyn2, cluster = TRUE, attrs = c(rank="same")))
    } else {
      sgL = list(list(graph=sgDyn1, cluster = TRUE, attrs = c(rank="same")),
                 list(graph=sgDyn2, cluster = TRUE, attrs = c(rank="same")))
    }
    
    graph.par(list(nodes=list(col=dyn1col, lty="solid", lwd=2, fontsize=14,cex=1.1),graph=list(main="Transition structure")))
    graph.obj <- Rgraphviz::layoutGraph(graph.obj, subGList= sgL)
    
    if(n.static>0) graph::nodeRenderInfo(graph.obj)[["fill"]][staticnames] = statcol
    graph::nodeRenderInfo(graph.obj)[["fill"]][dyn1names] = dyn1col
    graph::nodeRenderInfo(graph.obj)[["fill"]][dyn2names] = dyn2col
    
    layout(matrix(c(1,1,1,1,1,3,3,
                    1,1,1,1,1,2,2,
                    1,1,1,1,1,3,3), nrow = 3, ncol = 7, byrow = TRUE))
    Rgraphviz::renderGraph(graph.obj)
    graph.par(list(nodes=list(col=dyn1col, lty="solid", lwd=1, fontsize=14,cex=1.1),graph=list(main="nodes (t):")))
    if(n.static>0) {
      graph::nodeRenderInfo(legendG)[["fill"]]["stat"] = statcol
      graph::nodeRenderInfo(legendG)[["fill"]]["i+1"] = dyn2col
      graph::edgeRenderInfo(legendG)[["lwd"]]["stat~i"] = 0
      graph::edgeRenderInfo(legendG)[["lwd"]]["i~i+1"] = 0
    } else {
      graph::nodeRenderInfo(legendG)[["fill"]]["i+1"] = dyn2col
      graph::edgeRenderInfo(legendG)[["lwd"]]["i~i+1"] = 0
    }
    Rgraphviz::renderGraph(legendG)
  }
  par(old.par)
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
#'@param n.dynamic number of dynamic variables in one time slice of a DBN
#'@param n.static number of static variables in one time slice of a DBN; note that for function to work correctly all static variables have to be in the first n.static columns of the matrix
#'@return plots the graph which includes edges from graph1 and graph2, however edges which are different in graph1 compared to graph2 are coloured according to the type of a mistake: false positive with red, false negative with dashed grey, error in direction with magenta
#'@examples
#'dbnscore<-scoreparameters("bge",DBNdata,
#'dbnpar = list(samestruct=TRUE, slices=5, stationary=TRUE),
#'DBN=TRUE,bgnodes=c(1,2,3))
#'\dontrun{
#'orderDBNfit<-iterativeMCMC(dbnscore,chainout = TRUE, mergetype = "skeleton",scoreout=TRUE,alpha=0.4)
#'plotdiffs.DBN(orderDBNfit$max$DAG,DBNmat,struct="trans",n.dynamic=12,n.static=3)
#'plotdiffs.DBN(orderDBNfit$max$DAG,DBNmat,struct="init",n.dynamic=12,n.static=3)
#'}
#'@export
#'@author Polina Suter
plotdiffs.DBN<-function(eDBN,trueDBN,struct=c("init","trans"),
                        n.dynamic,n.static=0) {
  old.par <- par(no.readonly = TRUE)
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
      colnames(legadj)<-c("stat","1")
      legadj[1,2]<-1
      legendG<-m2graph(legadj)
      staticnames<-nodelabs[1:n.static]
      legcol<-c(statcol,dyn1col)
      names(legcol)<-c("stat","1")
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
      sg1 = subGraph(dynamicnames, graph.obj)
      sgL = list(list(graph=sg1, cluster = TRUE))
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
    
    graph.par(list(graph=list(main="Comparison of DBN initial structures")))
    
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
    
    if(n.static>0) {
      layout(matrix(c(1,1,1,1,4,
                      1,1,1,1,2,
                      1,1,1,1,3), nrow = 3, ncol = 5, byrow = TRUE))
      Rgraphviz::renderGraph(graph.plot)
      graph::plot(legendG,attrs=list(graph=list(rankdir="TB"),
                                     edge=list(lwd=0)),
                  nodeAttrs=list(fillcolor=legcol),
                  main="nodes (t):",cex.main=1.5)
    } else {
      Rgraphviz::renderGraph(graph.plot)
    }
    par(mar = c(0,0,0,0))
    plot(1:10,1:10,type="n", axes = FALSE, xlab = "", ylab = "")
    op <- par(cex = 1.7)
    legend("center",legend=c(tpname,
                             fpname,
                             fnname,
                             edname),lty=c(tplty,fplty,fnlty,edlty),
           col=c(tpcol,fpcol,fncol,edcol),bty="n",title="edges:",cex=0.7) 
  }
  
  if(struct=="trans") {
    
    n<-n.static+2*n.dynamic
    
    if(n.static>0) {
      legadj<-matrix(0,nrow=3,ncol=3)
      legadj[1,2]<-1
      legadj[2,3]<-1
      colnames(legadj)<-c("stat","i","i+1")
      legendG<-m2graph(legadj)
      legcol<-c(statcol,dyn1col,dyn2col)
      names(legcol)<-c("stat","i","i+1")
      
      colvector = c(rep(statcol,n.static), rep(dyn1col,n.dynamic),rep(dyn2col,n.dynamic))
      names(colvector)<-nodelabs
    } else {
      legadj<-matrix(0,nrow=2,ncol=2)
      legadj[1,2]<-1
      colnames(legadj)<-c("i","i+1")
      legendG<-m2graph(legadj)
      legcol<-c(dyn1col,dyn2col)
      names(legcol)<-c("i","i+1")
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
    
    graph.par(list(graph=list(main="Comparison of DBN transition structures")))
    layout(matrix(c(1,1,1,1,1,4,4,
                    1,1,1,1,1,2,2,
                    1,1,1,1,1,3,3), nrow = 3, ncol = 7, byrow = TRUE))
    
    Rgraphviz::renderGraph(graph.plot)
    
    graph::plot(legendG,attrs=list(graph=list(rankdir="TB"),edge=list(lwd=0)),
                nodeAttrs=list(fillcolor=legcol),main="nodes (t):",cex.main=1.5)
    
    par(mar = c(0,0,0,0))
    plot(1:10,1:10,type="n", axes = FALSE, xlab = "", ylab = "")
    op <- par(cex = 1.7)
    legend("center",legend=c(tpname,
                             fpname,
                             fnname,
                             edname),lty=c(tplty,fplty,fnlty,edlty),
           col=c(tpcol,fpcol,fncol,edcol),bty="n",title="edges:",cex=0.7) 
  }
  par(old.par)
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
#'@return plots the graph which includes edges from graph1 and graph2, however edges which are different in graph1 compared to graph2 are coloured according to the type of a mistake: false positive with red, false negative with dashed grey, error in direction with magenta
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
                    name2="graph2",clusters=NULL) {
  old.par <- par(no.readonly = TRUE)
  if(!is.matrix(graph1)) {
    adj<-graph2m(graph1)
  } else {
    adj<-graph1
  }
  if(!is.matrix(graph2)) {
    adjt<-graph2m(graph2)
  } else {
    adjt<-graph2
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
  
  graph.par(list(nodes=list(lty="solid", lwd=1, fontsize=18),
                 graph=list(main=paste("similarities/differences between",name1,"and",name2), cex.main=1.5)))
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
    
  
    layout(matrix(c(1,1,2), nrow = 1, ncol = 3, byrow = TRUE))
    Rgraphviz::renderGraph(graph.plot)
    
    par(mar = c(0,0,0,0))
    plot(1:10,1:10,type="n", axes = FALSE, xlab = "", ylab = "")
    op <- par(cex = 1.7)
    legend("bottom",legend=c(tpname,fpname,fnname,edname),lwd=c(tplwd,fplwd,fnlwd,edlwd),
           col=c(tpcol,fpcol,fncol,edcol),bty="n",cex=0.6) 
  }
  par(old.par)
}









