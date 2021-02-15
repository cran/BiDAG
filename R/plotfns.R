#'Plotting posterior probabilities of single edges
#'
#'This function plots posterior probabilities of all possible edges in the graph as a function of MCMC iterations. It can be used for convergence diagnostics of MCMC
#'sampling algorithms order MCMC and partition MCMC.
#' @param MCMCtrace an object of class MCMCres
#' @param cutoff number representing a threshold of posterior probability below which lines will not be plotted
#' @param pdag logical, when true DAGs in a sample will be first coverted to CPDAGs
#' @param onlyedges (optional) binary matrix, only edges corresponding to entries which equal 1 will be plotted
#' @param highlight (optional) binary matrix, edges corresponding to entries which equal 1 are highlighted with "red"
#'@return plots posterior probabilities of edges in the graph as a function of MCMC iterations
#'@examples
#'score100<-scoreparameters("bde", Asia[1:100,])
#'orderfit100<-orderMCMC(score100,plus1=TRUE)
#'\dontrun{
#'score5000<-scoreparameters("bde", Asia)
#'orderfit5000<-orderMCMC(score5000,plus1=TRUE)
#'plotpedges(orderfit100, pdag=TRUE)
#'plotpedges(orderfit5000, pdag=TRUE)
#'}
#'@author Polina Suter
#'@export
plotpedges<-function(MCMCtrace,cutoff=0.2,pdag=FALSE,onlyedges=NULL,highlight=NULL) {
  MCMCtrace<-MCMCtrace$traceadd$incidence
  if(is.null(MCMCtrace)) {
    stop("no saved MCMC steps found! try chainout=TRUE when sampling")
  }
  #cols5<-c("#f2f0f7","#b3cde3","#8c96c6","#8856a7","#810f7c","#41ae76")
  cols5<-c("#8c96c6","#8c96c6","#8c96c6","#8c96c6","#8c96c6","#41ae76")
  
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
    cutoffelems<-which(onlyedges==1)
  } else {
    cutoffelems<-which(posteriors[[lchain]]>cutoff)
  }

  numelem<-length(cutoffelems)
  if(is.null(highlight)) {
    redi<-c()
  } else {
    redi<-which(highlight==1) 
  }

  postvec<-list()
  colvec<-vector()
  k<-1
  for(i in cutoffelems){
    postvec[[k]]<-lapply(posteriors,function(x)x[i])
    if(i%in%redi) {
      colvec[k]<-6
    } else {
      colvec[k]<-defcolrange(posteriors[[lchain]][i])
    }
    k<-k+1
  }
  colvec<-as.integer(colvec)
  plot(x=c(1:lchain),y=postvec[[1]],type="l",
       col=cols5[colvec[1]],xlab="MCMC step",ylab="estimated posterior",
       main="posterior probabilities of single edges",ylim=c(0,1),cex.main=1)
  for(i in 2:numelem){
    lines(x=c(1:lchain),postvec[[i]],type="l",col=cols5[colvec[i]])
  }
}


#'Comparing posterior probabilitites of single edges based on several samples
#'
#'This function can be used to compare posterior probabilities of edges in a graph 
#'based on two samples of graphs
#'@param pmat a list of square matrices, representing posterior probabilities of single edges in a Bayesian network
#'@param highlight numeric, defines maximum acceptable difference between posterior probabilities of an edge in two samples; points corresponding to higher differences are highlighted 
#'@param cut numeric value corresponding to a minimum posterior probabilitity which is included into calculation of squared correlation and MSE
#'@param main character string, a title for the plot; ignored if length of pmat higher than 2
#'@param xlab character string, a title for the x-axis; ignored if length of pmat higher than 2
#'@param ylab character string, a title for the y-axis; ignored if length of pmat higher than 2
#'@return plots posterior probabilitites of single edges from two MCMC runs; returns squared correlation and MSE of posterior probabilities higher than the value defined by the argument cut
#'@examples
#'Asiascore<-scoreparameters("bde", Asia)
#'\dontrun{
#'orderfit<-list()
#'orderfit[[1]]<-orderMCMC(Asiascore,MAP=FALSE)
#'orderfit[[2]]<-orderMCMC(Asiascore,MAP=FALSE)
#'pedges<-lapply(orderfit,edgep,pdag=TRUE)
#'plotpcor(pedges)
#'}
#'@author Polina Suter
#'@export
plotpcor<-function(pmat,highlight=0.3,cut=0.05,main="",
                   xlab="sample 1",
                   ylab="sample 2") {
  oldpar<-par(no.readonly = TRUE)
  nruns<-length(pmat)
  if(nruns==1) stop("the number of matrices in the list must be at least two!")
  if(nruns>2) {
  
    vecy<-list()
    for(i in 1:nruns) {
      vecy[[i]]<-as.vector(pmat[[i]])
    }

    par(mfrow=c(nruns,nruns)) 
    
    for(i in 1:nruns) {
      for(j in 1:nruns) {
        pcorcore(vecy[[i]],vecy[[j]],pmat[[i]],pmat[[j]],highlight,i,j)
      }
    }
    par(oldpar)
  } else {
  varnames<-colnames(pmat[[1]])
  vec1<-as.vector(pmat[[1]])
  vec2<-as.vector(pmat[[2]])
  
  if(highlight<1){
    diffmat<-abs(pmat[[1]]-pmat[[2]])
    pointstohighlight<-which(diffmat>highlight)
    diffedges<-which(diffmat>highlight,arr.ind=TRUE)
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
  if(length(diffedges)>0) {
    res$diffedges<-matrix(varnames[diffedges],ncol=2,nrow=nrow(diffedges))
    colnames(res$diffedges)<-c("from","to")
  }
  return(res)
  }


}





pcorcore<-function(vec1,vec2,mat1,mat2,highlight,i,j) {
  if(highlight<1){
    diffmat<-abs(mat1-mat2)
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
       xlab=paste("sample",i),ylab=paste("sample",j),
       xlim=c(0,1),ylim=c(0,1))
  lines(vec1,vec2,type="p",col="grey")
  if(highlight<1) lines(vec1high,vec2high,type="p",col="red")
}


