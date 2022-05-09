#' Converting a single BiDAG chain to mcmc object 
#' 
#' This function converts a single object of one of the BiDAG classes,
#' namely 'orderMCMC' and 'partitionMCMC' to an object of class 'mcmc'. This object can
#' be further used for convergence and mixing diagnostics implemented in the package coda
#'  
#' @param MCMCtrace object of class \code{orderMCMC} or \code{partitionMCMC} 
#' @param edges logical, when FALSE (default), then only DAG score trace is extracted; when TRUE, a trace of posterior probabilities is extracted for every edge (based on the last 200 sampled structures) resulting in n*(n-1)) trace vectors, where n is the number of nodes in the network
#' @param pdag logical, when edges=TRUE, defines is the DAGs ar converted to CPDAGs prior to computing posterior probabilities; ignored otherwise
#' @param p numeric, between 0 and 1; defines the minimum probability for including posterior traces in the returned objects (for probabilities close to 0 PRSF diagnostics maybe too conservative)
#' @param burnin numeric between \code{0} and \code{1}, indicates the percentage of the samples which will be discarded as 'burn-in' of the MCMC chain; the rest  of the samples will be used to calculate the posterior probabilities; 0.2 by default
#' @param window integer, defines a number of DAG samples for averaging and computing edges' posterior probabilities; ignored when edges=FALSE 
#' @param cumulative logical, indicates is posterior probabilities should be calculated based on a cumulative sample of DAGs, where 25\% of first samples are discarded
#' @return Object of class \code{mcmc} from the package \pkg{coda}
#'@examples
#'\dontrun{
#'myscore<-scoreparameters("bde",Asia)
#'ordersample<-sampleBN(myscore,"order")
#'order_mcmc<-bidag2coda(ordersample)
#'par(mfrow=c(1,2))
#'densplot(order_mcmc)
#'traceplot(order_mcmc)
#'}
#'@author Polina Suter
#'@export
bidag2coda<-function(MCMCtrace, edges=FALSE, pdag=TRUE, p=0.1,burnin=0.2,window=100,cumulative=FALSE) {
  
  thin<-MCMCtrace$info$iterations/MCMCtrace$info$samplesteps-1
  lmcmc<-length(MCMCtrace$trace)
  firstsample<-ceiling(lmcmc*burnin)
  
  if(!edges) {
    mcmcobj<-mcmc(data=as.matrix(MCMCtrace$trace[firstsample:lmcmc],nrow=1), start = firstsample, thin = thin)
  } else {
    MCMCtrace<-MCMCtrace$traceadd$incidence
    if(is.null(MCMCtrace)) {
      stop("no saved MCMC steps found! try chainout=TRUE when sampling")
    }
    lchain<-length(MCMCtrace)
    if(pdag==TRUE) {
      MCMCtrace<-lapply(MCMCtrace,dagadj2cpadj)
    }
    countmatrix<-MCMCtrace[[1]]
    posteriors<-as.vector(MCMCtrace[[1]])
    curstartstep<-1
    for (i in 2:lchain) {
      if(cumulative) {
        newstartstep<-ceiling(0.25*i)
        countmatrix<-countmatrix+MCMCtrace[[i]]
        if(newstartstep>curstartstep) {
          countmatrix<-countmatrix-MCMCtrace[[curstartstep]]
          curstartstep<-newstartstep
        }
        posteriors<-cbind(posteriors,as.vector(as.matrix(countmatrix)/(i-curstartstep)))
      } else {
        if(i<(window+1)) {
          countmatrix<-countmatrix+MCMCtrace[[i]]
          posteriors<-cbind(posteriors,as.vector(as.matrix(countmatrix)/i))
        } else {
          countmatrix<-countmatrix+MCMCtrace[[i]]-MCMCtrace[[i-window]]
          posteriors<-cbind(posteriors,as.vector(as.matrix(countmatrix)/window))
        }
      }
    }
    if(p>0) {
      edg<-which(posteriors[,ncol(posteriors)>p])
      posteriors<-posteriors[edg,]
    }
    mcmcobj<-mcmc(data=t(posteriors[,firstsample:lmcmc]), start = firstsample, thin = thin)
  }
  return(mcmcobj)
}

#' Converting multiple BiDAG chains to mcmc.list
#' 
#' This function converts a list of objects of classes
#' 'orderMCMC' and 'partitionMCMC' to an object of class 'mcmc.list'. This object can
#' be further used for convergence and mixing diagnostics implemented in the R-package coda.
#'  
#' @param MCMClist a list of objects of classes \code{orderMCMC} or \code{partitionMCMC} 
#' @param edges logical, when FALSE (default), then only DAG score traces are extracted from each element of MCMClist and all other parameters are ignored; when TRUE, a trace of posterior probabilities is computed for every edge resulting in a maximum of n*n trace vectors, where n is the number of nodes in the network
#' @param pdag logical, when edges=TRUE, defines is the DAGs ar converted to CPDAGs prior to computing posterior probabilities; ignored otherwise
#' @param p numeric, between 0 and 1; defines the minimum probability for including posterior traces in the returned objects (for probabilities close to 0, PRSF diagnostics maybe too conservative; the threshold above 0 is recommended)
#' @param burnin numeric between \code{0} and \code{1}, indicates the percentage of the samples which will be discarded as 'burn-in' of the MCMC chain; the rest  of the samples will be used to calculate the posterior probabilities; 0.2 by default
#' @param window integer, defines a number of DAG samples for averaging and computing edges' posterior probabilities; ignored when edges=FALSE 
#' @param cumulative logical, indicates is posterior probabilities should be calculated based on a cumulative sample of DAGs, where 25\% of first samples are discarded
#' @return Object of class \code{mcmc.list} from the package \pkg{coda}
#'@examples
#'\dontrun{
#'scoreBoston<-scoreparameters("bge",Boston)
#'ordershort<-list()
#'#run very short chains -> convergence issues
#'ordershort[[1]] <- sampleBN(scoreBoston, algorithm = "order", iterations=2000)
#'ordershort[[2]] <- sampleBN(scoreBoston, algorithm = "order", iterations=2000)
#'codashort_edges<-bidag2codalist(ordershort,edges=TRUE,pdag=TRUE,p=0.05,burnin=0.2,window=10)
#'gd_short<-gelman.diag(codashort_edges, transform=FALSE, autoburnin=FALSE, multivariate=FALSE)
#'length(which(gd_short$psrf[,1]>1.1))/(length(gd_short$psrf[,1]))
#'}
#'@author Polina Suter
#'@references Robert J. B. Goudie and Sach Mukherjee (2016). A Gibbs Sampler for Learning DAGs. J Mach Learn Res. 2016 Apr; 17(30): 1–39.
#'@export
bidag2codalist<-function(MCMClist, edges=FALSE, pdag=TRUE,p=0.1,burnin=0.2,window=10,cumulative=FALSE) {
  thin<-vector()
  lmcmc<-length(MCMClist[[1]]$trace)
  firstsample<-ceiling(lmcmc*burnin)
  for(i in 1:length(MCMClist)) {
    thin[i]<-MCMClist[[i]]$info$iterations/MCMClist[[i]]$info$samplesteps-1
  }
  if(edges==FALSE | p==0) {
    for(i in 1:length(MCMClist)) {
      MCMClist[[i]]<-bidag2coda(MCMClist[[i]],edges=edges,pdag=pdag,p=0)
    }
  } else {
    for(i in 1:length(MCMClist)) {
      MCMClist[[i]]<-bidag2codacore(MCMClist[[i]],edges=edges,pdag=pdag,burnin=burnin,window=window,cumulative=cumulative)
    }
    edg<-c()
    for(i in 1:length(MCMClist)) {
      edg<-union(edg,which(MCMClist[[i]][,ncol(MCMClist[[i]])]>p))
    }
    for(i in 1:length(MCMClist)) {
      MCMClist[[i]]<- MCMClist[[i]][edg,]
      MCMClist[[i]]<-mcmc(data=t(MCMClist[[i]][,firstsample:lmcmc]), start = firstsample, thin = thin[i])
    }
  }
  return(mcmc.list(MCMClist))
  
}


bidag2codacore<-function(MCMCtrace, edges=FALSE, pdag=TRUE,burnin=0.2,window=100,cumulative=FALSE) {

    MCMCtrace<-MCMCtrace$traceadd$incidence
    if(is.null(MCMCtrace)) {
      stop("no saved MCMC steps found! try chainout=TRUE when sampling")
    }
    lchain<-length(MCMCtrace)
    if(pdag==TRUE) {
      MCMCtrace<-lapply(MCMCtrace,dagadj2cpadj)
    }
    countmatrix<-MCMCtrace[[1]]
    posteriors<-as.vector(MCMCtrace[[1]])
    curstartstep<-1
    for (i in 2:lchain) {
      if(cumulative) {
        newstartstep<-ceiling(0.25*i)
        countmatrix<-countmatrix+MCMCtrace[[i]]
        if(newstartstep>curstartstep) {
          countmatrix<-countmatrix-MCMCtrace[[curstartstep]]
          curstartstep<-newstartstep
        }
        posteriors<-cbind(posteriors,as.vector(as.matrix(countmatrix)/(i-curstartstep)))
      } else {
        if(i<(window+1)) {
          countmatrix<-countmatrix+MCMCtrace[[i]]
          posteriors<-cbind(posteriors,as.vector(as.matrix(countmatrix)/i))
        } else {
          countmatrix<-countmatrix+MCMCtrace[[i]]-MCMCtrace[[i-window]]
          posteriors<-cbind(posteriors,as.vector(as.matrix(countmatrix)/window))
        }
      }
    }

  return(posteriors)
}