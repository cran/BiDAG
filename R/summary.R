#SUMMARY methods for classes:
#scoreparameters
#scorespace
#orderMCMC
#partitionMCMC
#iterativeMCMC
#itercomp
#samplecomp

#' Summary of object of class 'scoreparameters'
#'
#' @param object object of class 'scoreparameters'
#' @param ... ignored
#' 
#' @rdname scoreparameters
#' @method summary scoreparameters
#' @export
summary.scoreparameters <-function(object, ...){
  cat("object of class 'scoreparameters' \n")
  cat("number of nodes (variables):", object$n, "\n")
  cat("number of observations:", nrow(object$data), "\n")
  cat("score type:", object$type, "\n")
  
  if(object$DBN==TRUE) {
    cat("score object created for a DBN \n")
    if(object$bgn>0) {
      cat("static nodes (no parents are allowed):", object$static, "\n")
    }
  } else {
    if(object$bgn>0) {
      cat("background nodes (no parents are allowed):", object$bgnodes, "\n")
    }
  }
  
  if(!is.null(object$weightvector)) {
    cat("data is weighted \n")
  }
}




#' Summary of object of class 'scorespace'
#'
#' @param object object of class 'scorespace'
#' @param ... ignored
#' 
#' @rdname scorespace
#' @method summary scorespace
#' @export
summary.scorespace <-function(object, ...){
  cat("object of class 'scorespace'")
  cat("\n\n")
  n<-ncol(object$adjacency)
  possedges<-n*n-n
  cat("Core search space ($adjacency) contains ", sum(object$adjacency), " edges out of ", possedges,
      "edges in a full search space","\n")
  if(is.list(object$tables[[1]])) {
    cat("Search space is extended", "\n")
  } else {
    cat("Search space is NOT extended", "\n")
  }
  
  nbl<-sum(object$blacklist)
  cat(nbl, " edges from the full space were blacklisted \n")
}



#' Summary of object of class 'orderMCMC'
#'
#' @param object object of class 'orderMCMC'
#' @param ... ignored
#' 
#' @rdname orderMCMC
#' @method summary orderMCMC
#' @export
summary.orderMCMC <- function(object, ...) {
  cat("object of class 'orderMCMC'")
  cat("\n\n")
  cat("Results:","\n")
  cat("maximum score DAG with", ncol(object$DAG), "nodes and ", sum(object$DAG)," edges:", "\n")
  cat("maximum DAG score=", object$score,"\n")
  cat("scores of samped DAGs (trace):")
  cat(object$trace[1], "...", object$trace[length(object$trace)])
  cat("\n\n")
  cat("MCMC settings:","\n")
  cat(paste("algorithm:",object$info$algo,"\n"))
  cat(paste("number of MCMC iterations:",object$info$iterations,"\n"))
  cat(paste("number of MCMC sampling steps (length of trace):",object$info$samplesteps,"\n"))
  cat(paste("initial search space:",object$info$spacealgo,"\n"))
  cat(paste("sample/MAP: ",object$info$sampletype,"\n"))
  cat("\n")
  cat("Additional output:","\n")
  cat("\n")
  if(!is.null(object$traceadd)) {
    cat(paste("traceadd contains",length(object$traceadd$incidence),"sampled DAGs","\n"))
  }
  if(!is.null(object$scoretable)) {
    cat("scoretable, object of class scorespace \n" )
  }
  cat("\n")
}


#' Summary of object of class 'iterativeMCMC'
#'
#' @param object object of class 'iterativeMCMC'
#' @param ... ignored
#' 
#' @rdname iterativeMCMC
#' @method summary iterativeMCMC
#' @export
summary.iterativeMCMC <- function(object, ...) {
  cat("object of class 'iterativeMCMC'")
  cat("\n\n")
  cat("Results:","\n")
  cat("maximum score DAG with", ncol(object$DAG), "nodes and ", sum(object$DAG)," edges: \n")
  cat("DAG score=", object$score,"\n")
  cat("\n")
  cat(paste("algorithm:",object$info$algo,"\n"))
  cat(paste("number of search space expansion steps:", length(object$maxtrace)-1,"\n"))
  cat(paste("number of edges in the intial search space:",sum(object$startspace),"\n"))
  cat(paste("number of added edges:",sum(object$endspace)-sum(object$startspace),"\n"))
  cat(paste("total number of MCMC iterations:",object$info$iterations*length(object$maxtrace),"\n"))
  cat(paste("total number of MCMC sampling steps (length of trace):",object$info$samplesteps*length(object$maxtrace),"\n"))
  cat(paste("number of MCMC iterations per expansion step:",object$info$iterations,"\n"))
  cat(paste("number of MCMC sampling steps per expansion step:",object$info$samplesteps,"\n"))
  cat(paste("initial search space:",object$info$spacealgo,"\n"))
  cat(paste("sample/MAP: ",object$info$sampletype,"\n"))
  cat("\n")
  cat("Additional output:","\n")
  if(!is.null(object$traceadd)) {
    cat(paste("traceadd contains",length(object$traceadd$incidence[[1]])*length(object$maxtrace),"sampled DAGs \n"))
    cat("\n")
  }
  if(!is.null(object$scoretable)) {
    cat("scoretable, object of class 'scorespace' \n")
  }
  cat("\n")
}


#' Summary of object of class 'partitionMCMC'
#'
#' @param object object of class 'partitionMCMC'
#' @param ... ignored
#' 
#' @rdname partitionMCMC
#' @method summary partitionMCMC
#' @export
summary.partitionMCMC <- function(object, ...) {
  cat("object of class 'partitionMCMC'")
  cat("\n\n")
  cat("Results:","\n")
  cat("maximum score DAG with", ncol(object$DAG), "nodes and ", sum(object$DAG)," edges:", "\n")
  cat("maximum DAG score=", object$score,"\n")
  cat("scores of samped DAGs (trace):")
  cat(object$trace[1], "...", object$trace[length(object$trace)])
  cat("\n\n")
  cat("MCMC settings:","\n")
  cat(paste("algorithm:",object$info$algo,"\n"))
  cat(paste("number of MCMC iterations:",object$info$iterations,"\n"))
  cat(paste("number of MCMC sampling steps (length of trace):",object$info$samplesteps,"\n"))
  cat(paste("initial search space:",object$info$spacealgo,"\n"))
  cat(paste("sample/MAP: ",object$info$sampletype,"\n"))
  cat("\n")
  cat("Additional output:","\n")
  cat("\n")
  if(!is.null(object$traceadd)) {
    cat(paste("traceadd contains",length(object$traceadd$incidence),"sampled DAGs","\n"))
  }
  if(!is.null(object$scoretable)) {
    cat("scoretable, object of class scorespace \n" )
  }
  cat("\n")
}



#' Summary of object of class 'itercomp'
#'
#' @param object object of class 'itercomp'
#' @param ... ignored
#' 
#' @rdname itercomp
#' @method summary itercomp
#' @export
summary.itercomp <-function(object, ...){
cat("object of class 'itercomp'")
cat("\n\n")
n<-nrow(object)
colo<-colnames(object)
if(n>1) {
  cat("structure fit changes: first -> last expansion iteration: \n")
  for(i in 1:ncol(object)) {
    cat(colo[i],":", object[1,i], "->",object[n,i],"\n")
  }
}
}


#' Summary of object of class 'samplecomp'
#'
#' @param object object of class 'samplecomp'
#' @param ... ignored
#' 
#' @rdname samplecomp
#' @method summary samplecomp
#' @export
summary.samplecomp <-function(object, ...){
  cat("object of class 'samplecomp'")
  cat("\n\n")
  n<-nrow(object)
  colo<-colnames(object)
  keymetrics<-c("TPR","FDR","SHD")
  if(n>1) {
    cat("best thresholds p for key metrics: \n")
  
      besttpr<-which.max(object[,"TPR"])[1]
      cat("TPR",":","p =",object[besttpr,"p"], "TPR =",object[besttpr,"TPR"],"\n")
      bestfdr<-which.min(object[,"FDR"])[1]
      cat("FDR",":","p =",object[bestfdr,"p"], "FDR =",object[bestfdr,"FDR"],"\n")
      bestshd<-which.min(object[,"SHD"])[1]
      cat("SHD",":","p =",object[bestshd,"p"], "SHD =",object[bestshd,"SHD"],"\n")

  }
}
