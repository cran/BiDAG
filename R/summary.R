#' S3 methods for class 'scoreparameters'
#' 
#' Summary of object of class 'scoreparameters'
#'
#' @param object object of class 'scoreparameters'
#' @param ... ignored
#' 
# @rdname scoreparameters
#' @method summary scoreparameters
#' @export summary.scoreparameters
#' @export
summary.scoreparameters <-function(object,...){
  cat("object of class 'scoreparameters' \n")
  cat("number of nodes (variables):", object$n, "\n")
  cat("number of observations:", nrow(object$data), "\n")
  cat("score type:", object$type, "\n")
  
  if(object$DBN==TRUE) {
    cat("score object created for a DBN \n")
    if(object$bgn>0) {
      cat("static nodes (no parents are allowed):", object$bgnodes, "\n")
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



#' S3 methods for class 'MCMCscoretab'
#' 
#' Summary of object of class 'MCMCscoretab'
#'
#' @param object object of class 'MCMCscoretab'
#' @param ... ignored
#' 
# @rdname MCMCspace
#' @method summary MCMCscoretab
#' @export summary.MCMCscoretab
#' @export
summary.MCMCscoretab <-function(object,...){
  n<-ncol(object$adjacency)
  possedges<-n*n-n
  cat("Core search space ($adjacency) contains ", sum(object$adjacency), " edges out of ", possedges,
      "edges in a full search space","\n")
  if(length(object$scoretable[[1]])>1) {
    cat("Search space is extended", "\n")
  } else {
    cat("Search space is NOT extended", "\n")
  }
}


#' S3 methods for class 'MCMCres'
#' 
#' Summary of object of class 'MCMCres'
#'
#' @param object object of class 'MCMCres'
#' @param ... ignored
#' 
# @rdname MCMCres
#' @method summary MCMCres
#' @export summary.MCMCres
#' @export
summary.MCMCres <- function(object, ...) {
  cat("Results:","\n")
  cat("\n")
  cat("maximum score DAG with", ncol(object$DAG), "nodes and ", sum(object$DAG)," edges:", "\n")
  cat("DAG score=", object$score,"\n")
  cat("\n")
  cat("MCMC settings:","\n")
  cat("\n")
  cat(paste("algorithm:",object$info$algo,"\n"))
  cat(paste("number of MCMC iterations:",object$info$iterations,"\n"))
  cat(paste("number of MCMC sampling steps (length of trace):",object$info$samplesteps,"\n"))
  cat(paste("initial search space:",object$info$spacealgo,"\n"))
  cat(paste("sample/MAP: ",object$info$sampletype,"\n"))
  cat("\n")
  cat("Additional output:","\n")
  cat("\n")
  if(is.null(object$traceadd)) {
    cat("sampled matrices were not saved","\n")
    cat("to save sampled matrices set 'chainout' to TRUE","\n")
    
  } else {
    cat(paste("Additional output ($traceadd) contains",length(object$traceadd$incidence),"sampled DAGs","\n"))
  }
  if(!is.null(object$scoretable)) {
    cat("Additional output ($scoretable) includes the adjacency matrix of the core search space and the score tables corresponding to this search space","\n")
  }
  cat("\n")
}

#' S3 methods for class 'MCMCmult'
#' 
#' Summary of object of class 'MCMCmult'
#'
#' @param object object of class 'MCMCmult'
#' @param ... ignored
#' 
# @rdname MCMCmult
#' @method summary MCMCmult
#' @export summary.MCMCmult
#' @export
summary.MCMCmult <- function(object, ...) {
  cat("Results:","\n")
  cat("\n")
  cat("maximum score DAG with", ncol(object$DAG), "nodes and ", sum(object$DAG)," edges:", "\n")
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
  cat("\n")
  if(is.null(object$chain)) {
    cat("sampled matrices were not saved","\n")
    cat("to save sampled matrices set 'chainout' to TRUE","\n")
    
  } else {
    cat(paste("$chain contains",length(object$traceadd$incidence)*length(object$maxtrace),"sampled DAGs","\n"))
  }
  cat("\n")
  if(!is.null(object$scoretable)) {
    cat("$scoretable includes the adjacency matrix of the core search space and the score tables corresponding to this search space","\n")
  }
  cat("\n")
}
