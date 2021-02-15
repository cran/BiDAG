##############################
#print methods for S3 classes#
#   MCMCres                  #
#   MCMCmult                 #
#   scoreparameters          #
#   MCMCspace                #
##############################

#' S3 methods for class 'scoreparameters'
#' 
#' Print object of class 'scoreparameters'
#'
#' @param x object of class 'scoreparameters'
#' @param ... ignored
#' 
# @rdname scoreparameters
#' @export print.scoreparameters
#' @export
print.scoreparameters <-function(x,...){
  cat("$type \n")
  cat(x$type,"\n\n")
  
  cat("$labels.short \n")
  cat(x$labels.short,"\n\n")
  
  cat("$data \n")
  cat("...\n")
  cat("data contains",nrow(x$data),"observations of",ncol(x$data),"variables \n" )
  
  cat("$DBN \n")
  cat(x$DBN,"\n\n")
  
  cat("$MDAG \n")
  cat(x$MDAG, "\n\n")
  
  if(x$type=="bge") {
    cat("$am \n")
    cat(x$am,"\n\n")
    cat("$aw \n")
    cat(x$am,"\n\n")
    if(x$bgn>0) {
      cat("$bgnodes \n")
      cat(x$bgnodes,"\n\n")
    }
  } else if (x$type%in%c("bde","bdecat")) {
    cat("$chi \n")
    cat(x$chi,"\n\n")
    cat("$pf \n")
    cat(x$pf,"\n\n")
    if(x$bgn>0) {
      cat("$bgnodes \n")
      cat(x$bgnodes,"\n\n")
    }
  } else if (x$type=="dbn") {
    cat("$slices \n")
    cat(x$slices,"\n\n")
    
    if(x$bgn>0) {
      cat("$static \n")
      cat(x$static,"\n\n")
    }

  }
  
  if(!is.null(x$weightvector)) {
    cat("$weightvector \n")
    cat("...\n") 
    cat("...numeric vector of length",length(x$weightvector),"\n") 
  }

}



#' S3 methods for class 'MCMCres'
#' 
#' Prints object of class 'MCMCres'
#'
#' @param x object of class 'MCMCres'
#' @param ... ignored
#' 
# @rdname MCMCres
#' @export print.MCMCres
#' @export
print.MCMCres <-function(x,...){
  cat("object of class 'MCMCres', contains the result of running MCMC", "\n")
  cat("\n")
  cat("$DAG\n")
  cat("adjacency matrix of a DAG with", ncol(x$DAG), "nodes and ", sum(x$DAG)," edges", "\n")
  cat("\n")
  cat("$score\n")
  cat(x$score,"\n")
  cat("\n")
  cat("$maxorder\n")
  cat(x$maxorder[1],"...",x$maxorder[length(x$maxorder)],"\n")
  cat("\n")
  cat("$info\n")
  print(x$info)
  cat("\n")
  cat("$trace\n")
  cat("score trace of sampled DAGs\n")
  cat("\n")
  if(!is.null(x$traceadd)) {
    cat("$traceadd\n")
    cat("adjacency matrices of sampled DAGs, corresponding orders and order scores\n")
    cat("\n")
  }
  if(!is.null(x$scoretable)) {
    cat("$scoretable\n")
    cat("score tables corresponding to core search space $endspace\n")
  }
}

#' S3 methods for class 'MCMCmult'
#' 
#' Prints object of class 'MCMCmult'
#'
#' @param x object of class 'MCMCmult'
#' @param ... ignored
#' 
# @rdname MCMCmult
#' @export print.MCMCmult
#' @export
print.MCMCmult <-function(x,...){
  cat("object of class 'MCMCmult', contains the result of running iterative MCMC", "\n")
  cat("\n")
  cat("$DAG\n")
  cat("adjacency matrix of a DAG with", ncol(x$DAG), "nodes and ", sum(x$DAG)," edges", "\n")
  cat("\n")
  cat("$maxorder\n")
  cat(x$maxorder[1],"...",x$maxorder[length(x$maxorder)],"\n")
  cat("$score\n")
  cat(x$score,"\n")
  cat("\n")
  cat("$maxtrace:\n")
  cat("local maximums at each expansion step:", unlist(lapply(x$maxtrace,function(x)x$score)))
  cat("\n")
  cat("\n")
  cat("$info\n")
  print(x$info)
  cat("\n")
  cat("$trace\n")
  cat("score trace of sampled DAGs\n")
  cat("\n")
  if(!is.null(x$traceadd)) {
    cat("$traceadd\n")
    cat("adjacency matrices of sampled DAGs, corresponding orders and order scores\n")
    cat("\n")
  }
  if(!is.null(x$scoretable)) {
      cat("$scoretable\n")
      cat("score tables corresponding to core search space $endspace\n")
  }
}

#' S3 methods for class 'MCMCscoretab'
#' 
#' Prints object of class 'MCMCscoretab'
#'
#' @param x object of class 'MCMCscoretab'
#' @param ... ignored
#' 
# @rdname MCMCscoretab
#' @export print.MCMCscoretab
#' @export
print.MCMCscoretab <-function(x,...){
  cat("$adjacency, adjacency matrix representing a core search space:", "\n")
  cat("\n")
  print(x$adjacency)
  cat("\n")
  cat("$tables, a list containing computed local scores used in MCMC chain:", "\n")
  cat("\n")
  if(length(x$tables[[1]])>1) {
    cat("$tables[[1]][[1]]", "\n")
    cat("\n")
    cat(x$tables[[1]][[1]][1,])
    cat("...", "\n")
  } else {
    cat("$tables[[1]]", "\n")
    cat("\n")
    cat(x$tables[[1]][1,])
    cat("...", "\n")
  }
 
  cat("\n")
}


