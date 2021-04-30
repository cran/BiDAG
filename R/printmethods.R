##############################
#print methods for S3 classes#
#   orderMCMC                #
#   partitionMCMC            #
#   iterativeMCMC            #
#   scoreparameters          #
#   scorespace               #
#   intersim                 #
#   sampsim                  #
##############################

#' Print object of class 'scoreparameters'
#'
#' @param x object of class 'scoreparameters'
#' @param ... ignored
#' 
#' @rdname scoreparameters
#' @export
print.scoreparameters <-function(x, ...){
  cat("$type \n")
  cat(x$type,"\n")
  cat("$data \n")
  if(!x$DBN) {
    cat("data contains",nrow(x$data),"observations of",ncol(x$data),"variables \n" )
  } else {
    cat("data contains",nrow(x$data),"rows and",ncol(x$data),"columns \n") 
  }
  cat("$DBN \n")
  cat(x$DBN,"\n")
  if (x$DBN) {
    cat("$slices \n")
    cat(x$slices,"\n")
    if(x$bgn>0) {
      cat("$static \n")
      cat(x$static,"\n")
      cat("$firstslice \n")
      cat("...\n")
      cat("$otherslices \n")
      cat("...\n\n")
    }
  } else {
    if(x$type=="bge") {
      if(x$bgn>0) {
        cat("$bgnodes \n")
        cat(x$bgnodes,"\n")
      } 
      cat("$am \n")
      cat(x$am,"\n")
      cat("$aw \n")
      cat(x$am,"\n")
      cat("$means \n")
      cat(x$means[1],"...", x$means[length(x$means)] ,"\n")
      cat("$SigmaN \n")
      cat("matrix", ncol(x$SigmaN), "x",ncol(x$SigmaN), "\n")
      
    } else if (x$type%in%c("bde","bdecat")) {
      cat("$chi \n")
      cat(x$chi,"\n\n")
      cat("$pf \n")
      cat(x$pf,"\n\n")
      if(x$bgn>0) {
        cat("$bgnodes \n")
        cat(x$bgnodes,"\n\n")
      }
    } 
  }
  
  
  if(!is.null(x$weightvector)) {
    cat("$weightvector \n")
    cat("...\n") 
    cat("...numeric vector of length",length(x$weightvector),"\n") 
  }

}


#' Prints object of class 'orderMCMC'
#'
#' @param x object of class 'orderMCMC'
#' @param ... ignored
#' 
#' @rdname orderMCMC
#' @export
print.orderMCMC <-function(x, ...){
  cat("object of class 'orderMCMC', from Call:", "\n")
  print(x$info$fncall)
  cat("\n")
  cat("$DAG\n")
  cat("adjacency matrix of a DAG with", ncol(x$DAG), "nodes and ", sum(x$DAG)," edges", "\n")
  cat("$score\n")
  cat(x$score,"\n")
  cat("$maxorder\n")
  cat(x$maxorder[1],"...",x$maxorder[length(x$maxorder)],"\n")
  cat("$info\n")
  cat("... \n")
  cat("$trace\n")
  cat(x$trace[1], "...", x$trace[length(x$trace)])
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

#' Prints object of class 'partitionMCMC'
#'
#' @param x object of class 'partitionMCMC'
#' @param ... ignored
#' 
#' @rdname partitionMCMC
#' @export
print.partitionMCMC <-function(x, ...){
  cat("object of class 'partitionMCMC', from Call:", "\n")
  print(x$info$fncall)
  cat("\n")
  cat("$DAG\n")
  cat("adjacency matrix of a DAG with", ncol(x$DAG), "nodes and ", sum(x$DAG)," edges", "\n")
  cat("$score\n")
  cat(x$score,"\n")
  cat("$maxorder\n")
  cat(x$maxorder[1],"...",x$maxorder[length(x$maxorder)],"\n")
  cat("$info\n")
  cat("... \n")
  cat("$trace\n")
  cat(x$trace[1], "...", x$trace[length(x$trace)])
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


#' Prints object of class 'iterativeMCMC'
#'
#' @param x object of class 'iterativeMCMC'
#' @param ... ignored
#' 
#' @rdname iterativeMCMC
#' @export
print.iterativeMCMC <-function(x, ...){
  cat("object of class 'iterativeMCMC', from Call:", "\n")
  print(x$info$fncall)
  cat("\n")
  cat("$DAG\n")
  cat("adjacency matrix of a DAG with", ncol(x$DAG), "nodes and ", sum(x$DAG)," edges", "\n")
  cat("$maxorder\n")
  cat(x$maxorder[1],"...",x$maxorder[length(x$maxorder)],"\n")
  cat("$score\n")
  cat(x$score,"\n")
  cat("$maxtrace:\n")
  cat("local maximums at each expansion step:\n")
  cat(unlist(lapply(x$maxtrace,function(x)x$score)),"\n")
  cat("$info\n")
  cat("... \n")
  cat("$trace\n")
  cat(x$trace[[1]][1], "...", x$trace[[length(x$trace)]][length(x$trace[[1]][1])])
  cat("\n")
  if(!is.null(x$traceadd)) {
    cat("$traceadd\n")
    cat("adjacency matrices of sampled DAGs, corresponding orders and order scores\n")
    cat("\n")
  }
  if(!is.null(x$scoretable)) {
      cat("$scoretable\n")
      cat("scoretable, object of class 'scorespace'\n")
  }
}

#' Prints 'scorespace' object
#'
#' @param x object of class 'scorespace'
#' @param ... ignored
#' 
#' @rdname scorespace
#' @export
print.scorespace <-function(x, ...){
  cat("object of class 'scorespace'")
  n<-ncol(x$adjacency)
  cat("\n")
  cat("$adjacency \n")
  cat("matrix", n,"x",n, "\n\n")
  n<-length(x$tables)
  if(is.list(x$tables[[1]])) {
    cat("$tables", "\n")
    cat("[[1]][[1]]", "\n")
    cat(x$tables[[1]][[1]][1,],"\n")
    cat("...", "\n")
    lastt<-length(x$tables[[n]])
    cat("[[",n,"]][[",lastt,"]]\n", sep="")
    cat(x$tables[[n]][[lastt]][1,],"\n")
    cat("...", "\n")
    
  } else {
    cat("$tables[[1]]", "\n")
    cat(x$tables[[1]][1,], "...\n")
    cat("...\n")
    cat("$tables[[", length(x$tables),"]]\n",sep="")
    cat(x$tables[[length(x$tables)]][1,],"...\n")
  }
 
  cat("\n")
  
  cat("$blacklist \n")
  cat("matrix", n,"x",n, "\n")
}

#' Prints itercomp object.
#'
#' @param x object of class 'itercomp'
#' @param ... ignored
#' @rdname itercomp
#' @export 
print.itercomp <-function(x, ...){
  cat("object of class 'itercomp'")
  cat("\n\n")
  nit<-nrow(x)
  if(nit<5) {
    print(x[1:nit,,drop=FALSE])
  } else {
    print(x[1:2,])
    cat("... \n")
    print(x[c(nit-1,nit),])
  }
}

#' Prints samplecomp object.
#' 
#' @param x object of class 'samplecomp'
#' @param ... ignored
#' @rdname samplecomp
#' @export
print.samplecomp <-function(x, ...){
  cat("object of class 'samplecomp'")
  cat("\n\n")
  np<-nrow(x)
  if(is.null(np)) {
    print(x[1:length(x)])
  } else if(np<10) {
    print(x[1:np,,drop=FALSE])
  } else {
    print(x[1:5,])
    cat("... \n")
    print(x[c(np-4,np),])
  }
}



