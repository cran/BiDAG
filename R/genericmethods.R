#this script contains definition of classes in BiDAG
#and all SUMMARY AND PRINT functions for these classes

setClass("MCMCmax",
         slots = c(DAG = "matrix",
                   score = "numeric",
                   reach = "integer",
                   order= "vector"))
setClass("MCMCtrace",
         slots = c(incidence = "list",
                   DAGscores = "list",
                   partitionscores = "list",
                   orderscores = "list",
                   order= "list",
                   partiton="list"))
setClass("MCMCspace",
         slots = c(adjacency = "matrix",
                   scoretables = "list"))

setClass("MCMCres",
         slots = c(max = "MCMCmax",
                   trace = "MCMCtrace",
                   space = "MCMCspace",
                   info="list"))
setClass("MCMCmult",
         slots = c(max = "MCMCmax",
                   maxtrace = "list",
                   trace = "MCMCtrace",
                   space = "MCMCspace",
                   info="list"))

#SUMMARY METHODS

setMethod("summary", "MCMCres",
          function(object, ...) {
            cat(paste("algorithm:",object$info$algo,"\n"))
            cat(paste("number of MCMC iterations:",object$info$iterations,"\n"))
            cat(paste("initial search space:",object$info$spacealgo,"\n"))
            cat(paste("sample/MAP: ",object$info$sampletype,"\n"))
            cat(paste("number of sampling steps:",object$info$samplesteps,"\n"))
            cat(paste("score of the maximum score DAG:",round(object$max$score,2),"\n"))
            cat(paste("maximum reached at:",object$max$reach, "step out of", object$info$samplesteps,"\n"))
            
            if(is.null(object$chain)) {
              cat("MCMC trace was not saved","\n")
            } else {
              cat(paste("MCMC trace contains:",length(object$chain$incidence),"saved steps","\n"))
            }
            if(!is.null(object$space)) {
              cat("additional objects returned: the search space and the score tables","\n")
            }
          })


setMethod("summary", "MCMCmult",
          function(object, ...) {
            cat("    MCMC settings\n")
            cat(paste("algorithm:",object$info$algo,"\n"))
            cat(paste("sample/MAP: ",object$info$sampletype,"\n"))
            cat(paste("initial search space:",object$info$spacealgo,"\n"))
            cat(paste("number of MCMC iterations:",object$info$iterations,"\n"))
            cat(paste("number of sample steps:",object$info$samplesteps,"\n"))
            cat("\n")
            
            cat("   Search space optimization\n")
            cat(paste("number of space expansion steps:",length(object$maxtrace),"\n"))
            cat(paste("number of edges in the intial search space:",sum(object$startspace),"\n"))
            cat(paste("number of added edges:",sum(object$endspace)-sum(object$startspace),"\n"))
            cat("\n")
            
            cat("   MAP estimate\n")
            cat(paste("score of the maximum score DAG:",round(object$max$score,2),"\n"))
            cat("\n")
            
            cat("   Additional objects\n")
            if(is.null(object$chain)) {
            } else {
              cat(paste("MCMC trace, contains ",length(object$maxtrace)," x ",length(object$chain$incidence[[1]]),"saved steps","\n"))
            }
            if(!is.null(object$scoretable)) {
              cat("score tables","\n")
            }
          })


#PRINT methods



#'@export
#
print.MCMCmax <-function(x,...){
  
  cat("Number of edges in the maximum score DAG: ", sum(x$DAG), "\n")
  cat("Maximum DAG score: ", x$score, "\n")
  cat("Order: ", x$order, "\n")
  
}

#'@export
#a generic method (does not need description) print for scoreparameters class
print.scoreparameters <-function(x,...){
  if (x$type=="bge") {
    cat("Score type: BGe", "\n")
  } else if (x$type=="bde") {
    cat("Score type: BDe", "\n")
    cat("Prior pseudo counts:", x$chi,"\n")
    cat("Edge penalization factor:", x$pf,"\n")
  } else if (x$type=="bdecat") {
    cat("Score type: BDe for categorical data", "\n")
    cat("Prior pseudo counts:", x$chi,"\n")
    cat("Edge penalization factor:", x$pf,"\n")
  } else if (x$type=="dbn") {
    cat("Score type: DBN", "\n")
    cat("DBN score type: ", "\n")
    } else {
    cat("Score type: user defined", "\n")
  }
  if(is.null(x$weightvector)) {
    cat("Data is not weighted\n")
  } else {
    "Data is weighted according to the weight vector"
  }
  cat("Score constant vector:", x$scoreconstvec,"\n")
}

#'@export
#a standard print method for function compareDAGs
print.compDAGs<-function(x,...) {
  cat("TP: ", x[2], "\n")
  cat("TPR: ", x[5], "\n")
  cat("SHD: ", x[1], "\n")
  cat("FP: ", x[3], "\n")
  cat("FN: ", x[4], "\n")
}

#'@export
#a standard print method for function MCMCtrace
print.MCMCtrace <-function(x,...){
  cat("Number of saved steps: ", length(x$incidence), "\n")
}





