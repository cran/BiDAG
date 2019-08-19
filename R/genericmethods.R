setClass("MCMCtrace",
         slots = c(incidence = "list",
                   DAGscores = "list",
                   orderscores = "list",
                   orders= "list"))
setClass("MCMCmult",
         slots = c(incidence = "list",
                   DAGscores = "list",
                   orderscores = "list",
                   orders= "list"))

setClass("MCMCtracepart",
         slots = c(incidence = "list",
                   DAGscores = "list",
                   partitionscores = "list",
                   order= "list",
                   partiton="list"))
#'@export
#
print.MCMCmult <-function(x,...){
  cat("number of plus1 iterations", length(x), "\n")
  cat("number of saved DAGs", length(x)*length(x[[1]][[1]]), "\n")
  
}

#'@export
#
print.MCMCspace <-function(x,...){
  cat("Number of edges in the search space: ", sum(x$adjacency), "\n")
  if (all(lapply(x$scoretable,length)==1)) {
    cat("Search space is not extended\n")
  } else {
    cat("Search space is extended (plus1)\n")
  }
}


#'@export
#
print.MCMCmax <-function(x,...){
  
  cat("Number of edges in the maximum score DAG: ", sum(x$DAG), "\n")
  cat("Maximum DAG score: ", x$score, "\n")
  cat("Order: ", x$order, "\n")
  
}


setMethod("plot", signature(x = "MCMCmult"),
          function(x, y, ...){
            nchains<-length(x$DAGscores)
            scorevecmin<-unlist(x$DAGscores[[1]])
            scorevecprevmax<-unlist(x$DAGscores[[nchains-1]])
            minprev<-min(scorevecprevmax)
            scorevecmax<-unlist(x$DAGscores[[nchains]])
            vecl<-length(scorevecmin)
            # burnin<-ceiling(vecl*0.2)
            # score20<-min(scorevecmax[burnin:vecl])
            scoremax<-max(scorevecmax)
            scoremin<-min(scorevecmin)
            scoremaxmin<-min(scorevecmax)
            par(mfrow=c(1,2))
            plot(scorevecmin,type="l",col=1,xlab="iteration",ylab="logscore",
                 ylim=c(scoremin,scoremax),main="DAG scores: all plus iterations",cex.main=0.7)
            for(i in 2:nchains) {
              lines(unlist(x$DAGscores[[i]]),col=i)
            }
            plot(c(scorevecprevmax,scorevecmax),type="l",col="blue",xlab="iteration",ylab="logscore",
                 ylim=c(minprev,scoremax),main="last 2 plus1 iterations",cex.main=0.7)
            
          })


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
  cat("SHD: ", x[1], "\n")
  cat("TP: ", x[2], "\n")
  cat("FP: ", x[3], "\n")
}

#'@export
#a standard print method for function MCMCtrace
print.MCMCtrace <-function(x,...){
  cat("Number of saved steps: ", length(x$incidence), "\n")
}


#'@export
#a standard print method for function MCMCtrace
print.MCMCres <-function(x,...){
print(x$max)
if(!is.null(x$chain)) {
  print(x$chain)
}

if(!is.null(x$space)) {
  print(x$space)
}
  
}

