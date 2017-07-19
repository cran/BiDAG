#'@export
#a generic method (does not need description) print for scoreparameters class
print.scoreparameters <-function(x,...){
  if (x$type=="bge") {
    cat("Scoretype: BGe", "\n")
  } else {
    cat("Scoretype: BDe", "\n")
    cat("Prior pseudo counts:", x$chi,"\n")
    cat("Edges prenalization factor:", x$pf,"\n")
  }
  
  if(is.null(x$weightvector)) {
    cat("Data is not weighted\n")
  } else {
    "Data is weighted according to the weightvector"
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
#
print.MCMCchain <-function(x,...){
  cat("Number of saved steps: ", length(x$incidence), "\n")
}

#'@export
#
print.MCMCmult <-function(x,...){
  
}

#'@export
#
print.MCMCspace <-function(x,...){
  cat("Number of edges in the searchspace: ", sum(x$adjacency), "\n")
  if (all(lapply(x$scoretable,length)==1)) {
    cat("Search space is not extended (not plus1)\n")
  } else {
    cat("Search space is extended (plus1)\n")
  }
}


#'@export
#
print.MCMCmax <-function(x,...){
  
  cat("Number of edges in the maximum score DAG: ", sum(x$DAG), "\n")
  cat("maximum DAG score: ", x$score, "\n")
  cat("Order: ", x$order, "\n")
  
}
