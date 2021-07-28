#'Deriving an adjecency matrix of a full DBN
#'
#'This function transforms a compact 2-slice adjacency matrix of DBN into full T-slice adjacency matrix
#'
#' @param DBNmat a square matrix, representing initial and transitional structure of a DBN; the size of matrix is 2*dyn+b
#' @param slices integer, number of slices in an unrolled DBN
#' @param b integer, number of static variables
#' @return an adjacency matrix of an unrolled DBN
#' @examples
#' compact2full(DBNmat, slices=5, b=3)
#' @export
compact2full<-function(DBNmat,slices,b=0) {
  dyn<-(ncol(DBNmat)-b)/2
  if(slices<3) {
    return(DBNmat)
  } else {

    if(all(is.character(colnames(DBNmat)))){
      baseall<-colnames(DBNmat)
      basenames<-colnames(DBNmat)[1:dyn+b]
    } else {
      if(b!=0) {
        staticnames<-paste("s",1:b,sep="")
        basenames<-paste("v",1:dyn,sep="")
        baseall<-c(staticnames,basenames)
      } else {
        basenames<-paste("v",1:dyn,sep="")
        baseall<-basenames
      }
    }
    
    for(i in 3:slices) {
    baseall<-c(baseall,paste(basenames,".",i,sep=""))
    }
    
    nbig<-slices*dyn+b
    DBNbig<-matrix(0,nrow=nbig,ncol=nbig)
    print(baseall)
    colnames(DBNbig)<-baseall
    rownames(DBNbig)<-baseall
    
    
    DBNbig[1:(dyn+b),1:dyn+b]<-DBNmat[1:(dyn+b),1:dyn+b] #copying initial structure
   
    intStruct<-DBNmat[1:dyn+dyn+b,1:dyn+dyn+b] #internal structure
    transStruct<-DBNmat[1:dyn+b,1:dyn+dyn+b] #transitional structure
    if(b>0) {
      bgStrct<-DBNmat[1:b,1:dyn+dyn+b] #edges from static variables
    }
    for(i in 1:(slices-1)) {
      if(b>0) {
        DBNbig[1:b,1:dyn+i*dyn+b]<-bgStrct
      }
      DBNbig[1:dyn+(i-1)*dyn+b,1:dyn+i*dyn+b]<-transStruct
      DBNbig[1:dyn+i*dyn+b,1:dyn+i*dyn+b]<-intStruct
    }
    return(DBNbig)
  }
}
#'Deriving a compact adjacency matrix of a DBN 
#'
#'This function transforms an unrolled adjacency matrix of DBN into a compact representation
#'
#' @param DBNmat a square matrix, representing the structure of an unrolled DBN; the size of matrix is slices*dyn+b; all static variables are assumed to be in the first b rows and columns of the matrix
#' @param b integer, number of static variables; 0 by default
#' @examples
#' full2compact(DBNunrolled,b=3)
#'@export
full2compact<-function(DBNmat,b=0) {
   dyn<-(ncol(DBNmat)-b)/2
    DBNcompact<-DBNmat[1:(2*dyn+b),1:(2*dyn+b)]
    return(DBNcompact)
}

#turns internal representation into user-friendly
DBNtransform<-function(DBNmat,param) {
  newDBNmat<-matrix(0,nrow=param$n+param$nsmall,ncol=param$n+param$nsmall)
  colnames(newDBNmat)<-param$labels.short
  rownames(newDBNmat)<-param$labels.short
  newDBNmat[param$usrinitstr$rows,param$usrinitstr$cols]<-DBNmat[param$intstr$rows,param$intstr$cols]
  newDBNmat[param$usrintstr$rows,param$usrintstr$cols]<-DBNmat[param$intstr$rows,param$intstr$cols]
  newDBNmat[param$usrtrans$rows,param$usrtrans$cols]<-DBNmat[param$trans$rows,param$trans$cols]
  return(newDBNmat)
}
#turns internal representation into user-friendly
DBNtransform.init<-function(DBNmat,param) {
  if(param$bgn>0) {
      newDBNmat<-matrix(0,nrow=param$bgn+param$nsmall,ncol=param$bgn+param$nsmall)
      colnames(newDBNmat)<-param$labels.short[1:param$n]
      rownames(newDBNmat)<-param$labels.short[1:param$n]
      newDBNmat[,1:param$bgn]<-DBNmat[,1:param$bgn+param$nsmall]
      newDBNmat[,1:param$nsmall+param$bgn]<-DBNmat[,1:param$nsmall]
      DBNmat<-newDBNmat
      newDBNmat[1:param$bgn,]<-DBNmat[1:param$bgn+param$nsmall,]
      newDBNmat[1:param$nsmall+param$bgn,]<-DBNmat[1:param$nsmall,]
      return(newDBNmat)
  } else {
    return(DBNmat)
  }
}
#turns user-friendly representation into internal
DBNbacktransform<-function(DBNmat,param) {
  if(!is.null(colnames(DBNmat))) {
    oldnodelabels<-colnames(DBNmat)
    newnodelabels<-oldnodelabels
    newnodelabels[param$intstr$cols]<-oldnodelabels[param$usrtrans$cols]
  }
  newDBNmat<-matrix(0,nrow=param$n+param$nsmall,ncol=param$n+param$nsmall)
  newDBNmat[param$intstr$rows,param$intstr$cols]<-1*(DBNmat[param$usrintstr$rows,param$usrintstr$cols]|DBNmat[param$usrinitstr$rows,param$usrinitstr$cols])
  newDBNmat[param$trans$rows,param$trans$cols]<-DBNmat[param$usrtrans$rows,param$usrtrans$cols]
  if(!param$split) {
    return(newDBNmat)
  } else {
    res<-list()
    initDBNmat<-DBNmat[1:param$n,1:param$n]
    newinitDBNmat<-DBNmat[1:param$n,1:param$n]
    
    newinitDBNmat[,1:param$bgn+param$nsmall]<-initDBNmat[,1:param$bgn]
    newinitDBNmat[,1:param$nsmall]<-initDBNmat[,1:param$nsmall+param$bgn]
    initDBNmat<-newinitDBNmat
    newinitDBNmat[1:param$bgn+param$nsmall,]<-initDBNmat[1:param$bgn,]
    newinitDBNmat[1:param$nsmall,]<-initDBNmat[1:param$nsmall+param$bgn,]
    res$init<-newinitDBNmat
    
    transDBNmat<-matrix(0,nrow=param$otherslices$n,ncol=param$otherslices$n)
    DBNmat<-DBNcut(DBNmat,dyn=param$nsmall,b=param$bgn)
    transDBNmat[param$intstr$rows,param$intstr$cols]<-DBNmat[param$usrintstr$rows,param$usrintstr$cols]
    transDBNmat[param$trans$rows,param$trans$cols]<-DBNmat[param$usrtrans$rows,param$usrtrans$cols]
    res$trans<-transDBNmat
    
    return(res)
  }
}
DBNcut<-function(adj,dyn,b){
  adj[,1:(dyn+b)]<-0
  return(adj)
}
DBNinit<-function(adj,dyn,b){
  adj<-adj[1:(b+dyn),1:(b+dyn)]
  if(b>0) {
    adj[,1:b]<-0
  }
  return(adj)
}

#Combining initial and transition DBN structures in one matrix
mergeDBNstr<-function(initStruct,transStruct) {
  n<-ncol(initStruct)
  if(!is.matrix(initStruct)) {
    initStruct<-graph2m(initStruct)
  }
  if(!is.matrix(transStruct)) {
    transStruct<-graph2m(transStruct)
  }
  n<-ncol(initStruct)
  transStruct[1:n,1:n]<-initStruct
  return(transStruct)
}

#Combining orders for a DBN
mergeDBNord<-function(initorder,transorder) {
 return(c(initorder,transorder))
}

#Combining order scores for a DBN
mergeDBNscore<-function(initscore,transscore) {
  return(initscore+transscore)
}

#this function produces common result for DBN structure learning when samestruct=FALSE
mergeDBNres<-function(result.init,result.trans,scorepar,algo) {
  
  res<-list()
  
  maxtrans<-DBNtransform(result.trans$DAG,scorepar)
  maxinit<-DBNtransform.init(result.init$DAG,scorepar)
  res$DAG<-mergeDBNstr(maxinit,maxtrans)
  res$order<-mergeDBNord(result.init$order,result.trans$order)
  res$score<-mergeDBNscore(result.init$score,result.trans$score)
  
  if(!is.null(result.init$traceadd)) {
    
    result.init$traceadd$incidence<-lapply(result.init$traceadd$incidence,DBNtransform.init,param=scorepar)
    result.trans$traceadd$incidence<-lapply(result.trans$traceadd$incidence,DBNtransform,param=scorepar)
    result.trans$traceadd$incidence<-lapply(result.trans$traceadd$incidence,DBNcut,dyn=scorepar$nsmall,b=scorepar$bgn)
    res$traceadd$incidence<-mapply(mergeDBNstr,result.init$traceadd$incidence,result.trans$traceadd$incidence,SIMPLIFY = FALSE)
    res$trace<-mapply(mergeDBNscore,result.init$trace,result.trans$trace)
    
  if(algo=="order"){
    res$traceadd$orders<-mapply(mergeDBNord,result.init$traceadd$orders,result.trans$traceadd$orders,SIMPLIFY = FALSE)
    res$traceadd$orderscores<-mapply(mergeDBNscore,result.init$traceadd$orderscores,result.trans$traceadd$orderscores)
  } else if (algo=="partition") {
    res$traceadd$order<-mapply(mergeDBNord,result.init$traceadd$order,result.trans$traceadd$order,SIMPLIFY = FALSE)
    res$traceadd$partitionscores<-mapply(mergeDBNscore,result.init$traceadd$partitionscores,result.trans$traceadd$partitionscores)
  }
  }
  
  attr(res,"class")<-"MCMCres"
  
  return(res)
}


#this function produces common result for DBN iterative structure learning when samestruct=FALSE
mergeDBNres.it<-function(result.init,result.trans,scorepar) {
  
  res<-list()
  
  res$init<-result.init
  res$trans<-result.trans
  
  maxtrans<-DBNtransform(result.trans$DAG,scorepar)
  maxinit<-DBNtransform.init(result.init$DAG,scorepar)
  
  
  for(i in 1:length(res$trans$maxtrace)) {
    res$trans$maxtrace[[i]]$DAG<-DBNtransform(res$trans$maxtrace[[i]]$DAG,scorepar)
    res$trans$maxtrace[[i]]$DAG<-DBNcut(res$trans$maxtrace[[i]]$DAG,dyn=scorepar$nsmall,b=scorepar$bgn)
  }
  
  for(i in 1:length(res$init$maxtrace)) {
    res$init$maxtrace[[i]]$DAG<-DBNtransform.init(res$init$maxtrace[[i]]$DAG,scorepar)
    res$init$maxtrace[[i]]$DAG<-DBNinit(res$init$maxtrace[[i]]$DAG,dyn=scorepar$nsmall,b=scorepar$bgn)
  }
  
  res$DAG<-mergeDBNstr(maxinit,maxtrans)
  res$order<-mergeDBNord(result.init$order,result.trans$order)
  res$score<-mergeDBNscore(result.init$score,result.trans$score)

  endinit<-DBNtransform.init(result.init$endspace,scorepar)
  endtrans<-DBNtransform(result.trans$endspace,scorepar)
  startinit<-DBNtransform.init(result.init$startspace,scorepar)
  starttrans<-DBNtransform(result.trans$startspace,scorepar)
    
  res$endspace<-mergeDBNstr(endinit,endtrans)
  res$startspace<-mergeDBNstr(startinit,starttrans)
  
  return(res)
}






