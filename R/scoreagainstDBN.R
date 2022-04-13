#'Score against DBN
#'
#'Scoring observations against a DBN structure
#'
#'@param scorepar object of class 'scoreparameters'
#'@param incidence adjacency matrix of a DAG
#'@param datatoscore matrix or vector containing observations to be scored
#'@param marginalise (logical) should marginal score be used?
#'@param onlymain (logical) should static nodes be included in the score?
#'@param datainit optional, in case of unbalanced design, the mean score of available samples for T0 are computed
#'@return vector of log-scores
#'@author Polina Suter
#'@export
scoreagainstDBN<-function(scorepar, incidence, datatoscore=NULL,
                         marginalise=FALSE, onlymain=FALSE, datainit=NULL){
  backDBN<-DBNbacktransform_l(incidence,scorepar,coln=TRUE)
  if(scorepar$split) {
    initMat<-backDBN$init
    transMat<-backDBN$trans
  } else {
    initMat<-DBNinit(backDBN,scorepar$nsmall,scorepar$bgn)
    transMat<-backDBN
  }

  if(!is.null(datatoscore)) {
    datasplit<-splitDBNdata(datatoscore,scorepar)
  } else {
    datasplit<-NULL
  }
  if(scorepar$stationary) {
    if(scorepar$slices==2) {
    totscore<-scoreagainstDAG(scorepar$firstslice, initMat, datatoscore=datasplit$init,onlymain=TRUE)+
      scoreagainstDAG(scorepar$otherslices, transMat, datatoscore=datasplit$trans,onlymain=TRUE)
    } else {
      totscore<-0
      if(is.null(datainit)) {
        totscore<-totscore+scoreagainstDAG(scorepar$firstslice, initMat, datatoscore=datasplit[[scorepar$slices]],onlymain=TRUE)
      } else {
        for(i in 1:length(datainit)) {
          totscore<-totscore+scoreagainstDAG(scorepar$firstslice, initMat, datatoscore=datainit[[i]],onlymain=TRUE)
        }
        totscore<-totscore/length(datainit)
      }
        for(i in 1:(scorepar$slices-1)) {
          nas<-which(apply(datasplit[[i]],1,function(x)all(!is.na(x))))
          if((nrow(datasplit[[i]])-length(nas))>0) {
            addscore<-rep(0,nrow(datasplit[[i]]))
            addscore[nas]<-scoreagainstDAG(scorepar$otherslices, transMat, datatoscore=datasplit[[i]][nas,],onlymain=TRUE)
          } else {
            addscore<-scoreagainstDAG(scorepar$otherslices, transMat, datatoscore=datasplit[[i]],onlymain=TRUE)
          }
          totscore<-totscore+addscore
        }
    }
  } else {
    totscore<-0
    if(is.null(datainit)) {
      totscore<-totscore+scoreagainstDAG(scorepar$paramsets[[scorepar$slices]], initMat, datatoscore=datasplit[[scorepar$slices]],onlymain=TRUE)
    } else {
      for(i in 1:length(datainit)) {
        totscore<-totscore+scoreagainstDAG(scorepar$paramsets[[scorepar$slices]], initMat, datatoscore=datainit[[i]],onlymain=TRUE)
      }
      totscore<-totscore/length(datainit)
      for(i in 1:(scorepar$slices-1)) {
        nas<-which(apply(datasplit[[i]],1,function(x)all(!is.na(x))))
        if((nrow(datasplit[[i]])-length(nas))>0) {
          addscore<-rep(0,nrow(datasplit[[i]]))
          addscore[nas]<-scoreagainstDAG(scorepar$paramsets[[i]], transMat, datatoscore=datasplit[[i]][nas,],onlymain=TRUE)
        } else {
          addscore<-scoreagainstDAG(scorepar$paramsets[[i]], transMat, datatoscore=datasplit[[i]],onlymain=TRUE)
        }
        totscore<-totscore+addscore
      }
    }
  }
  return(totscore)
}

#'Score against DBN
#'
#'Scoring observations against a DBN structure
#'
#'@param scorepar object of class 'scoreparameters'
#'@param incidence adjacency matrix of a DAG
#'@param datatoscore matrix or vector containing observations to be scored
#'@param marginalise (logical) should marginal score be used?
#'@param onlymain (logical) should static nodes be included in the score?
#'@param datainit optional, in case of unbalanced design, the mean score of available samples for T0 are computed
#'@return vector of log-scores
#'@author Polina Suter
scoreagainstDBN3<-function(scorepar, incidence, datatoscore=NULL,
                          marginalise=FALSE, onlymain=FALSE, datainit=NULL){
  backDBN<-DBNbacktransform_l(incidence,scorepar,coln=TRUE)
  if(scorepar$split) {
    initMat<-backDBN$init
    transMat<-backDBN$trans
  } else {
    initMat<-DBNinit(backDBN,scorepar$nsmall,scorepar$bgn)
    transMat<-backDBN
  }

  if(!is.null(datatoscore)) {
    datasplit<-splitDBNdata(datatoscore,scorepar)
  } else {
    datasplit<-NULL
  }
  if(scorepar$stationary) {
    if(scorepar$slices==2) {
      totscore<-scoreagainstDAG(scorepar$firstslice, initMat, datatoscore=datasplit$init,onlymain=TRUE)+
        scoreagainstDAG(scorepar$otherslices, transMat, datatoscore=datasplit$trans,onlymain=TRUE)
    } else {
      totscore<-0
      if(is.null(datainit)) {
        totscore<-totscore+scoreagainstDAG(scorepar$firstslice, initMat, datatoscore=datasplit[[scorepar$slices]],onlymain=TRUE)
      } else {
        for(i in 1:length(datainit)) {
          totscore<-totscore+scoreagainstDAG(scorepar$firstslice, initMat, datatoscore=datainit[[i]],onlymain=TRUE)
        }
        totscore<-totscore/length(datainit)
      }
      for(i in 1:(scorepar$slices-1)) {
        nas<-which(apply(datasplit[[i]],1,function(x)all(!is.na(x))))
        if((nrow(datasplit[[i]])-length(nas))>0) {
          addscore<-rep(0,nrow(datasplit[[i]]))
          addscore[nas]<-scoreagainstDAG(scorepar$otherslices, transMat, datatoscore=datasplit[[i]][nas,],onlymain=TRUE)
        } else {
          addscore<-scoreagainstDAG(scorepar$otherslices, transMat, datatoscore=datasplit[[i]],onlymain=TRUE)
        }
        totscore<-addscore
      }
    }
  } else {
    totscore<-0
    if(is.null(datainit)) {
      totscore<-totscore+scoreagainstDAG(scorepar$paramsets[[scorepar$slices]], initMat, datatoscore=datasplit[[scorepar$slices]],onlymain=TRUE)
    } else {
      for(i in 1:length(datainit)) {
        totscore<-totscore+scoreagainstDAG(scorepar$paramsets[[scorepar$slices]], initMat, datatoscore=datainit[[i]],onlymain=TRUE)
      }
      totscore<-totscore/length(datainit)
      for(i in 1:(scorepar$slices-1)) {
        nas<-which(apply(datasplit[[i]],1,function(x)all(!is.na(x))))
        if((nrow(datasplit[[i]])-length(nas))>0) {
          addscore<-rep(0,nrow(datasplit[[i]]))
          addscore[nas]<-scoreagainstDAG(scorepar$paramsets[[i]], transMat, datatoscore=datasplit[[i]][nas,],onlymain=TRUE)
        } else {
          addscore<-scoreagainstDAG(scorepar$paramsets[[i]], transMat, datatoscore=datasplit[[i]],onlymain=TRUE)
        }
        totscore<-addscore
      }
    }
  }
  return(totscore)
}
splitDBNdata<-function(datatoscore,param,addinit=NULL) {
  datasplit<-list()
    if(param$slices==2) {
    if(param$bgn>0) {
      datasplit$init<-datatoscore[,c(1:param$nsmall+param$bgn,1:param$bgn)]
      datasplit$trans<-datatoscore[,c(1:param$nsmall+param$nsmall+param$bgn,1:param$bgn,1:param$nsmall+param$bgn)]
    } else {
      datasplit$init<-datatoscore[,1:param$nsmall]
      datasplit$trans<-datatoscore[,c(1:param$nsmall+param$nsmall,1:param$nsmall)]
    }
    } else {
      for (i in 1:(param$slices-1)) {
        datasplit[[i]]<-datatoscore[,c(1:param$nsmall+i*param$nsmall,1:param$nsmall+(i-1)*param$nsmall)]
      }
      if(is.null(addinit)) {
        datasplit[[param$slices]]<-datatoscore[,1:param$nsmall]
      } else {
        datasplit[[param$slices]]<-addinit
      }
    }
  return(datasplit)
}
DBNbacktransform_l<-function(DBNmat,param,coln=FALSE) {
  if(!is.null(colnames(DBNmat))) {
    oldnodelabels<-colnames(DBNmat)
    newnodelabels<-oldnodelabels
    newnodelabels[param$intstr$cols]<-oldnodelabels[param$usrtrans$cols]
    if(param$bgn==0) newnodelabels[param$trans$rows]<-oldnodelabels[param$usrinitstr$rows] else {
      newnodelabels[c(param$intstr$rows[1:param$bgn],param$trans$rows)]<-oldnodelabels[param$usrinitstr$rows]
    }
  }
  newDBNmat<-matrix(0,nrow=param$n+param$nsmall,ncol=param$n+param$nsmall)
  newDBNmat[param$intstr$rows,param$intstr$cols]<-1*(DBNmat[param$usrintstr$rows,param$usrintstr$cols]|DBNmat[param$usrinitstr$rows,param$usrinitstr$cols])
  newDBNmat[param$trans$rows,param$trans$cols]<-DBNmat[param$usrtrans$rows,param$usrtrans$cols]
  if(!param$split) {
    if(coln) colnames(newDBNmat)<-rownames(newDBNmat)<-newnodelabels
    return(newDBNmat)
  } else {
    res<-list()
    initDBNmat<-DBNmat[1:param$n,1:param$n]
    newinitDBNmat<-DBNmat[1:param$n,1:param$n]

    if(param$bgn>0) {
    newinitDBNmat[,1:param$bgn+param$nsmall]<-initDBNmat[,1:param$bgn]
    }
    newinitDBNmat[,1:param$nsmall]<-initDBNmat[,1:param$nsmall+param$bgn]
    initDBNmat<-newinitDBNmat
    if(param$bgn>0) {
    newinitDBNmat[1:param$bgn+param$nsmall,]<-initDBNmat[1:param$bgn,]
    }
    newinitDBNmat[1:param$nsmall,]<-initDBNmat[1:param$nsmall+param$bgn,]
    res$init<-newinitDBNmat

    transDBNmat<-matrix(0,nrow=param$n+param$nsmall,ncol=param$n+param$nsmall)
    DBNmat<-DBNcut_l(DBNmat,dyn=param$nsmall,b=param$bgn)
    transDBNmat[param$intstr$rows,param$intstr$cols]<-DBNmat[param$usrintstr$rows,param$usrintstr$cols]
    transDBNmat[param$trans$rows,param$trans$cols]<-DBNmat[param$usrtrans$rows,param$usrtrans$cols]
    res$trans<-transDBNmat

    return(res)
  }
}
DBNcut_l<-function(adj,dyn,b){
  adj[,1:(dyn+b)]<-0
  return(adj)
}
DBNunitedata<-function(dbndata,dyn,b){
nsmall<-dyn
n<-b+dyn
bgn<-b
slices<-(ncol(dbndata)-bgn)/dyn
    if(all(is.character(colnames(dbndata)))){
      nodeslabels<-colnames(dbndata)
    } else {
      if(b>0) {
        staticnames<-sapply(c(1:bgn), function(x)paste("s",x,sep=""))
        dynamicnames<-rep(sapply(c(1:nsmall), function(x)paste("v",x,sep="")),slices)
        for(i in 2:slices) {
          dynamicnames[1:nsmall+(i-1)*nsmall]<-paste(dynamicnames[1:nsmall+(i-1)*nsmall],".",i,sep="")
        }
        nodeslabels<-c(staticnames,dynamicnames)
      } else {
        nodeslabels<-rep(sapply(c(1:n), function(x)paste("v",x,sep="")),slices)
        for(i in 2:slices) {
          nodeslabels[1:nsmall+(i-1)*nsmall]<-paste(nodeslabels[1:nsmall+(i-1)*nsmall],".",i,sep="")
        }
      }
      colnames(dbndata)<-nodeslabels
    }
   labels.short<-nodeslabels[1:(n+nsmall)]
   # other slices we layer the data,
   datalocal <- dbndata[,1:(2*nsmall+bgn)]
   collabels<-colnames(datalocal)
   if (bgn>0){
     bgdata<-dbndata[,1:b]
     if(slices > 2){ # layer on later time slices
       for(jj in 1:(slices-2)){
         datatobind<-cbind(bgdata,dbndata[,nsmall*jj+1:(2*nsmall)+bgn])
         colnames(datatobind)<-collabels
         datalocal <- rbind(datalocal,datatobind)
       }
     }
   } else {
     if(slices > 2){ # layer on later time slices
       for(jj in 1:(slices-2)){
         datatobind<-dbndata[,n*jj+1:(2*n)]
         colnames(datatobind)<-collabels
         datalocal <- rbind(datalocal,datatobind)
       }
     }
   }
return(datalocal)
}
