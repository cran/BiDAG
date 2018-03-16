newspaceskel<-function(n,startspace,currspace,softlimit,hardlimit,posterior,blacklist,MCMCchain=NULL,mergetype="skeleton") {
  switch(mergetype,
         "dag" = { 
           mdag<-dag.threshold(n,MCMCchain,pbarrier=posterior,pdag=FALSE)
           newadj<-1*(!blacklist&(startspace|mdag))
           toomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(toomanyneib)>0){newadj[,toomanyneib]<-(1*(!blacklist&mdag))[,toomanyneib]}
           toomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(toomanyneib)>0){newadj[,toomanyneib]<-currspace[,toomanyneib]}
         },
         "cpdag" = { 
           mcp<-dag.threshold(n,MCMCchain,pbarrier=posterior,pdag=TRUE)
           newadj<-1*(!blacklist&(startspace|mcp))
           toomanyneib<-which(apply(newadj,2,sum)>softlimit)
           if(length(toomanyneib)>0) {
             mdag<-dag.threshold(n,MCMCchain,pbarrier=posterior,pdag=FALSE)
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|mdag)))[,toomanyneib]
           }
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-mdag[,tootoomanyneib]}
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-currspace[,tootoomanyneib]}
         },
         "skeleton" = { 
           mskel<-1*(dag.threshold(n,MCMCchain,pbarrier=posterior,pdag=FALSE)|t(dag.threshold(n,MCMCchain,pbarrier=posterior,pdag=FALSE)))
           newadj<-1*(!blacklist&(startspace|mskel))
           toomanyneib<-which(apply(newadj,2,sum)>4)
           if(length(toomanyneib)>0) {
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|dag.threshold(n,MCMCchain,pbarrier=posterior,pdag=TRUE))))[,toomanyneib]
           }
           toomanyneib<-which(apply(newadj,2,sum)>softlimit)
           if(length(toomanyneib)>0) {
             mdag<-dag.threshold(n,MCMCchain,pbarrier=posterior,pdag=FALSE)
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|mdag)))[,toomanyneib]
           }
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-mdag[,tootoomanyneib]}
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-currspace[,tootoomanyneib]}
         }
  )
  return(newadj)
}

newspacemap<-function(n,startspace,currspace,softlimit,hardlimit,blacklist, maxN,MCMCchain=NULL,mergetype="skeleton") {
  switch(mergetype,
         "dag" = { 
           maxdag<-MCMCchain[[maxN]]
           newadj<-1*(!blacklist&(startspace|maxdag))
           toomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(toomanyneib)>0){newadj[,toomanyneib]<-(1*(!blacklist&maxdag))[,toomanyneib]}
           toomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(toomanyneib)>0){newadj[,toomanyneib]<-currspace[,toomanyneib]}
         },
         "cpdag" = { 
           maxcp<-dagadj2cpadj(MCMCchain[[maxN]])
           newadj<-1*(!blacklist&(startspace|maxcp))
           toomanyneib<-which(apply(newadj,2,sum)>softlimit)
           if(length(toomanyneib)>0) {
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|MCMCchain[[maxN]])))[,toomanyneib]
           }
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-MCMCchain[[maxN]][,tootoomanyneib]}
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-currspace[,tootoomanyneib]}
         },
         "skeleton" = { 
           maxskel<-1*(MCMCchain[[maxN]]|t(MCMCchain[[maxN]]))
           newadj<-1*(!blacklist&(startspace|maxskel))
           toomanyneib<-which(apply(newadj,2,sum)>4)
           if(length(toomanyneib)>0) {
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|dagadj2cpadj(MCMCchain[[maxN]]))))[,toomanyneib]
           }
           toomanyneib<-which(apply(newadj,2,sum)>softlimit)
           if(length(toomanyneib)>0) {
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|MCMCchain[[maxN]])))[,toomanyneib]
           }
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-MCMCchain[[maxN]][,tootoomanyneib]}
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-currspace[,tootoomanyneib]}
         }
  )
  return(newadj)
}

definestartspace<-function(alpha,param,cpdag=FALSE,n,algo="pc") {
  if(is.null(alpha)) { if(n<50) {alpha<-0.4} else {alpha<-max(20/n,0.01)}}
  
  if(param$type=="bde") {
    if(cpdag){
      pc.skel<-pc(suffStat = list(d1=param$d1,d0=param$d0,data=param$data),
                  indepTest = weightedbinCItest, alpha = alpha, labels = colnames(param$data),
                  verbose = FALSE)
      
    } else {
      pc.skel<-pcalg::skeleton(suffStat = list(d1=param$d1,d0=param$d0,data=param$data),
                               indepTest = weightedbinCItest, alpha = alpha, labels = colnames(param$data),
                               verbose = FALSE)
    }
    
  } else if(param$type=="bge") {
    if(is.null(param$weightvector)) {
      cormat<-cor(param$data)
      N<-nrow(param$data)
    } else { N<-sum(param$weightvector)
    cormat<-cov.wt(param$data,wt=param$weightvector,cor=TRUE)$cor}
    if(cpdag){
      pc.skel<-pcalg::pc(suffStat = list(C = cormat, n = N),
                         indepTest = gaussCItest,
                         alpha=alpha,labels=colnames(param$data),skel.method="stable",verbose = FALSE)
    } else {
      pc.skel<-pcalg::skeleton(suffStat = list(C = cormat, n = N),
                               indepTest = gaussCItest,
                               alpha=alpha,labels=colnames(param$data),method="stable",verbose = FALSE)
    }
  }
  
  g<-pc.skel@graph
  startspace<-1*(dag2adjacencymatrix(g))
}



