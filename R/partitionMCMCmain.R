#implements partition MCMC scheme for structure learning problem, searches in plus1 neighbourhood of a search space defined
#by startspace or iterativeMCMC function
#authors Polina Suter, Jack Kuipers, partly derived from <doi:10.1080/01621459.2015.1133426>
partitionMCMCplus1sample<-function(param,startspace,blacklist=NULL,moveprobs,numit,stepsave,
                                   startorder=NULL,scoretable=NULL,DAG,gamma=1,verbose=TRUE,
                                   scoreout=FALSE){

  MCMCtraces<-list()
  
  n<-param$n
  nsmall<-param$nsmall
  matsize<-ifelse(param$DBN,n+nsmall,n)
  
  if(!param$DBN) {
    if(param$bgn!=0) {
      updatenodes<-c(1:n)[-param$bgnodes]
    } else {
      updatenodes<-c(1:n)
    }
  } else {
    updatenodes<-c(1:nsmall)
  }
  
  if(is.null(startspace) & is.null(scoretable)) {
    if(verbose) cat("defining a search space with iterativeMCMC \n")
    searchspace<-iterativeMCMC(scorepar=param,moveprobs=NULL,plus1it=NULL,
                                     iterations=NULL,stepsave=NULL,softlimit=9,hardlimit=14,alpha=NULL,
                                     verbose=verbose,chainout=FALSE,scoreout=TRUE,
                                     gamma=gamma,cpdag=FALSE,mergetype="skeleton",
                                     blacklist=blacklist)
    startspace<-searchspace$scoretable$adjacency
    if(param$DBN) {
      maxDAG<-DBNbacktransform(searchspace$DAG,param)
    } else {
      maxDAG<-searchspace$DAG
      scoretable<-searchspace$scoretable$tables
    }


    if(!is.null(param$bgnodes)) {
      forpart<-DAGtopartition(param$nsmall,maxDAG[updatenodes,updatenodes])
    } else {
      forpart<-DAGtopartition(n,maxDAG)
    }
    
    if (is.null(blacklist)) {
      blacklist<-matrix(0,nrow=matsize,ncol=matsize)
    }
    
    diag(blacklist)<-1
    if(!is.null(param$bgnodes)) {
      for(i in param$bgnodes) {
        blacklist[,i]<-1
      }
    }
    blacklistparents<-list()
    for  (i in 1:matsize) {
      blacklistparents[[i]]<-which(blacklist[,i]==1)
    }
    
  } else {


      if (is.null(blacklist)) {
        blacklist<-matrix(0,nrow=matsize,ncol=matsize)
      }
      
      diag(blacklist)<-1
      if(!is.null(param$bgnodes)) {
        for(i in param$bgnodes) {
          blacklist[,i]<-1
        }
      }
     
    if(is.null(scoretable)) { 
      startspace<-1*(startspace&!blacklist)
    } else {
      startspace<-scoretable$adjacency
      blacklist<-scoretable$blacklist
      scoretable<-scoretable$tables
    }
    
      blacklistparents<-list()
      for  (i in 1:matsize) {
        blacklistparents[[i]]<-which(blacklist[,i]==1)
      }
      
    
    if (is.null(DAG)) {
      forpart<-list()
      if(is.null(param$bgnodes)) {
        forpart$permy<-c(1:n)
        forpart$party<-c(n)
        updatenodes<-c(1:n)
      } else {
        forpart$permy<-c(1:n)[-param$bgnodes]
        forpart$party<-c(param$nsmall)
        updatenodes<-c(1:n)[-param$bgnodes]
      }
    } else {
      
      if(!is.null(param$bgnodes)) {
        forpart<-DAGtopartition(param$nsmall,DAG[updatenodes,updatenodes])
      } else {
        forpart<-DAGtopartition(n,DAG)
      }
      }
  }
  
  
  if(verbose) cat("search space identified, score tables to be calculated \n")
  permy<-forpart$permy
  party<-forpart$party
  starttable<-Sys.time()
  ptab<-listpossibleparents.PC.aliases(startspace,isgraphNEL=FALSE,n,updatenodes=updatenodes)
  parenttable<-ptab$parenttable
  aliases<-ptab$aliases
  numberofparentsvec<-ptab$numberofparentsvec
  numparents<-ptab$numparents
  plus1lists<-PLUS1(n,aliases,updatenodes,blacklistparents)
  rowmapsallowed<-parentsmapping(parenttable,numberofparentsvec,n,updatenodes)
  if (is.null(scoretable)) {
    scoretable<-scorepossibleparents.PLUS1(parenttable,plus1lists,n,param,updatenodes,
                                           rowmapsallowed,numparents,numberofparentsvec)
    }
  scoretab<-list()
  for (i in updatenodes) {
  scoretab[[i]]<-matrix(sapply(scoretable[[i]], unlist),nrow=nrow(scoretable[[i]][[1]]))}
  posetparenttable<-poset(ptab$parenttable,ptab$numberofparentsvec,rowmapsallowed,n,updatenodes=updatenodes)
  plus1allowedpart<-plus1allowed.partition(posetparenttable,scoretable,numberofparentsvec,rowmapsallowed,
                                n,plus1lists=plus1lists,numparents,updatenodes)
  needednodetable<-partitionlist(parenttable,ptab$numberofparentsvec,n,updatenodes=updatenodes)
  numberofpartitionparentsvec<-partitionlistnumberofparents(needednodetable,ptab$numberofparentsvec,
                                                            n,updatenodes)
  needednodebannedrow<-partitionmapneedednodebannedrow(ptab$numparents,ptab$numberofparentsvec,n,
                                                       updatenodes)
  rowmapsneeded<-neededparentsmapping(parenttable,ptab$numberofparentsvec,needednodetable,
                                       numberofpartitionparentsvec,needednodebannedrow,n,updatenodes)
  neededposetparents<-needed.poset(ptab$parenttable,ptab$numberofparentsvec,
                                   needednodebannedrow,rowmapsneeded,n,updatenodes)
  neededposetparenttable<-lapply(neededposetparents,function(x)x$table)
  neededposetparentsvec<-lapply(neededposetparents,function(x)x$sizes)
  plus1neededpart<-plus1needed.partition(ptab$numparents,parenttable,neededposetparenttable,neededposetparentsvec,
                                           ptab$numberofparentsvec,rowmapsallowed,needednodebannedrow,scoretable,
                                           plus1lists,n,updatenodes)
  endtable<-Sys.time()
  if(verbose) cat("score tables calculated, partition MCMC starts \n")
    partres<-partitionMCMCplus1(n,param$nsmall,permy,party,numit,stepsave,parenttable,scoretable,scoretab,
                                 aliases,plus1neededpart,plus1allowedpart,plus1lists,rowmapsneeded,rowmapsallowed,
                                 needednodetable,ptab$numberofparentsvec,
                                 numberofpartitionparentsvec,needednodebannedrow,
                                 neededposetparentsvec,moveprobs,param$bgnodes,matsize=matsize)
  endmcmc<-Sys.time()

  
  
    if(param$DBN) {
      MCMCchain<-lapply(partres$incidence,function(x)DBNtransform(x,param=param))
      MCMCtraces$incidence<-MCMCchain
    } else {
      MCMCtraces$incidence<-lapply(partres$incidence,function(x)assignLabels(x,param$labels))
    }
    MCMCtraces$DAGscores<-partres$DAGscores
    MCMCtraces$partitionscores<-partres$partitionscores
    MCMCtraces$order<-partres$order
    MCMCtraces$partition<-partres$partition
    
   
  
  maxobj<-storemaxMCMC(partres,param)
  maxN<-which.max(unlist(partres[[2]]))

  
  
  result<-list()
  result$info<-list()
  tabletime<-endtable-starttable
  if(units(tabletime)=="mins") {
    tabletime<-as.numeric(tabletime*60)
  }
  mcmctime<-endmcmc-endtable
  if(units(mcmctime)=="mins") {
    mcmctime<-as.numeric(mcmctime*60)
  }
  result$info$times<-c(tabletime,mcmctime)
  result$trace<-MCMCtraces$DAGscores
  MCMCtraces$DAGscores<-NULL
  result$traceadd<-MCMCtraces
  result$DAG<-maxobj$DAG
  result$CPDAG<-graph2m(dag2cpdag(m2graph(result$DAG)))
  result$score<-maxobj$score
  result$maxorder<-maxobj$order
  if(scoreout) {
    result$scoretable<-list()
    result$scoretable$adjacency<-startspace
    result$scoretable$tables<-scoretable
    result$scoretable$blacklist<-blacklist
    attr( result$scoretable,"class")<-"MCMCscoretab"
  }
  attr(result,"class")<-"MCMCres"
  
  return(result)
}
