#implements partition MCMC scheme for structure learning problem, searches in plus1 neighbourhood of a search space defined
#by startspace or iterativeMCMCsearch function

partitionMCMCplus1sample<-function(n,param,startspace,blacklist=NULL,moveprobs,numit,stepsave,
                                   startorder=NULL,scoretable=NULL,DAG,gamma=1,verbose=TRUE){

  if(!is.null(param$bgnodes)) {
    mainnodes<-c(1:n)[-param$bgnodes]
    bglogscore<-sum(bgNodeScore(n,param$bgnodes,param)[param$bgnodes])
    startorder<-orderbgn(startorder,param$bgnodes)
  } else {
    mainnodes<-c(1:n)
    bglogscore<-0
  }
  
  if (is.null(blacklist)) {
    blacklist<-matrix(rep(0,n*n),ncol=n)
  } 
  diag(blacklist)<-1
  blacklistparents<-list()
  for  (i in 1:n) {
    blacklistparents[[i]]<-which(blacklist[,i]==1)
  }
  
  if(is.null(startspace)) {
    if(verbose){print("defining a search space with iterativeMCMCsearch")}
    searchspace<-iterativeMCMCsearch(n,scoreparam=param,moveprobs=NULL,plus1it=NULL,
                                     iterations=NULL,stepsave=NULL,softlimit=9,hardlimit=14,alpha=NULL,
                                     verbose=verbose,chainout=FALSE,scoreout=TRUE,
                                     gamma=gamma,cpdag=FALSE,mergetype="skeleton",
                                     blacklist=blacklist)
    
    maxDAG<-searchspace$max$DAG
    startspace<-searchspace$space$adjacency
    scoretable<-searchspace$space$scoretable
    forpart<-DAGtopartition(n,maxDAG,param$bgnodes)
  } else {
    startspace<-1*(startspace&!blacklist)
    if (is.null(DAG)) {
      forpart<-list()
      if(is.null(param$bgnodes)) {
      forpart$permy<-c(1:n)
      forpart$party=c(n)
      mainnodes<-c(1:n)
      } else {
        forpart$permy<-c(1:n)[-param$bgnodes]
        forpart$party=c(param$nsmall)
        mainnodes<-c(1:n)[-param$bgnodes]
      }
    } else {
      
      if(!is.null(param$bgnodes)) {
        forpart<-DAGtopartition(param$nsmall,DAG[mainnodes,mainnodes])
      } else {
        forpart<-DAGtopartition(n,DAG)
      }
      }
    }
  if(verbose){print("search space identified, score tables to be calculated")}
  permy<-forpart$permy
  party<-forpart$party
  starttime<-Sys.time()
  ptab<-listpossibleparents.PC.aliases(startspace,isgraphNEL=FALSE,n,updatenodes=mainnodes)
  parenttable<-ptab$parenttable
  aliases<-ptab$aliases
  numberofparentsvec<-ptab$numberofparentsvec
  numparents<-ptab$numparents
  plus1lists<-PLUS1(n,aliases,mainnodes,blacklistparents)
  rowmapsallowed<-parentsmapping(parenttable,numberofparentsvec,n,mainnodes)
  if (is.null(scoretable)) {
    scoretable<-scorepossibleparents.PLUS1(parenttable,plus1lists,n,param,mainnodes,rowmapsneeded,numparents,numberofparentsvec)
    }
  scoretab<-list()
  for (i in mainnodes) {
  scoretab[[i]]<-matrix(sapply(scoretable[[i]], unlist),nrow=nrow(scoretable[[i]][[1]]))}
  posetparenttable<-poset(ptab$parenttable,ptab$numberofparentsvec,rowmapsallowed,n,updatenodes=mainnodes)
  plus1allowedpart<-plus1allowed.partition(posetparenttable,scoretable,numberofparentsvec,rowmapsallowed,
                                n,plus1lists=plus1lists,numparents,mainnodes)
  needednodetable<-partitionlist(parenttable,ptab$numberofparentsvec,n,updatenodes=mainnodes)
  numberofpartitionparentsvec<-partitionlistnumberofparents(needednodetable,ptab$numberofparentsvec,
                                                            n,mainnodes)
  needednodebannedrow<-partitionmapneedednodebannedrow(ptab$numparents,ptab$numberofparentsvec,n,
                                                       mainnodes)
  rowmapsneeded<-neededparentsmapping(parenttable,ptab$numberofparentsvec,needednodetable,
                                       numberofpartitionparentsvec,needednodebannedrow,n,mainnodes)
  neededposetparents<-needed.poset(ptab$parenttable,ptab$numberofparentsvec,
                                   needednodebannedrow,rowmapsneeded,n,mainnodes)
  neededposetparenttable<-lapply(neededposetparents,function(x)x$table)
  neededposetparentsvec<-lapply(neededposetparents,function(x)x$sizes)
  plus1neededpart<-plus1needed.partition(ptab$numparents,parenttable,neededposetparenttable,neededposetparentsvec,
                                           ptab$numberofparentsvec,rowmapsallowed,needednodebannedrow,scoretable,
                                           plus1lists,n,mainnodes)
  if(verbose){print("score tables calculated, partition MCMC starts")}
    partresult<-partitionMCMCplus1(n,param$nsmall,permy,party,numit,stepsave,parenttable,scoretable,scoretab,
                                 aliases,plus1neededpart,plus1allowedpart,plus1lists,rowmapsneeded,rowmapsallowed,
                                 needednodetable,ptab$numberofparentsvec,
                                 numberofpartitionparentsvec,needednodebannedrow,
                                 neededposetparentsvec,moveprobs,param$bgnodes,bglogscore)
  endtime<-Sys.time()
  if(verbose){print(endtime-starttime)}
  result<-list()
  result$chain<-partresult
  attr(result$chain,"class")<-"MCMCtracepartr"
  return(result)
}
