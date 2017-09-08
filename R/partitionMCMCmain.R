#implements partition MCMC scheme for structure learning problem, searches in plus1 neighbourhood of a search space defined
#by startspace or iterativeMCMCsearch function

partitionMCMCplus1sample<-function(n,param,startspace,blacklist=NULL,moveprobs,numit,stepsave,
                                   startorder=NULL,scoretable=NULL,DAG,gamma=1,verbose=TRUE){
  
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
                                     gamma=gamma,cpdag=FALSE,mergetype="skeleton",blacklist=blacklist)
    
    maxDAG<-searchspace$max$DAG
    startspace<-searchspace$space$adjacency
    scoretable<-searchspace$space$scoretable
    forpart<-DAGtopartition(n,maxDAG)
  } else {
    startspace<-1*(startspace&!blacklist)
    if (is.null(DAG)) {
      forpart<-list()
      forpart$permy<-c(1:n)
      forpart$party=c(n)
    } else {forpart<-DAGtopartition(n,DAG)}
    }
  if(verbose){print("search space identified, score tables to be calculated")}
  permy<-forpart$permy
  party<-forpart$party
  starttime<-Sys.time()
  ptab<-listpossibleparents.PC.aliases(startspace,isgraphNEL=FALSE,n)
  parenttable<-ptab$parenttable
  aliases<-ptab$aliases
  numberofparentsvec<-ptab$numberofparentsvec
  numparents<-ptab$numparents
  updatenodes<-c(1:n)
  plus1lists<-PLUS1(n,aliases,updatenodes,blacklistparents)
  rowmapsneeded<-parentsmapping(parenttable,numberofparentsvec,n,updatenodes)
  if (is.null(scoretable)) {
    scoretable<-scorepossibleparents.PLUS1(parenttable,plus1lists,n,param,updatenodes,rowmapsneeded,numparents,numberofparentsvec)
    }
  scoretab<-list()
  for (i in 1:n) {
  scoretab[[i]]<-matrix(sapply(scoretable[[i]], unlist),nrow=nrow(scoretable[[i]][[1]]))}
  posetparenttable<-poset(ptab$parenttable,ptab$numberofparentsvec,rowmapsneeded,n)
  plus1neededpart<-poset.scores(posetparenttable,scoretable,numberofparentsvec,rowmapsneeded,
                                n,plus1lists=plus1lists,numparents)
  needednodetable<-partitionlist(parenttable,ptab$numberofparentsvec,n)
  numberofpartitionparentsvec<-partitionlistnumberofparents(needednodetable,ptab$numberofparentsvec,n)
  needednodebannedrow<-partitionmapneedednodebannedrow(ptab$numparents,ptab$numberofparentsvec,n)
  rowmapsallowed<-neededparentsmapping(parenttable,ptab$numberofparentsvec,needednodetable,
                                       numberofpartitionparentsvec,needednodebannedrow,n)
  neededposetparents<-needed.poset(ptab$parenttable,ptab$numberofparentsvec,
                                   needednodebannedrow,rowmapsallowed,n)
  neededposetparenttable<-lapply(neededposetparents,function(x)x$table)
  neededposetparentsvec<-lapply(neededposetparents,function(x)x$sizes)
  plus1allowedpart<-plus1allowed.partition(ptab$numparents,parenttable,neededposetparenttable,neededposetparentsvec,
                                           ptab$numberofparentsvec,rowmapsneeded,needednodebannedrow,scoretable,
                                           plus1lists,n)
  if(verbose){print("score tables calculated, partition MCMC starts")}
  partresult<-partitionMCMCplus1(n,permy,party,numit,stepsave,parenttable,scoretable,scoretab,
                                 aliases,plus1neededpart,plus1allowedpart,plus1lists,rowmapsneeded,rowmapsallowed,
                                 needednodetable,ptab$numberofparentsvec,
                                 numberofpartitionparentsvec,needednodebannedrow,neededposetparentsvec,moveprobs)
  endtime<-Sys.time()
  if(verbose){print(endtime-starttime)}
  result<-list()
  result$chain<-partresult
  attr(result$chain,"class")<-"MCMCchain"
  return(result)
}
