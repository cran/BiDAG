#iterative MCMC search
iterativeMCMCplus1<-function(n,param,iterations,stepsave,plus1it=NULL,MAP=TRUE, posterior=0.5,
                             startorder=c(1:n),moveprobs,softlimit=9,hardlimit=12,plus1=FALSE,chainout=FALSE,
                             scoreout=FALSE,startspace=NULL,blacklist=NULL,gamma=1,verbose=FALSE,alpha=NULL,
                             cpdag=FALSE,mergecp="skeleton",addspace=NULL,scoretable=NULL) {
    updatenodeslist<-list()
    MCMCchains<-list()
    updatenodes<-c(1:n)
    maxorder<-startorder
    
    if (is.null(blacklist)) {
      blacklist<-matrix(rep(0,n*n),ncol=n)
    }
    diag(blacklist)<-1
    blacklistparents<-list()
    for  (i in 1:n) {
      blacklistparents[[i]]<-which(blacklist[,i]==1)
    }
    
    if (is.null(startspace)) startspace<-definestartspace(alpha,param,cpdag,n,algo="pc")

    startskeleton<-1*(startspace&!blacklist)
    
    if(!is.null(addspace)) {startskel<-1*((addspace|startskeleton)&!blacklist)
    } else { startskel<-startskeleton }
    
    ptab<-listpossibleparents.PC.aliases(startskel,isgraphNEL=FALSE,n,updatenodes)

    if (verbose) {
      print("skeleton ready")
      flush.console()
      }
    if (plus1==FALSE){
      starttime<-Sys.time()
      parenttable<-ptab$parenttable # basic parenttable without plus1 lists
      aliases<-ptab$aliases #aliases for each node since all nodes in parent tables are named as 1,2,3,4. not real parent names
      numberofparentsvec<-ptab$numberofparentsvec
      numparents<-ptab$numparents
      rowmaps<-parentsmapping(parenttable,numberofparentsvec,n,updatenodes)
      scoretable<-scorepossibleparents.alias(parenttable,aliases,n,param,updatenodes,rowmaps,numparents,numberofparentsvec)
      posetparenttable<-poset(parenttable,numberofparentsvec,rowmaps,n,updatenodes)

      if(MAP==TRUE){
        maxmatrices<-posetscoremax(posetparenttable,scoretable,numberofparentsvec,rowmaps,
                                   n,plus1lists=NULL,updatenodes)

        MCMCresult<-orderMCMCbasemax(n,startorder,iterations,stepsave,moveprobs,parenttable,
                                     scoretable,aliases,numparents,rowmaps,maxmatrices,numberofparentsvec,gamma=gamma)
        endtime<-Sys.time()
        if(verbose) {
          print(endtime-starttime)
          flush.console()
          }
      } else {
        bannedscore<-poset.scores(posetparenttable,scoretable,numberofparentsvec,rowmaps,
                                  n,plus1lists=NULL,ptab$numparents)

        MCMCresult<-orderMCMCbase(n,startorder,iterations,stepsave,moveprobs,parenttable,
                                  scoretable,aliases,numparents,rowmaps,
                                  bannedscore,numberofparentsvec,gamma=gamma)
        endtime<-Sys.time()
        if(verbose) {
          print(endtime-starttime)
          flush.console()
          }
      }

      maxres<-list()
      maxN<-which.max(unlist(MCMCresult[[2]]))
      maxres$DAG<-MCMCresult[[1]][[maxN]]
      maxres$order<-MCMCresult[[4]][[maxN]]
      maxres$score<-MCMCresult[[2]][[maxN]]

      if (scoreout){
        if(chainout){output<-4}
        else{output<-3}
      } else{
        if(chainout) {output<-2}
        else {output<-1}
      }
      colnames(maxres$DAG)<-colnames(param$data)

      switch(as.character(output),
             "1"={ # return only maximum DAG and order
               result<-list()
               result$max<-maxres
               attr(result$max,"class")<-"MCMCmax"
               return(result)
             },
             "2"={ # return max DAG, order, last search space incidence and all scoretables
               result<-list()
               result$max<-maxres
               attr(result$max,"class")<-"MCMCmax"
               result$chain<-MCMCresult
               return(result)
             },
             "3"={ # return all MCMC all saved MCMC steps: incidence, DAGscore, orderscore and order and max result
               result<-list()
               result$max<-maxres
               attr(result$max,"class")<-"MCMCmax"
               result$space<-list()
               result$space$adjacency<-startskeleton
               result$space$scoretable<-scoretable
               return(result)
             },
             "4"={ # return all MCMC all saved MCMC steps,max result,last search space and scoretables
               result<-list()
               result$max<-maxres
               attr(result$max,"class")<-"MCMCmax"
               result$space<-list()
               result$space$adjacency<-startskeleton
               result$space$scoretable<-scoretable
               result$chain<-MCMCresult
               return(result)
             }
      )

    } else {
      starttime<-Sys.time()
      parenttable<-ptab$parenttable # basic parenttable without plus1 lists
      aliases<-ptab$aliases #aliases for each node since all nodes in parent tables are done as 1,2,3,4... not real parent names
      numberofparentsvec<-ptab$numberofparentsvec
      numparents<-ptab$numparents
      plus1lists<-PLUS1(n,aliases,updatenodes,blacklistparents)
      rowmaps<-parentsmapping(parenttable,numberofparentsvec,n,updatenodes)
      if(is.null(scoretable)) {
      scoretable<-scorepossibleparents.PLUS1(parenttable,plus1lists,n,param,updatenodes,rowmaps,numparents,numberofparentsvec) }
      posetparenttable<-poset(parenttable,numberofparentsvec,rowmaps,n,updatenodes)
      if(MAP==TRUE){
        maxmatrices<-posetscoremax(posetparenttable,scoretable,numberofparentsvec,
                                   rowmaps,n,plus1lists=plus1lists,updatenodes)
      } else {
        bannedscore<-poset.scores(posetparenttable,scoretable,ptab$numberofparentsvec,rowmaps,
                                  n,plus1lists=plus1lists,ptab$numparents,updatenodes)
      }
      oldadj<-startskeleton
      if(!is.null(plus1it)){
        for (i in 1:plus1it){#the number of iterations to run MCMC chain to expand the search space
          if(i>1){
            newptab<-listpossibleparents.PC.aliases(newadj,isgraphNEL=FALSE,n,updatenodes)
            parenttable[updatenodes]<-newptab$parenttable[updatenodes] # basic parenttable without plus1 lists
            aliases[updatenodes]<-newptab$aliases[updatenodes] #aliases for each node since all nodes in parent tables are done as 1,2,3,4... not real parent names
            numberofparentsvec[updatenodes]<-newptab$numberofparentsvec[updatenodes]
            numparents[updatenodes]<-newptab$numparents[updatenodes]
            newplus1lists<-PLUS1(n,aliases,updatenodes,blacklistparents)
            plus1lists$mask[updatenodes]<- newplus1lists$mask[updatenodes]
            plus1lists$parents[updatenodes]<- newplus1lists$parents[updatenodes]
            plus1lists$aliases[updatenodes]<- newplus1lists$aliases[updatenodes]

            rowmaps[updatenodes]<-parentsmapping(parenttable,numberofparentsvec,n,updatenodes)[updatenodes]
            scoretable[updatenodes]<-scorepossibleparents.PLUS1(parenttable,plus1lists,n,param,updatenodes,rowmaps,numparents,numberofparentsvec)[updatenodes]
            posetparenttable[updatenodes]<-poset(parenttable,numberofparentsvec,rowmaps,n,updatenodes)[updatenodes]

            if (MAP) {
              newmaxmatrices<-posetscoremax(posetparenttable,scoretable,numberofparentsvec,
                                            rowmaps,n,plus1lists=plus1lists,updatenodes)
              maxmatrices$maxmatrix[updatenodes]<- newmaxmatrices$maxmatrix[updatenodes]
              maxmatrices$maxrow[updatenodes]<- newmaxmatrices$maxrow[updatenodes]
            } else {
              newbannedscore<-poset.scores(posetparenttable,scoretable,numberofparentsvec,rowmaps,
                                           n,plus1lists=plus1lists,numparents,updatenodes)
              bannedscore[updatenodes]<-newbannedscore[updatenodes]
            }
            if(verbose) {
              print(paste("MCMC plus1 iteration",i))
              flush.console()
              }
          } else {
            if(verbose) {
              print(paste("score tables calculated, MCMC plus1 starts"))
              flush.console()
              }
            }
          if(MAP) {
            MCMCresult<-orderMCMCplus1max(n,startorder,iterations,stepsave,moveprobs,parenttable,
                                          scoretable,aliases,numparents,rowmaps,plus1lists,maxmatrices,numberofparentsvec,
                                          gamma=gamma)
            endtime<-Sys.time()
            if(verbose) {
              print(endtime-starttime)
              flush.console()}
          } else {
            starttime<-Sys.time()
            MCMCresult<-orderMCMCplus1(n,startorder,iterations,stepsave,moveprobs,parenttable,
                                       scoretable,aliases,numparents,rowmaps,plus1lists,
                                       bannedscore,numberofparentsvec,gamma=gamma)
            endtime<-Sys.time()
            print(endtime-starttime)
            if(verbose) {
              print(endtime-starttime)
              flush.console()}
          }

          MCMCchains[[i]]<-MCMCresult
          maxN<-which.max(unlist(MCMCresult[[2]]))
          if(i>1){
            if (MCMCresult[[2]][[maxN]]>maxscore){
              maxDAG<-MCMCresult[[1]][[maxN]]
              maxorder<-MCMCresult[[4]][[maxN]]
              maxscore<-MCMCresult[[2]][[maxN]]
              maxit<-i
            }
          } else {
            maxDAG<-MCMCresult[[1]][[maxN]]
            maxorder<-MCMCresult[[4]][[maxN]]
            maxscore<-MCMCresult[[2]][[maxN]]
            maxit<-1
          }
         
          if (i<plus1it) {
            if (MAP) {
              maxcp<-dagadj2cpadj(MCMCresult[[1]][[maxN]]) #find CPDAG
              if (plus1it>1){ #plus1it parameter defining the number of plus1 iterations
                newadj<-newspacemap(n,startskeleton,oldadj,softlimit,hardlimit,blacklist,maxN=maxN,MCMCchain=MCMCresult[[1]],mergetype=mergecp)
              } } else {
                newadj<- newspaceskel(n,startskeleton,oldadj,softlimit,hardlimit,posterior, blacklist,MCMCchain=MCMCresult[[1]],mergetype=mergecp)
              }
            updatenodes<-which(apply(newadj==oldadj,2,all)==FALSE)
            updatenodeslist[[i]]<-updatenodes
            oldadj<-newadj
          }
          startorder<-MCMCresult[[4]][[maxN]]
        }
      } else {
        i<-1
        while (length(updatenodes)>0){
          if(i>1){
            newptab<-listpossibleparents.PC.aliases(newadj,isgraphNEL=FALSE,n,updatenodes)
            parenttable[updatenodes]<-newptab$parenttable[updatenodes] # basic parenttable without plus1 lists
            aliases[updatenodes]<-newptab$aliases[updatenodes] #aliases for each node since all nodes in parent tables are done as 1,2,3,4... not real parent names
            numberofparentsvec[updatenodes]<-newptab$numberofparentsvec[updatenodes]
            numparents[updatenodes]<-newptab$numparents[updatenodes]
            newplus1lists<-PLUS1(n,aliases,updatenodes,blacklistparents)
            plus1lists$mask[updatenodes]<- newplus1lists$mask[updatenodes]
            plus1lists$parents[updatenodes]<- newplus1lists$parents[updatenodes]
            plus1lists$aliases[updatenodes]<- newplus1lists$aliases[updatenodes]

            rowmaps[updatenodes]<-parentsmapping(parenttable,numberofparentsvec,n,updatenodes)[updatenodes]
            scoretable[updatenodes]<-scorepossibleparents.PLUS1(parenttable,plus1lists,n,param,updatenodes,rowmaps,numparents,numberofparentsvec)[updatenodes]
            posetparenttable[updatenodes]<-poset(parenttable,numberofparentsvec,rowmaps,n,updatenodes)[updatenodes]
            if (MAP) {
              newmaxmatrices<-posetscoremax(posetparenttable,scoretable,numberofparentsvec,
                                            rowmaps,n,plus1lists=plus1lists,updatenodes)
              maxmatrices$maxmatrix[updatenodes]<- newmaxmatrices$maxmatrix[updatenodes]
              maxmatrices$maxrow[updatenodes]<- newmaxmatrices$maxrow[updatenodes]
            } else {
              newbannedscore<-poset.scores(posetparenttable,scoretable,numberofparentsvec,rowmaps,
                                           n,plus1lists=plus1lists,numparents,updatenodes)
              bannedscore[updatenodes]<-newbannedscore[updatenodes]
            }
            if(verbose) {
              print(paste("MCMC plus1 iteration",i))
              flush.console()
              }
          } else {
            if(verbose) {
            print(paste("score tables calculated, MCMC plus1 starts"))
            flush.console()
          }
        }
          if(MAP) {
            MCMCresult<-orderMCMCplus1max(n,startorder,iterations,stepsave,moveprobs,parenttable,
                                          scoretable,aliases,numparents,rowmaps,plus1lists,maxmatrices,numberofparentsvec,
                                          gamma=gamma)
            endtime<-Sys.time()
            if(verbose) {
              print(endtime-starttime)
              flush.console()
              }
          } else {
            MCMCresult<-orderMCMCplus1(n,startorder,iterations,stepsave,moveprobs,parenttable,
                                       scoretable,aliases,numparents,rowmaps,plus1lists,
                                       bannedscore,numberofparentsvec,gamma=gamma)
            endtime<-Sys.time()
            if(verbose){
              print(endtime-starttime)
              flush.console()
              }
          }

          MCMCchains[[i]]<-MCMCresult
          maxN<-which.max(unlist(MCMCresult[[2]]))
          if(i>1){
            if (MCMCresult[[2]][[maxN]]>maxscore){
              maxDAG<-MCMCresult[[1]][[maxN]]
              maxorder<-MCMCresult[[4]][[maxN]]
              maxscore<-MCMCresult[[2]][[maxN]]
              maxit<-i
            }
          } else {
            maxDAG<-MCMCresult[[1]][[maxN]]
            maxorder<-MCMCresult[[4]][[maxN]]
            maxscore<-MCMCresult[[2]][[maxN]]
            maxit<-1
          }
          if (MAP) {
            newadj<- newspacemap(n,startskeleton,oldadj,softlimit,hardlimit,blacklist,maxN=maxN,MCMCchain=MCMCresult[[1]],mergetype=mergecp)
          } else {
            newadj<- newspaceskel(n,startskeleton,oldadj,softlimit,hardlimit,posterior,blacklist,MCMCchain=MCMCresult[[1]],mergetype=mergecp)
          }
          updatenodes<-which(apply(newadj==oldadj,2,all)==FALSE)
          updatenodeslist[[i]]<-updatenodes
          oldadj<-newadj
          startorder<-MCMCresult[[4]][[maxN]]
          i<-i+1
        }

      }
      
      maxres<-list()
      maxres$DAG<-maxDAG
      maxres$order<-maxorder
      maxres$score<-maxscore
      maxres$it<-maxit

      colnames(maxres$DAG)<-colnames(param$data)

      result<-list()
      if (scoreout){
        if(chainout){output<-4}
        else{output<-3}
      } else {
        if(chainout) {output<-2}
        else {output<-1}
      }

      switch(as.character(output),
             "1"={ # return only maximum DAG and order
               result$max<-maxres
               attr(result$max,"class")<-"MCMCmax"
               return(result)
             },
             "2"={ # return all MCMC all saved MCMC steps: incidence, DAGscore, orderscore and order and max result
               result$max<-maxres
               attr(result$max,"class")<-"MCMCmax"
               if (i==1) {
                 result$chain<-MCMCresult
                 attr(result$chain,"class")<-"MCMCchain"
               } else {
                 result$chain<-MCMCchains
                 attr(result$chain,"class")<-"MCMCmult"}
               return(result)
             },
             "3"={ # return max DAG, order, last search space incidence and all scoretables
               result$max<-maxres
               attr(result$max,"class")<-"MCMCmax"
               result$space<-list()
               result$space$adjacency<-oldadj
               result$space$scoretable<-scoretable
               attr(result$space,"class")<-"MCMCspace"
               return(result)
             },
             "4"={ # return all MCMC all saved MCMC steps,max result,last search space and scoretables
               result$max<-maxres
               attr(result$max,"class")<-"MCMCmax"
               result$space<-list()
               result$space$adjacency<-oldadj
               result$space$scoretable<-scoretable
               attr(result$space,"class")<-"MCMCspace"
               if(i==1)
               {result$chain<-MCMCresult
               attr(result$chain,"class")<-"MCMCchain"} else {
                 result$chain<-MCMCchains
                 attr(result$chain,"class")<-"MCMCmult"}
               return(result)
             }
      )
    }

}
