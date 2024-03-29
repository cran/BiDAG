#implements order MCMC on a defined search space, sampling version
#partly derived from <doi:10.1080/01621459.2015.1133426>
orderMCMCbase<-function(n,nsmall,startorder,iterations,stepsave,moveprobs,parenttable,scoretable,aliases,numparents,
                        rowmaps,scoresmatrices,numberofparentsvec,gamma=1,bgnodes,matsize,chainout=FALSE,compress=TRUE) {
  
  result<-list()
  
  if(!is.null(bgnodes)) {
    mainnodes<-c(1:matsize)[-bgnodes]
  } else {mainnodes<-c(1:matsize)}
  
  currentpermy<-startorder #starting order represented as a permutation
  currentorderscores<-orderscoreBase(matsize,currentpermy[1:nsmall],c(1:nsmall),parenttable,aliases,numparents,
                                     rowmaps,scoretable,scoresmatrices,currentpermy) #starting score
  currenttotallogscore<-sum(currentorderscores$totscores[mainnodes]) #log total score of all DAGs in the starting order
  currentDAG<-samplescoreplus1(matsize,mainnodes,currentorderscores,plus1lists=NULL,scoretable,
                               scoresmatrices,parenttable,numberofparentsvec,aliases) #score of a single DAG sampled from the starting order

  L1 <- list() # stores the adjacency matrix of a DAG sampled from the orders
  L2 <- vector() # stores its log BGe score
  L3 <- vector() # stores the log BGe score of the entire order
  L4 <- list() # stores the orders as permutations

  zlimit<- floor(iterations/stepsave) + 1 # number of outer iterations
  length(L1) <- zlimit
  length(L2) <- zlimit
  length(L3) <- zlimit
  length(L4) <- zlimit

  L1[[1]]<-currentDAG$incidence #starting DAG adjacency matrix
  L2[1]<-currentDAG$logscore #starting DAG score
  L3[1]<-currenttotallogscore #starting order score
  L4[[1]]<-currentpermy[1:nsmall] #starting order
  maxdag<-currentDAG$incidence
  maxscore<-L2[1]

  moveprobsstart<-moveprobs

  for (z in 2:zlimit){ #the MCMC chain loop with 'iteration' steps is in two parts
    for (count in 1:stepsave){ #since we only save the results to the lists each 'stepsave'


      chosenmove<-sample.int(4,1,prob=moveprobs)
      if(chosenmove<4){  # if it is 3 then we stay still
        
        proposedpermy<-currentpermy #sample a new order by swapping two elements
        switch(as.character(chosenmove),
               "1"={ # swap any two elements at random
                 sampledelements<-sample.int(nsmall,2,replace=FALSE) #chosen at random
                 
               },
               "2"={ # swap any adjacent elements
                 k<-sample.int(nsmall-1,1) #chose the smallest at random
                 sampledelements<-c(k,k+1)
                 
               },
               "3"={ # swap any adjacent elements
                 sampledpos<-sample.int(nsmall,1)
               },
               {# if neither is chosen, we have a problem
                 cat('The move sampling has failed!')
               }) #end switch

        if(chosenmove<3){
          
          proposedpermy[sampledelements]<-currentpermy[rev(sampledelements)] #proposed new order ???
          scorepositions<-c(min(sampledelements):max(sampledelements))
          rescorenodes<-proposedpermy[scorepositions] #we only need to rescore these nodes between the swapped elements to speed up the calculation
          proposedorderrescored<-orderscoreBase(matsize,rescorenodes,scorepositions,parenttable,aliases,numparents,rowmaps,
                                                scoretable,scoresmatrices,proposedpermy)
          proposedtotallogscore<-currenttotallogscore-sum(currentorderscores$totscores[rescorenodes])+sum(proposedorderrescored$totscores[rescorenodes]) #and the new log total score by updating only the necessary nodes
          scoreratio<-exp((proposedtotallogscore-currenttotallogscore)*gamma) #acceptance probability

          if(runif(1)<scoreratio){ #Move accepted then set the current order and scores to the proposal
            
            currentpermy<-proposedpermy
            currentorderscores$therow[rescorenodes]<-proposedorderrescored$therow[rescorenodes]
            currentorderscores$totscores[rescorenodes]<-proposedorderrescored$totscores[rescorenodes]
            currenttotallogscore<-proposedtotallogscore

          }
        } else {
          
              neworder<-positionscorebase(matsize,nsmall,currentorderscores,sampledpos,currentpermy,aliases,
                                           rowmaps,numparents,scoretable,scoresmatrices)
              currentpermy<-neworder$order
              currentorderscores<-neworder$score
              currenttotallogscore<-sum(neworder$totscores)
          }

      }
    }
    currentDAG<-samplescoreplus1(matsize,mainnodes,currentorderscores,plus1lists=NULL,scoretable,scoresmatrices,
                                 parenttable,numberofparentsvec,aliases)
    if(chainout) {
    if(compress) {
      L1[[z]]<-Matrix(currentDAG$incidence,sparse=TRUE) #store compressed adjacency matrix of a sampled DAG each 'stepsave'
    } else {
      L1[[z]]<-currentDAG$incidence #store adjacency matrix of a sampled DAG each 'stepsave'
    }
    }
    
    L2[z]<-currentDAG$logscore #and log score of a sampled DAG
    L3[z]<-currenttotallogscore #and the current order score
    L4[[z]]<-currentpermy[1:nsmall] #and store current order
    if(L2[z]>maxscore) {
      maxscore<-L2[z]
      maxdag<-currentDAG$incidence
    } 
  }
  result$incidence<-L1
  result$DAGscores<-L2
  result$orderscores<-L3
  result$orders<-L4
  result$maxscore<-maxscore
  result$maxdag<-Matrix(maxdag,sparse=TRUE)
  
  return(result)
}

