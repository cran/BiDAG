# this function scores by propagating through the partition poset graph

posetpartitionscoreparents<-function(numberofparents,neededposetparenttable,neededparentsvec,
                                     numberofparentsvec,rowmaps,needednodebannedrow,scoretable,n) {
  
  scorematrices<-list()

  for (j in 1:n)  {
    np<-numberofparents[j]
    if(np==1) {
      scorematrices[[j]]<-as.matrix(scoretable[[j]][2,1])
    } else if (np>1){
      binomcoefs<-choose(np,c(np:1))*(2^c(np:1)-1)
      nrows<-nrow(neededposetparenttable[[j]])
      P_local <- numeric(nrows)
      nrowold<-length(rowmaps[[j]]$forward) # size of the other poset graph
      # at the top of the graph we only allow one possible parent set

      for (i in nrows:(nrows-np+1)) {
        k <- needednodebannedrow[[j]][i] # the banned nodes row, there should be (n-1) banned nodes
        conjugatescore<-scoretable[[j]][rowmaps[[j]]$backwards[nrowold-rowmaps[[j]]$forward[k]+1],1]
        P_local[i]<-conjugatescore
      }

      # for each other level of the poset graph we need to add the parents divided by the difference in levels

      cutoff<-1

      for(level in 0:(np-2)) {

        for (i in (nrows-np):cutoff)  { # find the parents in the poset graph
          k <- needednodebannedrow[[j]][i] # the banned nodes row
          parentnodes <- neededposetparenttable[[j]][i,c(1:neededparentsvec[[j]][i])]
          maxparents<-max(P_local[parentnodes])
          # take the sum of the parent scores and divide by the relevant factor
          parentsum<-log(sum(exp(P_local[parentnodes]-maxparents)))+maxparents -
            log(np-rev(numberofparentsvec[[j]])[k]-level+1)
          # take the conjugate score
          conjugatescore<-scoretable[[j]][rowmaps[[j]]$backwards[nrowold-rowmaps[[j]]$forward[k]+1],1]
          # find max and combine
          maxoverall<-max(parentsum,conjugatescore)
          P_local[i]<- log(exp(parentsum-maxoverall)+exp(conjugatescore-maxoverall)) + maxoverall
        }

        cutoff<-cutoff+binomcoefs[level+1]
      }
      scorematrices[[j]]<-as.matrix(P_local)
    }
  }
  return(scorematrices)
}

#builds a banned score table

plus1allowed.partition<-function(numparents,parenttable,neededposetparenttable,neededparentsvec,
                                 numberofparentsvec,rowmapsallowed,needednodebannedrow,scoretable,plus1lists,n) {
  
  revnumberofparentsvec<-lapply(numberofparentsvec,rev)
  rowmaps<-rowmapsallowed
  if (is.null(plus1lists)) {
    scorematrices<-list()
    for (j in 1:n){
      np<-numparents[j]
      if(np==1) {
        scorematrices[[j]]<-as.matrix(scoretable[[j]][2,1])
      } else if (np>1) {
        binomcoefs<-choose(np,c(np:1))*(2^c(np:1)-1)
        nrows<-nrow(neededposetparenttable[[j]])
        P_local <- numeric(nrows)
        nrowold<-length(rowmaps[[j]]$forward) # size of the other poset graph
        # at the top of the graph we only allow one possible parent set
        for (i in nrows:(nrows-np+1)) {
          k <- needednodebannedrow[[j]][i] # the banned nodes row, there should be (n-1) banned nodes
          conjugatescore<-scoretable[[j]][rowmaps[[j]]$backwards[nrowold-rowmaps[[j]]$forward[k]+1],1]
          P_local[i]<-conjugatescore
        }
        # for each other level of the poset graph we need to add the parents divided by the difference in levels
        cutoff<-1
        for(level in 0:(np-2)) {
          for (i in (nrows-np):cutoff)  { # find the parents in the poset graph
            k <- needednodebannedrow[[j]][i] # the banned nodes row
            parentnodes <- neededposetparenttable[[j]][i,c(1:neededparentsvec[[j]][i])]
            maxparents<-max(P_local[parentnodes])
            # take the sum of the parent scores and divide by the relevant factor
            parentsum<-log(sum(exp(P_local[parentnodes]-maxparents)))+maxparents - log(np-revnumberofparentsvec[[j]][k]-level+1)
            # take the conjugate score
            conjugatescore<-scoretable[[j]][rowmaps[[j]]$backwards[nrowold-rowmaps[[j]]$forward[k]+1],1]
            # find max and combine
            maxoverall<-max(parentsum,conjugatescore)
            P_local[i]<- log(exp(parentsum-maxoverall)+exp(conjugatescore-maxoverall)) + maxoverall
          }
          cutoff<-cutoff+binomcoefs[level+1]
        }
        scorematrices[[j]]<-as.matrix(P_local)
      }
    }
    return(scorematrices)
  } else{
    scorematrices.allowed<-list()
    for (j in 1:n){
      np<-numparents[j] #number of possible parents for node j
      binomcoefs<-choose(np,c(np:1))*(2^c(np:1)-1)
      ll<-length(plus1lists$parents[[j]])+1
      nrows<-nrow(neededposetparenttable[[j]])
      nrowspar<-nrow(parenttable[[j]])
      if(np>0) {P_local <- matrix(nrow=nrows,ncol=ll)}
      for (li in 1:ll){
        if (np==0){
          scorematrices.allowed[[j]]<-NULL
          break #we don't have allowed table, all scores already in the scoretable plus1 lists
        }
        else if(np==1) {
          #we need just 1 additional parent set which is not in the scoretable plus1 node+the only parent in scoretable
          P_local[1,li]<-as.matrix(scoretable[[j]][[li]][2,1])
          scorematrices.allowed[[j]]<-P_local
        } else if (np>1){
          nrowold<-length(rowmaps[[j]]$forward) # size of the other poset graph
          # at the top of the graph we only allow one possible parent set
          for (i in nrows:(nrows-np+1)){
            k <- needednodebannedrow[[j]][i] # the banned nodes row, there should be (n-1) banned nodes
            conjugatescore<-scoretable[[j]][[li]][rowmaps[[j]]$backwards[nrowold-rowmaps[[j]]$forward[k]+1],1]
            P_local[i,li]<-conjugatescore
          }

          # for each other level of the poset graph we need to add the parents divided by the difference in levels
          cutoff<-1
          for(level in 0:(np-2)){
            for (i in (nrows-np):cutoff)  { # find the parents in the poset graph
              k <- needednodebannedrow[[j]][i] # the banned nodes row
              parentnodes <- neededposetparenttable[[j]][i,c(1:neededparentsvec[[j]][i])]
              maxparents<-max(P_local[parentnodes])
              # take the sum of the parent scores and divide by the relevant factor
              parentsum<-log(sum(exp(P_local[parentnodes]-maxparents)))+maxparents - log(np-revnumberofparentsvec[[j]][k]-level+1)
              # take the conjugate score
              conjugatescore<-scoretable[[j]][[li]][rowmaps[[j]]$backwards[nrowold-rowmaps[[j]]$forward[k]+1],1]
              # find max and combine
              maxoverall<-max(parentsum,conjugatescore)
              P_local[i,li]<- log(exp(parentsum-maxoverall)+exp(conjugatescore-maxoverall)) + maxoverall
            }

            cutoff<-cutoff+binomcoefs[level+1]
          }
          scorematrices.allowed[[j]]<-as.matrix(P_local)
        }

      }
    }
    return(scorematrices.allowed)
  }

}

parentlistnonempty<-function(elements,n){
  
  matrixofparents<-matrix(NA,nrow=(2^length(elements)-1),ncol=n)
  
  cutoff<-0
  
  for (r in 1:length(elements)){
    possparents<-combinations(length(elements),r,elements)
    heighty<-nrow(possparents)
    matrixofparents[(1:heighty)+cutoff,1:r]<-possparents
    cutoff<-cutoff+heighty
  }
  return(matrixofparents)
}

partitionlistnumberofparents<-function(needednodetable,numberofparentsvec,n){
  numberofpartitionparentsvec<-list()
  
  for (i in 1:n){
    
    if(length(numberofparentsvec[[i]])==1) {numberofpartitionparentsvec[[i]]<-0}
    else {
      nrows<-length(numberofparentsvec[[i]])
      npar<-numberofparentsvec[[i]][nrows]
      P_local<-c()
      for (j in 1:(nrows-1))  { # last row has no possibilities
        np<-npar-numberofparentsvec[[i]][j]
        P_local<-c(P_local,rep(c(1:np),choose(np,c(1:np))))
      }
      numberofpartitionparentsvec[[i]]<-P_local
    }
  }
  return(numberofpartitionparentsvec)
}

partitionmapneedednodebannedrow<-function(numberofparents,numberofparentsvec,n) {
  needednodebannedrow<-list()
  for (i in 1:n) { 
    j<-numberofparents[i]
    needednodebannedrow[[i]]<-rep(1:2^j,2^(j-numberofparentsvec[[i]])-1)
  }
  return(needednodebannedrow)
}

partitionlist<-function(parenttable,numberofparentsvec,n) {
  needednodetable<-list()
  for  (i in 1:n) {
    cutoff<-0
    nrows<-nrow(parenttable[[i]])
    ncols<-ncol(parenttable[[i]])
    matrixofparents<-matrix(NA,nrow=(3^ncols-2^ncols),ncol=ncols)
    if(nrows>1) {
      for (j in 1:(nrows-1)){ # last row has no possibilities
        parentnodes <- parenttable[[i]][j,1:numberofparentsvec[[i]][j]]
        elements<-setdiff(c(1:ncols),parentnodes)
        newpart<-parentlistnonempty(elements,ncols)
        heighty<-nrow(newpart)
        matrixofparents[(1:heighty)+cutoff,]<-newpart
        cutoff<-cutoff+heighty
      }
    }
    needednodetable[[i]]<-as.matrix(matrixofparents)
  }
  return(needednodetable)
}

neededparentsmapping<-function(parenttable,numberofparentsvec,needednodetable,numberofpartitionparentsvec,needednodebannedrow,n) {
  
  maps<-list()
  for (i in 1:n){
    nrows<-nrow(needednodetable[[i]])
    ncols<-ncol(needednodetable[[i]])
    maps[[i]]<-list()
    P_local <- numeric(3^ncols) # we leave some elements empty
    P_localinv <- numeric(3^ncols) # here too
    if(nrows>1){
      for (j in 1:nrows)  {# the needed nodes
        needednodes <- needednodetable[[i]][j,c(1:numberofpartitionparentsvec[[i]][j])]
        k <- needednodebannedrow[[i]][j] # the banned nodes row
        if(k>1) { # if there is at least one banned node
          bannednodes <- parenttable[[i]][k,c(1:numberofparentsvec[[i]][k])]
          P_local[j]<-sum(3^bannednodes)/3+2*sum(3^needednodes)/3+1
        } else {P_local[j]<-2*sum(3^needednodes)/3+1}
        P_localinv[P_local[j]]<-j # the inverse mapping
      }
    }
    maps[[i]]$forward<-P_local
    maps[[i]]$backwards<-P_localinv
  }
  return(maps)
}

needed.poset<-function(parenttable,numberofparentsvec,needednodebannedrow,neededrowmaps,n) { 
  posetparents<-list()
  for (i in 1:n) {
    nrows<-length(needednodebannedrow[[i]])
    numpar<-numberofparentsvec[[i]][nrow(parenttable[[i]])]
    if(nrows>1) {
      posetneededtable<-matrix(NA,nrow=nrows,ncol=numpar)
      offsets<-rep(1,nrows)
      if(nrows>1) {
        for(j in (nrows:(2^numpar))){
          k <- needednodebannedrow[[i]][j] # the banned nodes row, there should be at least one banned node
          bannednodes <- parenttable[[i]][k,c(1:numberofparentsvec[[i]][k])]
          # we can either delete one of the banned nodes
          children1<-neededrowmaps[[i]]$backwards[neededrowmaps[[i]]$forward[j]-3^bannednodes/3]
          # move one of the banned nodes to a needed noded
          children2<-neededrowmaps[[i]]$backwards[neededrowmaps[[i]]$forward[j]+3^bannednodes/3]
          children<-c(children1,children2)
          posetneededtable[cbind(children,offsets[children])]<-j
          offsets[children]<-offsets[children]+1
        }
      }
    }
    posetparents[[i]]<-list()
    if (numpar==1) {
      posetneededtable<-matrix(NA,nrow=1,ncol=1)
      posetneededtable[1,1]<-1
      posetparents[[i]]$sizes<-c(1)
      posetparents[[i]]$table<-posetneededtable
    }
    else if (numpar>1) {
      posetparents[[i]]$sizes<-offsets-1
      posetparents[[i]]$table<-posetneededtable}
    else {
      posetparents[[i]]$sizes<-vector("integer",length=0)
      posetparents[[i]]$table<-NULL}
  }
  return(posetparents)
}