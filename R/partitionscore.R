partitionscoreplus1<-function(n,scorenodes,parenttable,aliases,scoretable,plus1lists,rowmapsneeded,rowmapsallowed,
                              scoresneeded,scoresallowed,permy,party,posy){

  partitionscores<-rep(0,n)
  scoresvec<-vector("numeric")
  therow1<-vector("numeric")
  therow2<-vector("numeric")
  m<-length(party)
  allowedlists1<-list()
  allowedlists2<-list()
  for (i in scorenodes){
    tablesize<-dim(parenttable[[i]]) # just to remove some arguments
    position<-which(permy==i)
    partyelement<-posy[position]

    if(partyelement==m){# no parents are allowed
      partitionscores[i]<-scoretable[[i]][[1]][1,1]
      therow1[i]<-0
      therow2[i]<-0
    }
    else if(tablesize[1]==0) { #no parents in PC-table, then only plus1 lists with requirednodes are good
      requirednodes<-permy[which(posy==(partyelement+1))]
      allowedneeded<-intersect(requirednodes,aliases[[i]])
      therow1[i]<-0
      therow2[i]<-1
      allowedlists2[[i]]<-which(plus1lists$parents[[i]]%in%requirednodes)+1
      scoresvec<-scoresneeded[[i]][therow2[i],allowedlists2[[i]]]
      maxallowed<-max(scoresvec)
      partitionscores[i]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
    }

    else    {
      bannednodes<-permy[which(posy<=partyelement)]
      allowednodes<-permy[which(posy>partyelement)]
      requirednodes<-permy[which(posy==(partyelement+1))]
      bannedpool<-which(aliases[[i]]%in%bannednodes) #new
      allowedneeded<-which(aliases[[i]]%in%requirednodes)

      if(length(bannedpool)==tablesize[2]){ #no parents from the main table allowed
        therow1[i]<-0
        therow2[i]<-nrow(parenttable[[i]])
        allowedlists2[[i]]<-which(plus1lists$parents[[i]]%in%requirednodes)+1 #plus1 lists for nodes which are in required set
        scoresvec<-scoresneeded[[i]][therow2[i],allowedlists2[[i]]]
        maxallowed<-max(scoresvec)
        partitionscores[i]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
      }
      else if(length(bannedpool)<tablesize[2]&&length(allowedneeded)==0){ #no nodes from main tables are in required nodes set,
        #but sume are in the allowed set
        therow1[i]<-0
        therow2[i]<-rowmapsneeded[[i]]$backwards[sum(2^bannedpool)/2+1]
        allowedlists2[[i]]<-which(plus1lists$parents[[i]]%in%requirednodes)+1
        scoresvec<-scoresneeded[[i]][therow2[i],allowedlists2[[i]]]
        maxallowed<-max(scoresvec)
        partitionscores[i]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
      }
      else if (length(aliases[[i]])==1){ #only 1 parent in possible parent set
        therow1[i]<-1
        therow2[i]<-1
        allowedlists2[[i]]<-c(which(plus1lists$parents[[i]]%in%requirednodes)+1)
        allowedlists1[[i]]<-c(1,which(plus1lists$parents[[i]]%in%allowednodes)+1)
        scoresvec<-scoresneeded[[i]][therow2[i],allowedlists2[[i]]]
        scoresvec<-c(scoresvec,scoresallowed[[i]][therow1[i],allowedlists1[[i]]])
        maxallowed<-max(scoresvec)
        partitionscores[i]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
      }
      else {
        therow2[i]<-rowmapsneeded[[i]]$backwards[sum(2^bannedpool)/2+1]
        if (length(bannedpool)==0){
          therow1[i]<-rowmapsallowed[[i]]$backwards[2*sum(3^allowedneeded)/3+1]
        } else {
          therow1[i]<-rowmapsallowed[[i]]$backwards[sum(3^bannedpool)/3+2*sum(3^allowedneeded)/3+1]}
        allowedlists2[[i]]<-which(plus1lists$parents[[i]]%in%requirednodes)+1
        allowedlists1[[i]]<-c(1,which(plus1lists$parents[[i]]%in%allowednodes)+1)
        scoresvec<-scoresneeded[[i]][therow2[i],allowedlists2[[i]]]
        scoresvec<-c(scoresvec,scoresallowed[[i]][therow1[i],allowedlists1[[i]]])
        maxallowed<-max(scoresvec)
        partitionscores[i]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
      }
    }
  }

  scores<-list()
  scores$therow1<-therow1
  scores$therow2<-therow2
  scores$allowedlists1<-allowedlists1
  scores$allowedlists2<-allowedlists2
  scores$totscores<-partitionscores
  #print(scores$totscores[scorenodes])
  return(scores)


}



