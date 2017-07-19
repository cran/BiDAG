
dettwobytwo <- function(D) {
  D[1,1]*D[2,2]-D[1,2]*D[2,1]
}

# The determinant of a 3 by 3 matrix
detthreebythree <- function(D){
  D[1,1]*(D[2,2]*D[3,3]-D[2,3]*D[3,2])-D[1,2]*(D[2,1]*D[3,3]-D[2,3]*D[3,1])+D[1,3]*(D[2,1]*D[2,3]-D[2,2]*D[3,1])
}

# The log of the BGe/BDe score, but simplified as much as possible
# see arXiv:1402.6863 
DAGcorescore<-function(j,parentnodes,n,param) {

switch(param$type,
"bge" = {
    TN<-param$TN
    awpN<-param$awpN
    scoreconstvec<-param$scoreconstvec
    
      lp<-length(parentnodes) #number of parents
      awpNd2<-(awpN-n+lp+1)/2
      A<-TN[j,j]
      
      switch(as.character(lp),
             "0"={# just a single term if no parents
               corescore <- scoreconstvec[lp+1] -awpNd2*log(A)
             },
             
             "1"={# no need for matrices
               D<-TN[parentnodes,parentnodes]
               logdetD<-log(D)
               B<-TN[j,parentnodes]
               logdetpart2<-log(A-B^2/D)
               corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
             },
             
             "2"={# can do matrix determinant and inverse explicitly
               # but this is numerically unstable for large matrices!
               # so we use the same approach as for 3 parents
               D<-TN[parentnodes,parentnodes]
               detD<-dettwobytwo(D)
               logdetD<-log(detD)
               B<-TN[j,parentnodes]
               #logdetpart2<-log(A-(D[2,2]*B[1]^2+D[1,1]*B[2]^2-2*D[1,2]*B[1]*B[2])/detD) #also using symmetry of D
               logdetpart2<-log(dettwobytwo(D-(B)%*%t(B)/A))+log(A)-logdetD
               corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
             },
             
             {# otherwise we use cholesky decomposition to perform both
               D<-as.matrix(TN[parentnodes,parentnodes])
               choltemp<-chol(D)
               logdetD<-2*log(prod(choltemp[(lp+1)*c(0:(lp-1))+1]))
               B<-TN[j,parentnodes]
               logdetpart2<-log(A-sum(backsolve(choltemp,B,transpose=TRUE)^2))
               corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
             })},
"bde" = {
  lp<-length(parentnodes) # number of parents
  noparams<-2^lp # number of binary states of the parents
  chi<-param$chi
  scoreconstvec<-param$scoreconstvec
  
  switch(as.character(lp),
         "0"={# no parents
           N1<-sum(param$d1[j,],na.rm=TRUE)
           N0<-sum(param$d0[j,],na.rm=TRUE)
           NT<-N0+N1
           corescore <- scoreconstvec[lp+1] + lgamma(N0+chi/(2*noparams)) + lgamma(N1+chi/(2*noparams)) - lgamma(NT+chi/noparams)
         },
         "1"={# one parent
           corescore<-scoreconstvec[lp+1]  
           summys<-param$data[parentnodes,]
           for(i in 1:noparams-1){
             totest<-which(summys==i)
             N1<-sum(param$d1[j,totest],na.rm=TRUE)
             N0<-sum(param$d0[j,totest],na.rm=TRUE)
             NT<-N0+N1
             corescore <- corescore + lgamma(N0+chi/(2*noparams)) + lgamma(N1+chi/(2*noparams)) - lgamma(NT+chi/noparams)
           }
         },     
         { # more parents
           summys<-colSums(2^(c(0:(lp-1)))*param$data[parentnodes,])
           N1s<-collectC(summys,param$d1[j,],noparams)
           N0s<-collectC(summys,param$d0[j,],noparams)   
           NTs<-N1s+N0s       
           corescore <- scoreconstvec[lp+1] + sum(lgamma(N0s+chi/(2*noparams))) + 
           sum(lgamma(N1s+chi/(2*noparams))) - sum(lgamma(NTs+chi/noparams))
         })
  

  }
)
  return(corescore)
}