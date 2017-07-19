weightedbinCItest <- function(x,y,S,suffStat) {
  ## p-value:
  CIcoretest(x, y, S, suffStat)
}



CIcoretest<-function(j,k,parentnodes,suffStat){
  
  lp<-length(parentnodes) # number of parents
  noparams<-2^lp # number of binary states of the parents
  
  switch(as.character(lp),
         "0"={# no parents
           N1<-sum(suffStat$d1[j,])
           N0<-sum(suffStat$d0[j,])
           NT<-N0+N1
           M1<-sum(suffStat$d1[k,])
           M0<-sum(suffStat$d0[k,])
           MT<-M0+M1
           P11<-sum(suffStat$d1[j,]*suffStat$data[k,])
           P10<-sum(suffStat$d1[j,]*(1-suffStat$data[k,]))
           P01<-sum(suffStat$d0[j,]*suffStat$data[k,])
           P00<-sum(suffStat$d0[j,]*(1-suffStat$data[k,]))
           PT<-P11+P10+P01+P00
           
           # calculate the statistic in the tedious way    
           
           if(P11>0){
             part1<-P11*(log(P11)+log(PT)-log(N1)-log(M1))
           } else{
             part1<-0
           }
           if(P10>0){
             part2<-P10*(log(P10)+log(PT)-log(N1)-log(M0))
           } else{
             part2<-0
           }
           if(P01>0){
             part3<-P01*(log(P01)+log(PT)-log(N0)-log(M1))
           } else{
             part3<-0
           }
           if(P00>0){
             part4<-P00*(log(P00)+log(PT)-log(N0)-log(M0))
           } else{
             part4<-0
           }
           
           Gsquared<-2*(part1+part2+part3+part4)
           
         },
         "1"={# one parent
           Gsquared<-0 
           summys<-suffStat$data[parentnodes,]
           
           for(i in 1:noparams-1){
             totest<-which(summys==i)
             N1<-sum(suffStat$d1[j,totest])
             N0<-sum(suffStat$d0[j,totest])
             NT<-N0+N1
             M1<-sum(suffStat$d1[k,totest])
             M0<-sum(suffStat$d0[k,totest])
             MT<-M0+M1
             P11<-sum(suffStat$d1[j,totest]*suffStat$data[k,totest])
             P10<-sum(suffStat$d1[j,totest]*(1-suffStat$data[k,totest]))
             P01<-sum(suffStat$d0[j,totest]*suffStat$data[k,totest])
             P00<-sum(suffStat$d0[j,totest]*(1-suffStat$data[k,totest]))
             PT<-P11+P10+P01+P00
             
             # calculate the statistic in the tedious way    
             
             if(P11>0){
               part1<-P11*(log(P11)+log(PT)-log(N1)-log(M1))
             } else{
               part1<-0
             }
             if(P10>0){
               part2<-P10*(log(P10)+log(PT)-log(N1)-log(M0))
             } else{
               part2<-0
             }
             if(P01>0){
               part3<-P01*(log(P01)+log(PT)-log(N0)-log(M1))
             } else{
               part3<-0
             }
             if(P00>0){
               part4<-P00*(log(P00)+log(PT)-log(N0)-log(M0))
             } else{
               part4<-0
             }
             
             Gsquared<-Gsquared+2*(part1+part2+part3+part4)
             
           }
         },     
         { # more parents
           Gsquared<-0  
           summys<-colSums(2^(c(0:(lp-1)))*suffStat$data[parentnodes,])
           
           for(i in 1:noparams-1){
             totest<-which(summys==i)
             N1<-sum(suffStat$d1[j,totest])
             N0<-sum(suffStat$d0[j,totest])
             NT<-N0+N1
             M1<-sum(suffStat$d1[k,totest])
             M0<-sum(suffStat$d0[k,totest])
             MT<-M0+M1
             P11<-sum(suffStat$d1[j,totest]*suffStat$data[k,totest])
             P10<-sum(suffStat$d1[j,totest]*(1-suffStat$data[k,totest]))
             P01<-sum(suffStat$d0[j,totest]*suffStat$data[k,totest])
             P00<-sum(suffStat$d0[j,totest]*(1-suffStat$data[k,totest]))
             PT<-P11+P10+P01+P00
             
             # calculate the statistic in the tedious way    
             
             if(P11>0){
               part1<-P11*(log(P11)+log(PT)-log(N1)-log(M1))
             } else{
               part1<-0
             }
             if(P10>0){
               part2<-P10*(log(P10)+log(PT)-log(N1)-log(M0))
             } else{
               part2<-0
             }
             if(P01>0){
               part3<-P01*(log(P01)+log(PT)-log(N0)-log(M1))
             } else{
               part3<-0
             }
             if(P00>0){
               part4<-P00*(log(P00)+log(PT)-log(N0)-log(M0))
             } else{
               part4<-0
             }
             
             Gsquared<-Gsquared+2*(part1+part2+part3+part4)
             
           }
         })
  
  pvally<-pchisq(Gsquared, noparams, lower.tail = FALSE)
  
  return(pvally)
}