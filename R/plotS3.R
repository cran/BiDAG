#PLOT methods for classes:
#orderMCMC
#partitionMCMC
#iterativeMCMC
#itercomp
#samplecomp


#' Plotting object of class 'iterativeMCMC'
#'
#' @param x object of class 'iterativeMCMC'
#' @param ... ignore
#' @param main name of the graph; "iterative MCMC, DAG scores" by default
#' @param xlab name of x-axis; "MCMC step"
#' @param ylab name of y-axis; "DAG logscore"
#' @param type type of line in the plot; "l" by default
#' @param col colour of line in the plot; "blue" by default
#' 
#' @rdname iterativeMCMC
#' @method plot iterativeMCMC
#' @export
plot.iterativeMCMC <-function(x,...,main="iterative MCMC, DAG scores", xlab="MCMC step", ylab="DAG logscore", type="l", col="blue"){
  old.par <- par(no.readonly = TRUE)
  x<-x$trace
  nchains<-length(x)
  nsteps<-length(x[[1]])
  scorevec<-Reduce('c',x)
  
  scoremax<-max(scorevec)
  scoremin<-min(scorevec)
  # Add extra space to right of plot area; change clipping to figure
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(scorevec,type=type,xlab=xlab,ylab=ylab, ...,
       ylim=c(scoremin,scoremax),main=main,
       col=col)
  legend(length(scorevec)*1.04, scoremin+(scoremax-scoremin)*0.2,col="grey", lty=2,angle=90, 
         legend=c("search \nspace \nexpansion"),bty="n")
  
  par(xpd=FALSE)
  for(i in 1:(nchains-1)) {
    abline(v=i*nsteps,col="grey",lty=2,ylim=c(scoremin,scoremax))
  }
  par(xpd=TRUE)
  par(old.par)
}



#' Plotting object of class 'orderMCMC'
#'
#' @param x object of class 'orderMCMC'
#' @param ... other parameters to be passed through to plotting functions
#' @param burnin number between \code{0} and \code{1}, indicates the percentage of the samples which will be discarded as `burn-in' of the MCMC chain; the rest  of the samples will be used to calculate the posterior probabilities; 0.2 by default
#' @param main name of the graph; "DAG logscores" by default
#' @param xlab name of x-axis; "iteration"
#' @param ylab name of y-axis; "logscore"
#' @param type type of line in the plot; "l" by default
#' @param col colour of line in the plot; "#0c2c84" by default
#' 
#' @rdname orderMCMC
#' @method plot orderMCMC    
#' @export
plot.orderMCMC <-function(x, ..., burnin = 0.2, main="DAG logscores", xlab="iteration", ylab="logscore", type="l", col="#0c2c84"){
  old.par <- par(no.readonly = TRUE)
  scorevec<-x$trace
  vecl<-length(scorevec)
  burnin<-ceiling(vecl*burnin)
  score20<-min(scorevec[burnin:vecl])
  scoremax<-max(scorevec)
  scoremin<-min(scorevec)
  
  par(mfrow=c(1,2))
  par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
  plot(scorevec,col=col, xlab=xlab, ylab=ylab, type=type, ...,
       ylim=c(scoremin,scoremax+(scoremax-scoremin)*0.02), main=main)
  par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
  plot(x=c(burnin:vecl),y=scorevec[burnin:vecl], type=type, col=col, xlab=xlab, ylab=ylab,
       ylim=c(score20,scoremax), main="excluding burn-in")
  par(old.par)
}



#' Plotting object of class 'partitionMCMC'
#'
#' @param x object of class 'partitionMCMC'
#' @param ... other parameters to be passed through to plotting functions
#' @param burnin number between \code{0} and \code{1}, indicates the percentage of the samples which will be discarded as `burn-in' of the MCMC chain; the rest  of the samples will be used to calculate the posterior probabilities; 0.2 by default
#' @param main name of the graph; "DAG logscores" by default
#' @param xlab name of x-axis; "iteration"
#' @param ylab name of y-axis; "logscore"
#' @param type type of line in the plot; "l" by default
#' @param col colour of line in the plot; "#0c2c84" by default
#' 
#' @rdname partitionMCMC
#' @method plot partitionMCMC
#' @export
plot.partitionMCMC <-function(x, ..., burnin = 0.2, main="DAG logscores", xlab="iteration", ylab="logscore", type="l", col="#0c2c84"){
  old.par <- par(no.readonly = TRUE)
  scorevec<-x$trace
  vecl<-length(scorevec)
  burnin<-ceiling(vecl*burnin)
  score20<-min(scorevec[burnin:vecl])
  scoremax<-max(scorevec)
  scoremin<-min(scorevec)
  
  par(mfrow=c(1,2))
  par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
  plot(scorevec,col=col, xlab=xlab, ylab=ylab, type=type, ...,
       ylim=c(scoremin,scoremax+(scoremax-scoremin)*0.02), main=main)
  par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
  plot(x=c(burnin:vecl),y=scorevec[burnin:vecl], type=type, col=col, xlab=xlab, ylab=ylab,
       ylim=c(score20,scoremax), main="excluding burn-in")
  par(old.par)
}


#' Plotting object of class 'itercomp'
#'
#' @param x object of class 'itercomp'
#' @param ... other parameters to be passed through to plotting functions
#' @param vars a tuple of variables which will be used for 'x' and 'y' axes; possible values: "SHD", "TP", "FP", "TPR", "FPR", "FPRn", "FDR", "score" 
#' @param type type of line in the plot;"b" by default
#' @param col colour of line in the plot; "blue" by default
#' @param showit (optional) vector of integers specifying indices of search expansion iterations to be labelled; by default no iterations are labelled
#' @rdname itercomp
#' @method plot itercomp
#' @export
plot.itercomp <-function(x, ..., vars = c("FP", "TP"), type="b", col="blue", showit=c()){
  old.par <- par(no.readonly = TRUE)
  nit<-nrow(x)
  if(nit>1) {
  scales<-vector()
  for(i in colnames(x)) {
    localrange<-max(x[,i])-min(x[,i])
    scales[i]<-localrange/nit*0.15
  }
  
  plot(x[,vars[1]], x[,vars[2]], type=type, col=col, xlab=vars[1],ylab=vars[2], ...)
  if(length(showit)>0) {
    for(i in showit) {
      text(x[i,vars[1]]+scales[vars[1]], x[i,vars[2]]-0.5*scales[vars[2]], i)
    }
  }
  }
  par(old.par)
}

#' Plotting object of class 'samplecomp'
#'
#' @param x object of class 'samplecomp'
#' @param ... other parameters to be passed through to plotting functions
#' @param vars a tuple of variables which will be used for 'x' and 'y' axes; possible values: "SHD", "TP", "FP", "TPR", "FPR", "FPRn", "FDR"
#' @param type type of line in the plot; "b" by default
#' @param col colour of line in the plotl; "blue" by default
#' @param showp logical, defines if points are labelled with the posterior threshold corresponding to the assessed model
#' 
#' @rdname samplecomp
#' @method plot samplecomp
#' @export
plot.samplecomp <-function(x, ..., vars = c("FP", "TP"), type="b", col="blue", showp=NULL){
  old.par <- par(no.readonly = TRUE)
  nit<-nrow(x)
  if(is.null(nit)) {
    message("plotting structure fit for the only threshold")
    plot(x[vars[1]], x[vars[2]], type="p", col=col, xlab=vars[1],ylab=vars[2], ...)
  } else if(nit>1) {
    scales<-vector()
    lims<-matrix(ncol=2,nrow=2)
    k<-1
    for(i in vars) {
      localrange<-max(x[,i])-min(x[,i])
      scales[i]<-localrange/nit*0.15
      lims[k,]<-c(min(x[,i])-localrange*0.15,max(x[,i])+localrange*0.15)
      k<-k+1
    }
    
    
    plot(x[,vars[1]], x[,vars[2]], type=type, col=col, xlab=vars[1],ylab=vars[2], ...,
         xlim=lims[1,],ylim=lims[2,])
    if(showp) {
      for(i in 1:nit) {
        text(x[i,vars[1]]+1.3*scales[vars[1]], x[i,vars[2]]-scales[vars[2]], x[i,"p"])
      }
    }
  }

  par(old.par)
}

