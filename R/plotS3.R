#' S3 methods for class 'MCMCmult'
#' 
#' Plotting object of class 'MCMCmult'
#'
#' @param x object of class 'MCMCmult'
#' @param ... ignored
#' 
# @rdname MCMCmult
#' @method plot MCMCmult
#' @export plot.MCMCmult
#' @export
plot.MCMCmult <-function(x,...){
  old.par <- par(no.readonly = TRUE)
  x<-x$trace
  nchains<-length(x)
  nsteps<-length(x[[1]])
  scorevec<-Reduce('c',x)
  
  scoremax<-max(scorevec)
  scoremin<-min(scorevec)
  # Add extra space to right of plot area; change clipping to figure
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(scorevec,type="l",xlab="MCMC step",ylab="DAG logscores",
       ylim=c(scoremin,scoremax),main="iterative MCMC, DAG scores",cex.main=1.2,
       col="blue")
  legend(length(scorevec)*1.04, scoremin+(scoremax-scoremin)*0.2,col="grey", lty=2,angle=90, 
         legend=c("search \nspace \nexpansion"),bty="n")
  
  par(xpd=FALSE)
  for(i in 1:(nchains-1)) {
    abline(v=i*nsteps,col="grey",lty=2,ylim=c(scoremin,scoremax))
  }
  par(xpd=TRUE)
  par(old.par)
}


#' S3 methods for class 'MCMCres'
#' 
#' Plotting object of class 'MCMCres'
#'
#' @param x object of class 'MCMCres'
#' @param ... ignored
#' 
# @rdname MCMCres
#' @method plot MCMCres
#' @export plot.MCMCres
#' @export
plot.MCMCres <-function(x,...){
  old.par <- par(no.readonly = TRUE)
  scorevec<-x$trace
  burnin<-0.2
  vecl<-length(scorevec)
  burnin<-ceiling(vecl*burnin)
  score20<-min(scorevec[burnin:vecl])
  scoremax<-max(scorevec)
  scoremin<-min(scorevec)
  
  par(mfrow=c(1,2))
  par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
  plot(scorevec,type="l",col="#0c2c84",xlab="iteration",ylab="logscore",
       ylim=c(scoremin,scoremax+(scoremax-scoremin)*0.02),main="DAG logscores",cex.main=1)
  par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
  plot(x=c(burnin:vecl),y=scorevec[burnin:vecl],type="l",col="#0c2c84",xlab="iteration",ylab="logscore",
  ylim=c(score20,scoremax),main="excluding burn-in (20% of samples)",cex.main=1)
  par(old.par)
}



