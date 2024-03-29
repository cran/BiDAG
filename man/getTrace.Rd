% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractors.R
\name{getTrace}
\alias{getTrace}
\title{Extracting trace from MCMC object}
\usage{
getTrace(x, which = 0)
}
\arguments{
\item{x}{object of class 'orderMCMC','partitionMCMC' or 'iterativeMCMC'}

\item{which}{integer, indication which trace is returned: DAG scores (which = 0), DAGs (which = 1),
orders (which = 2), order scores (which = 3)}
}
\value{
a list or a vector of objects representing MCMC trace, depends on parameter 'which'; by default, the trace of DAG scores is returned
}
\description{
This function extracts a trace of
\itemize{
\item DAG scores 
\item DAG adjacency matrices
\item orders
\item order scores 
}
from the result of the MCMC run. Note that the last three options
work only when the parameter 'scoreout' was set to TRUE.
}
\examples{
myscore<-scoreparameters("bge",Boston)
\dontrun{
orderfit<-sampleBN(myscore,algorithm="order")
DAGscores<-getTrace(orderfit,which=0)
DAGtrace<-getTrace(orderfit,which=1)
orderscores<-getTrace(orderfit,which=3)
}
}
