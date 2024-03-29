% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractors.R
\name{getMCMCscore}
\alias{getMCMCscore}
\title{Extracting score from MCMC object}
\usage{
getMCMCscore(x)
}
\arguments{
\item{x}{object of class 'orderMCMC','partitionMCMC' or 'iterativeMCMC'}
}
\value{
a score of a maximum-scoring DAG found/sampled in one MCMC run
}
\description{
This function extracts the score of a maximum DAG sampled in the MCMC run.
}
\examples{
myscore<-scoreparameters("bge", Boston)
\dontrun{
itfit<-learnBN(myscore,algorithm="orderIter")
getMCMCscore(itfit)
}
}
