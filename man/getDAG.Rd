% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractors.R
\name{getDAG}
\alias{getDAG}
\title{Extracting adjacency matrix (DAG) from MCMC object}
\usage{
getDAG(x, amat = TRUE, cp = FALSE)
}
\arguments{
\item{x}{object of class 'orderMCMC','partitionMCMC' or 'iterativeMCMC'}

\item{amat}{logical, when TRUE adjacency matrix is returned and object of class 'graphNEL' otherwise}

\item{cp}{logical, when TRUE the CPDAG (equivalence class) is returned and DAG otherwise; FALSE by default}
}
\value{
adjacency matrix of a maximum scoring DAG (or CPDAG) discovered/sampled in one MCMC run
}
\description{
This function extracts an adjacency matrix of
a maximum scoring DAG from the result of the MCMC run.
}
\examples{
myscore<-scoreparameters("bge", Boston)
\dontrun{
itfit<-learnBN(myscore,algorithm="orderIter")
maxEC<-getDAG(itfit,cp=TRUE)
}
}
