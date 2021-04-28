\name{orderMCMC class}
\alias{orderMCMC class}
\title{orderMCMC class structure}
\description{

  The structure of an object of S3 class \code{orderMCMC}.

}
\details{

  An object of class \code{orderMCMC} is a list containing at least the following
  components:

  \itemize{

    \item DAG: adjacency matrix of a maximum scoring DAG found/sampled in the MCMC scheme.
      
    \item CPDAG: adjacency matrix representing equivalence class of a maximum scoring DAG found/sampled in MCMC.
    
    \item score: score of a maximum scoring DAG found/sampled in MCMC.
    
    \item maxorder: order of a maximum scoring DAG found/sampled in MCMC.
    
    \item info: a list containing information about parameters and results of MCMC.
    
    \item trace: a vector containing log-scores of sampled DAGs.

  Optional components:

  \itemize{

    \item \code{traceadd}: list which  consists of three or four elements (depending on MCMC scheme used for sampling):
        \itemize{
        \item incidence: list containg adjacency matrices of sampled DAGs
        \item order: list of orders from which the DAGs were sampled
        \item orderscores: order log-scores 
        }
    \item \code{scoretable}: object of class \code{\link{scorespace class}}

  }
}
}
\author{Polina Suter}
\keyword{classes}

