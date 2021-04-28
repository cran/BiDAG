\name{iterativeMCMC class}
\alias{iterativeMCMC class}
\title{iterativeMCMC class structure}
\description{

  The structure of an object of S3 class \code{iterativeMCMC}.

}

\details{

  An object of class \code{iterativeMCMC} is a list containing at least the following
  components:

  \itemize{

    \item DAG: adjacency matrix of a maximum scoring DAG found/sampled in MCMC.
      
    \item CPDAG: adjacency matrix representing equivalence class of a maximum scoring DAG found/sampled in MCMC.
    
    \item score: score of a maximum scoring DAG found/sampled in MCMC.
    
    \item maxorder: order of a maximum scoring DAG found/sampled in MCMC.
    
    \item maxtrace: a list of graphs, their scores and orders used for search space expansion (either MAP graphs at each expansion iterations or consensus graphs corresponding to the threshold p)
    
    \item info: a list containing information about parameters and results of MCMC
    
    \item trace: a vector containing log-scores of sampled DAGs
    
    \item startspace: adjacency matrix representing the initial core space where MCMC was ran
    
    \item endspace: adjacency matrix representing the final core space where MCMC was ran

  Optional components:

  \itemize{

    \item \code{traceadd}: list which  consists of three elements:
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

