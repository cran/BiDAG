\name{MCMCscoretab}
\alias{MCMCscoretab}
\title{MCMCscoretab class structure}
\description{

  The structure of an object of S3 class \code{MCMCscoretab}.

}
\details{

  An object of class \code{MCMCscoretab} is a list containing at least the following
  components:

  \itemize{

    \item adjacency: adjacency martrix representing the core search space
    
    \item blacklist: adjacency martrix representing the blacklist used for computing score tables
      
    \item tables: a list of matrices (for core search space) or a list of lists of matrices (for extended search space) containing quantities needed for scoring orders and sampling DAGs in MCMC schemes; this list corresponds to adjacency and blacklist
    
    
}
}
\author{Polina Suter}
\keyword{classes}