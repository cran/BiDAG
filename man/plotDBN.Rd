% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotusingrgraphviz.R
\name{plotDBN}
\alias{plotDBN}
\title{Plotting a DBN}
\usage{
plotDBN(DBN, struct = c("init", "trans"), b = 0, shape = "circle", ...)
}
\arguments{
\item{DBN}{binary matrix (or a graph object) representing a 2-step DBN (compact or unrolled)}

\item{struct}{option used to determine if the initial or the transition structure should be plotted; acceptable values are init or trans}

\item{b}{number of static variables in the DBN, 0 by default; note that for function to work correctly all static variables have to be in the first b columns of the matrix}

\item{shape}{string, defining the shape of the box around each node; possible values are circle, ellipse, box}

\item{...}{optional parameters passed to \code{Rgraphviz} plotting functions e.g. \code{main}, \code{fontsize}}
}
\value{
plots the DBN defined by the adjacency matrix 'DBN' and number of static and dynamic variables. When 'struct' equals "trans" the transition structure is plotted,
otherwise initial structure is plotted
}
\description{
This function can be used for plotting initial and transition structures of a dynamic Bayesian network.
}
\examples{
plotDBN(DBNmat, "init", b=3)
plotDBN(DBNmat, "trans", b=3)

}
\author{
Polina Suter
}
