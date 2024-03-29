% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R, R/plotS3.R, R/printmethods.R,
%   R/summary.R
\name{orderMCMC}
\alias{orderMCMC}
\alias{plot.orderMCMC}
\alias{print.orderMCMC}
\alias{summary.orderMCMC}
\title{Structure learning with the order MCMC algorithm}
\usage{
orderMCMC(
  scorepar,
  MAP = TRUE,
  plus1 = TRUE,
  chainout = FALSE,
  scoreout = FALSE,
  moveprobs = NULL,
  iterations = NULL,
  stepsave = NULL,
  alpha = 0.05,
  cpdag = FALSE,
  gamma = 1,
  hardlimit = ifelse(plus1, 14, 20),
  verbose = FALSE,
  compress = TRUE,
  startspace = NULL,
  blacklist = NULL,
  startorder = NULL,
  scoretable = NULL
)

\method{plot}{orderMCMC}(
  x,
  ...,
  burnin = 0.2,
  main = "DAG logscores",
  xlab = "iteration",
  ylab = "logscore",
  type = "l",
  col = "#0c2c84"
)

\method{print}{orderMCMC}(x, ...)

\method{summary}{orderMCMC}(object, ...)
}
\arguments{
\item{scorepar}{an object of class \code{scoreparameters}, containing the data and score parameters, see constructor function \code{\link{scoreparameters}}}

\item{MAP}{logical, if TRUE (default) the search targets the MAP DAG (a DAG with maximum score),
if FALSE at each MCMC step a DAG is sampled from the order proportionally to its score}

\item{plus1}{logical, if TRUE (default) the search is performed on the extended search space}

\item{chainout}{logical, if TRUE the saved MCMC steps are returned, TRUE by default}

\item{scoreout}{logical, if TRUE the search space and score tables are returned, FALSE by default}

\item{moveprobs}{a numerical vector of 4 values in \code{\{0,1\}} corresponding to the probabilities of the following MCMC moves in the order space
\itemize{
\item exchanging 2 random nodes in the order
\item exchanging 2 adjacent nodes in the order
\item placing a single node elsewhere in the order
\item staying still
}}

\item{iterations}{integer, the number of MCMC steps, the default value is \eqn{6n^{2}\log{n}}}

\item{stepsave}{integer, thinning interval for the MCMC chain, indicating the number of steps between two output iterations, the default is \code{iterations/1000}}

\item{alpha}{numerical significance value in \code{\{0,1\}} for the conditional independence tests at the PC algorithm stage}

\item{cpdag}{logical, if TRUE the CPDAG returned by the PC algorithm will be used as the search
space, if FALSE (default) the full undirected skeleton will be used as the search space}

\item{gamma}{tuning parameter which transforms the score by raising it to this power, 1 by default}

\item{hardlimit}{integer, limit on the size of parent sets in the search space; by default 14 when MAP=TRUE and 20 when MAP=FALSE}

\item{verbose}{logical, if TRUE messages about the algorithm's progress will be printed, FALSE by default}

\item{compress}{logical, if TRUE adjacency matrices representing sampled graphs will be stored as a sparse Matrix (recommended); TRUE by default}

\item{startspace}{(optional) a square matrix, of dimensions equal to the number of nodes, which defines the search space for the order MCMC in the form of an adjacency matrix. If NULL, the skeleton obtained from the PC-algorithm will be used. If \code{startspace[i,j]} equals to 1 (0) it means that the edge from node \code{i} to node \code{j} is included (excluded) from the search space. To include an edge in both directions, both \code{startspace[i,j]} and \code{startspace[j,i]} should be 1.}

\item{blacklist}{(optional) a square matrix, of dimensions equal to the number of nodes, which defines edges to exclude from the search space. If \code{blacklist[i,j]} equals to 1 it means that the edge from node \code{i} to node \code{j} is excluded from the search space.}

\item{startorder}{(optional) integer vector of length n, which will be used as the starting order in the MCMC algorithm, the default order is random}

\item{scoretable}{(optional) object of class \code{scorespace} containing list of score tables calculated for example by the last iteration of the function \code{iterativeMCMC}. When not NULL, parameter \code{startspace} is ignored.}

\item{x}{object of class 'orderMCMC'}

\item{...}{ignored}

\item{burnin}{number between \code{0} and \code{1}, indicates the percentage of the samples which will be discarded as `burn-in' of the MCMC chain; the rest  of the samples will be used to calculate the posterior probabilities; 0.2 by default}

\item{main}{name of the graph; "DAG logscores" by default}

\item{xlab}{name of x-axis; "iteration"}

\item{ylab}{name of y-axis; "logscore"}

\item{type}{type of line in the plot; "l" by default}

\item{col}{colour of line in the plot; "#0c2c84" by default}

\item{object}{object of class 'orderMCMC'}
}
\value{
Object of class \code{orderMCMC}, which contains log-score trace of sampled DAGs as well 
as adjacency matrix of the maximum scoring DAG, its score and the order score. The output can optionally include DAGs sampled in MCMC iterations and the score tables. 
Optional output is regulated by the parameters \code{chainout} and \code{scoreout}. See \code{\link{orderMCMC class}} for a detailed class structure.
}
\description{
This function implements the order MCMC algorithm for the structure learning of Bayesian networks. This function can be used
for MAP discovery and for sampling from the posterior distribution of DAGs given the data.
Due to the superexponential size of the search space as the number of nodes increases, the 
MCMC search is performed on a reduced search space.
By default the search space is limited to the skeleton found through the PC algorithm by means of conditional independence tests 
(using the functions \code{\link[pcalg]{skeleton}} and \code{\link[pcalg]{pc}} from the `pcalg' package [Kalisch et al, 2012]).
It is also possible to define an arbitrary search space by inputting an adjacency matrix, for example estimated by partial correlations or other network algorithms.
Also implemented is the possibility to expand the default or input search space, by allowing each node in the network to have one additional parent.  This offers improvements in the learning and sampling of Bayesian networks.
}
\note{
see also extractor functions \code{\link{getDAG}}, \code{\link{getTrace}}, \code{\link{getSpace}}, \code{\link{getMCMCscore}}.
}
\examples{
\dontrun{
#find a MAP DAG with search space defined by PC and plus1 neighbourhood
Bostonscore<-scoreparameters("bge",Boston)
#estimate MAP DAG
orderMAPfit<-orderMCMC(Bostonscore)
summary(orderMAPfit)
#sample DAGs from the posterior distribution
ordersamplefit<-orderMCMC(Bostonscore,MAP=FALSE,chainout=TRUE)
plot(ordersamplefit)
}
}
\references{
P. Suter, J. Kuipers, G. Moffa, N.Beerenwinkel (2023) <doi:10.18637/jss.v105.i09>

Friedman N and Koller D (2003). A Bayesian approach to structure discovery in bayesian networks. Machine Learning 50, 95-125.

Kalisch M, Maechler M, Colombo D, Maathuis M and Buehlmann P (2012). Causal inference using graphical models with the R package pcalg. Journal of Statistical Software 47, 1-26.

Geiger D and Heckerman D (2002). Parameter priors for directed acyclic graphical models and the characterization of several probability distributions. The Annals of Statistics 30, 1412-1440.

Kuipers J, Moffa G and Heckerman D (2014). Addendum on the scoring of Gaussian acyclic graphical models. The Annals of Statistics 42, 1689-1691.

Spirtes P, Glymour C and Scheines R (2000). Causation, Prediction, and Search, 2nd edition. The MIT Press.
}
\author{
Polina Suter, Jack Kuipers, the code partly derived from the order MCMC implementation from Kuipers J, Moffa G (2017) <doi:10.1080/01621459.2015.1133426>
}
