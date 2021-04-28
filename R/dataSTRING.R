#' interactions dataset
#'
#' A data frame containing possible interactions between genes from \code{kirp} and \code{kirc}
#' data sets
#'
#' @source \url{https://string-db.org/}
#' @format A data frame with 179 rows and 3 columns; 
#' \itemize{
#' \item node1 character, name of a gene 
#' \item node2 character, name of a gene 
#' \item combined_score interaction score, reflecting confidence in the fact that interaction between gene1 and gene2 is possible
#'}
#' each row represents a possible interaction between two genes
#'
"interactions"

#' mapping dataset
#'
#' A data frame containing mapping between names of genes used in \code{kirp}/\code{kirc}
#' data sets and names used in STRING interactions list (see \code{\link{interactions}}).
#'
#' @source \url{https://string-db.org/}
#' @format A data frame with 46 rows and two columns:
#' \itemize{
#' \item queryItem character, name used for structure learning 
#' \item preferredName character, name used in STRING interactions data set
#'}
#'
"mapping"
