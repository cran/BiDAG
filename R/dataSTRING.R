#' intSTRING dataset
#'
#' A data frame containing possible interactions between genes from \code{kirp} and \code{kirc}
#' data sets
#'
#' @source \url{https://string-db.org/}
#' @format A data frame with 179 rows and 2 columns; each row represent a possible interaction between two genes
#'
"intSTRING"

#' mapSTRING dataset
#'
#' A data frame containing mapping between names of genes used in \code{kirp}/\code{kirc}
#' data sets and names used in STRING interactions list (see \code{\link{intSTRING}}).
#'
#' @source \url{https://string-db.org/}
#' @format A data frame with 46 rows and two columns:
#' \itemize{
#' \item queryItem character, name used for structure learning (for example see column names of \code{\link{kirp}})
#' \item preferredName character, name used in STRING interactions data set (see \code{\link{intSTRING}})
#'}
#'
"mapSTRING"
