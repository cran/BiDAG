#' KIPAN dataset
#'
#' Mutation data from TCGA pan-kidney cohort (KIPAN). 
#' Mutations are picked according to q-value computed 
#' by MutSigCV. Only mutations with q-value of less than
#' 0.05 for any of kidney cancer subtypes (KIRC, KIRP, KICH) are included. 
#'
# @source \url{https://portal.gdc.cancer.gov/}
#' @details { 
#' Each variable represents a gene. If in sample i gene j contains a mutation, than j-th element in 
#' row i equals 1, and 0 otherwise.
#' The rows are named according to sample names in TCGA.
#' The columns are named according to gene symbols.
#' }
#'@references \url{https://portal.gdc.cancer.gov/}
#'@references \url{http://firebrowse.org/iCoMut/?cohort=kich}
#'@references Lawrence, M. et al. Mutational heterogeneity in cancer and the search for new cancer-associated genes. Nature 499, 214-218 (2013). 
#'
"kipan"