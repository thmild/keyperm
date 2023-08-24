#' @aliases keyperm-package
#' @details This package contains an implementation of the permutation testing approach
#' to keyness as used in corpus linguistics.
#'
#' Keywords are words that appear more frequently in one corpus compared to 
#' another corpus. Usually this is assessed using test statistics, for example
#' the likelihood-ratio test on 2x2 contingency tables, resulting in scores for every term
#' that appears in the document. 
#'   
#' Conventionally, keyness scores are judged by reference to a limiting null distribution
#' under a token-by-token-sampling model. \code{keyperm} approximates the null distribution under
#' a document-by-document sampling model. 
#' 
#' The permutation distributions of a given keyness measure for 
#' each term is calculated by repeatedly shuffeling the copus labels. 
#' Number of documents per corpus is kept constant.
#' 
#' Apart from obtaining null distributions of common test statistics like 
#' LLR and Chi-Square, \code{keyperm} can also obtain null distributions of
#' of the logratio measure that is normally used as an effect size. 
#' 
#' Currently, the following types of scores are supported:
#' \describe{
#'     \item{\code{llr}}{The log-likelihood ratio}
#'     \item{\code{chisq}}{The Chi-Square-Statistic}
#'     \item{\code{diff}}{Difference of relative frequencies}
#'     \item{\code{logratio}}{Binary logarithm of the ratio of the relative frequencies, possibly using a laplace correction to avoid infinite values.}
#' }
#' 
#' The actual resampling procedure is implemented in an efficient manner using 
#' the Rcpp package and utilizing a special data structure (indexed frequency list). 
#' Currently, \code{keyperm} can generate indexed frequency list from term-document-matrices as implemented in
#' the package \code{tm}.
#' 
#' @keywords internal
#' 
#'
#' @example demo/example_reuters.R
"_PACKAGE"

#' @importFrom tm TermDocumentMatrix


