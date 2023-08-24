#' Convert results of permutation test for keyness to p-values
#' 
#' Calculate p-values from the results of \code{keyperm()} with 
#' \code{output = "counts"}. 
#' 
#' Valid (slightly conservative) p-values are calculated from an 
#' object of class \code{keyperm_results_counts} that is obtained  
#' by running \code{keyperm()} with \code{output = "counts"}. 
#' \code{keyperm_results_counts} is a matrix with three columns that 
#' contain the counts of generated permutations that resulted in a score
#' strictly less than, equal to and strictly greater that the observed score. 
#' 
#' For a one-sided p-value we use
#' 
#' \deqn{pvalue_greater = (no. greater + no. equal + 1)/(no. of perms + 1) }
#' or
#' \deqn{pvalue_less = (no. less + no. equal + 1)/(no. of perms + 1) }
#' 
#' Adding 1 in both the numerator and denominator amounts to including the observed
#' values. This results in a slightly conservative p-value, but guarantees that
#' the test is valid for any number of random permutations. It also means that
#' never a p-value of zero is returned but the minimum possible p-value is 
#' \eqn{1/(no. perms + 1)}. 
#' 
#' The two-sided p-value is calculated by
#' \deqn{pvalue_twosided = 2 * min(pvalue_less, pvalue_greater)} 
#' (values larger than 1 are set to 1). 
#' 
#' If \code{alternative} is not specified by the user, different defaults are
#' used depending on the scoretype (which is included as an attribute
#' in the \code{keyperm_results_counts} object). 
#' Since for \code{llr} and \code{chisq}, large values indicate a great 
#' deviation from equal frequencies without indicating the direction, 
#' \code{alternative == "greater"} is basically the only alternative of interest 
#' and is used as a default. 
#' For \code{diff} and \code{logratio} large absolute values indicate 
#' a great deviation from equal frequencies, and positive values correspond to 
#' higher frequencies in A, negative frequencies correspond to a higher frequency in B.
#' For these scoretypes, the default is \code{alternative = "two.sided"}. 
#' If only "positive" keywords for A with respect to B are desired, use \code{alternative = "less"}.
#' 
#' @param results results from permutation test.  
#'     Must be of class \code{keyperm_results_counts} 
#'     (obtained by setting \code{output = "counts"} in \code{keyperm()})
#' @param alternative direction of p-value to calculate, one of \code{"two.sided"},
#'     \code{"greater"}, \code{"less"}. Defaults depend on the scores used. See details. 
#' @return a numeric vector of p-values.
#' @export
p_value <- function(results, alternative = NULL) {
  if (!inherits(results, "keyperm_results_counts"))
    stop("only input of class keyperm_results_count is currently supported")
  if (is.null(alternative)) {
    alternative <- switch(attr(results, "scoretype"),
                          llr = "greater",
                          chisq = "greater",
                          diff = "two.sided",
                          logratio = "two.sided",
                          ratio = "two.sided")
  }  
  if (!(alternative %in% c("two.sided", "greater", "less")))
    stop("alternative must be one of \"two.sided\", \"greater\", or \"less\"")
  if (alternative == "greater") {
    pvals <- (results[,2] + results[,3] + 1) / (rowSums(results) + 1)
  } else if (alternative == "less") {
    pvals <- (results[,2] + results[,1] + 1) / (rowSums(results) + 1)
  } else {  # two-sided case: take correct tail and double p-value
    pvals_left <- (results[,2] + results[,1] + 1) / (rowSums(results) + 1)
    pvals_right <- (results[,2] + results[,3] + 1) / (rowSums(results) + 1)
    pvals <- 2 * pmin(pvals_left, pvals_right, 0.5)
  }  
  class(pvals) <- "keyperm_pvals"
  attr(pvals, "alternative") <- alternative
  pvals
}
