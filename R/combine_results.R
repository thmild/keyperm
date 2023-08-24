#' Combine results of permutation test for keyness 
#' 
#' Combine results of two runs of \code{keyperm()} with 
#' \code{output = "counts"}, possibly with different subsets of terms.
#' 
#' Results of two runs of \code{keyperm()} with \code{output = "counts"}, i.e. objects of
#' type \code{keyperm_results_counts} using can be combined 
#' using \code{combine_results()}. For this to make sense, \code{scoretype} needs to be 
#' the same in both results, but terms in both objects need not be the same.
#' 
#' There are at least two important uses of the function:
#' 
#' Parallelization: \code{keyperm()} is run several times with the same parameters 
#' on different cores, using \code{parallel::mclapply()} or a similar function.
#' 
#' Screening runs: \code{keyperm()} is first run using a small to medium number of permutations,
#' but considering all terms. Terms with p-values clearly exceeding some reasonable 
#' significance threshold are then excluded, and \code{keyperm()} is run a second time with a 
#' (preferably) large number of permutations but using only the remaining terms. The results of
#' both runs can then be combined into one object. The rationale behind this approach is that 
#' in many cases small p-values need to be determined with much greater accurary than larger ones
#' far away from significance, especially if a correction for multiple testing is to be applied 
#' or the p-values are used for ranking (although they should not...). 
#' 
#' @param results_1 Results from permutation test.  
#'     Must be of class \code{keyperm_results_counts} 
#'     (obtained by setting \code{output = "counts"} in \code{keyperm()})
#' @param results_2 Results from permutation test.  
#'     Must be of class \code{keyperm_results_counts} and have the same scoretype as
#'     \code{results_1}.  
#' @return An object of class \code{keyperm_results_counts}
#' 
#' @export
combine_results <- function(results_1, results_2) {
   if (!inherits(results_1, "keyperm_results_counts") || !inherits(results_2, "keyperm_results_counts"))
     stop("only input of class keyperm_results_count is currently supported")
    
   scoretype <- attr(results_1, "scoretype")
    if (attr(results_2, "scoretype") != scoretype)
      stop("scoretypes of results_1 and results_2 are not identical")
    
   new_terms <- union(dimnames(results_1)[[1]], dimnames(results_2)[[1]])
   not_in_1 <- setdiff(dimnames(results_2)[[1]], dimnames(results_1)[[1]])
   not_in_2 <- setdiff(dimnames(results_1)[[1]], dimnames(results_2)[[1]])
   pad_1 <- matrix(0, nrow = length(not_in_1), ncol = 3)
   dimnames(pad_1) <- list(not_in_1, c("less", "equal", "greater"))
   results_1 <- rbind(results_1, pad_1)
   results_1 <- results_1[new_terms,]
    
   pad_2 <- matrix(0, nrow = length(not_in_2), ncol = 3)
   dimnames(pad_2) <- list(not_in_2, c("less", "equal", "greater"))
   results_2 <- rbind(results_2, pad_2)
   results_2 <- results_2[new_terms,]

   result_combined <- results_1 + results_2
   class(result_combined) <- "keyperm_results_counts"  
   attr(result_combined, "scoretype") <- scoretype
   attr(result_combined, "output") <- "counts"
   
   result_combined
    
  }