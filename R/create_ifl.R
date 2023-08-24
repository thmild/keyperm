
#' Create an Indexed Frequency List 
#'
#' The keyperm package stores frequency lists in a special data structure called
#' indexed frequency list. This can currently be created from a tdm object as   
#' implemented in the tm package.
#'    
#' Indexed frequency lists are essentially frequency lists stored in a three-column format,
#' similar to the simple triplet matrix internally used by tm to store term-document-matrices.
#' The first column stores number of document \code{i}, second number of term \code{j} and the third the 
#' frequencies with which the term \code{j} occurs in document \code{i}. Zero occurences are omitted. 
#' All columns contain integers, and the frequency list is sorted by document.  
#'                   
#' The object returned is of class \code{indexed_frequency_list}. In addition to the actual frequency 
#' list it contains an index for fast access as well as pre-computed total number of tokens per
#' document and total occurences per term. 
#' 
#' @param tdm a tdm-matrix from the tm package. Currently, this is the only supported input, but others may be added in later versions.
#' @param subset_terms vector of terms to be considered. Can be integer (indices) or boolean. Terms not included still are counted for total number of token per document.
#' @param subset_docs vector of documents to be considered. Can be integer (indices) or boolean. Documents excluded do not contribute to total number of occurences of a term. 
#' @param corpus vector indicating which documents belong to corpus A (first corpus). Can be integer (indices) or boolean. Currently, only comparisons of two corpora are supported.
#' @return A list with class \code{indexed_frequency_list} containing the following components:
#' @export
create_ifl <- function(tdm, 
                       subset_terms = 1:dim(tdm)[1],
                       subset_docs = 1:dim(tdm)[2],
                       corpus) {

  # To suppress CRAN warning
  
  tm::TermDocumentMatrix(NULL)
  
  # save row and column totals so that they do not have to be recalculated 
  #
  # note that we include all terms in the sums as we want total number of tokens per document
  
  tdm_rs <- slam::row_sums(tdm[subset_terms, subset_docs]) 
  tdm_cs <- slam::col_sums(tdm[, subset_docs]) 
  
  ntotal <- sum(tdm_cs)
  
  if (any(tdm_rs == 0)) warning("Total frequency of at least one term is zero!")
  if (any(tdm_cs == 0)) warning("At least one document is empty!")
  
  # subsetting
  # to do: take care about duplicated indices
  
  tdm <- tdm[subset_terms, subset_docs]  
  
  # convert to data frame 
  
  tdm_df <- data.frame(term = tdm$i, doc = tdm$j, freq = tdm$v)
  tdm_df <- tdm_df[order(tdm_df$doc),] # sort by document
  
  doc_lookup <- stats::aggregate(term ~ doc, tdm_df, length)
  names(doc_lookup)[2] <- "nterms"
  
  # this is necessary so that we can handle empty documents 
  
  doc_lookup <- merge(x = data.frame(doc = 1:dim(tdm)[2]), 
                      y = doc_lookup,
                      by = "doc",
                      all = TRUE,
                      sort = TRUE)
  
  doc_lookup$nterms[is.na(doc_lookup$nterms)] <- 0
  
  doc_lookup$start <- cumsum(doc_lookup$nterms) - doc_lookup$nterms + 1 # for R indexing convention
  
  # if corpus is a logical vectors, TRUE denotes corpus A 
  # we convert to indices
  #
  # later versions may support a factor as input 
  # 
  # we need to be careful with to select only components from corpus that are
  # included in subset_docs 
  
  if (is.logical(corpus)) {
    corpus <- which(corpus[subset_docs]) # easy case
  } else if (is.integer(corpus)) {
    corpus2 <- rep(FALSE, length(corpus))
    corpus2[corpus] <- TRUE
    corpus <- which(corpus2[subset_docs]) 
  } else {
    stop("only logical or integer vectors accepted for corpus")
  }
  
  # create the actual indexed_frequency_list
  # details may change in future versions
  
  ifl <- list(freqlist = tdm_df,
              index = doc_lookup,
              rowsums = tdm_rs,
              colsums = tdm_cs,
              ntotal = ntotal,
              terms = dimnames(tdm)$Terms,
              docs = dimnames(tdm)$Docs,
              corp_A = corpus
  )    
  class(ifl) <- "indexed_frequency_list"
  ifl
}
