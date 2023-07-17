library(tm)
library(keyperm)

# load subcorpora "acq" and "crude" from Reuters

data(acq)
data(crude)

# convert to term-document-matrices and combine into single tdm

acq_tdm <- TermDocumentMatrix(acq, control = list(removePunctuation = TRUE))
crude_tdm <- TermDocumentMatrix(crude, control = list(removePunctuation = TRUE))
tdm <- c(acq_tdm, crude_tdm)

# generate a logical that indicates whether document comes from "acq" or "crude"

ndoc_A <- dim(acq_tdm)[2]
ndoc_B <- dim(crude_tdm)[2]
corpus <- rep(c(TRUE, FALSE), c(ndoc_A, ndoc_B))

# generate an indexed frequency list, the data structure used by keyperm

reuters_ifl <- create_ifl(tdm, corpus = corpus)

# calculate Log-Likelihood-Ratio scores for all terms and calculate
# p-values according to the (wrong) token-by-token sampling model

llr <- keyness_scores(reuters_ifl, type = "llr", laplace = 0)
head(round(pchisq(llr, df = 1, lower.tail = FALSE), digits = 4), n = 10)

# generate permutation distribution and p-values based on document-by-document sampling model

keyp <- keyperm(reuters_ifl, llr, type = "llr", 
                laplace = 0, output = "counts", nperm = 1000)
head(p_value(keyp, alternative = "greater"), n = 10)

# generate observed log-ratio values and (one-sided) p-values based
# on the permutation distribution (document-by-document sampling model)
# laplace-correction used (adding one occurence to both corpora)

logratio <- keyness_scores(reuters_ifl, type = "logratio", laplace = 1)
keyp2 <- keyperm(reuters_ifl, logratio, type = "logratio", 
                laplace = 1, output = "counts", nperm = 1000)
head(p_value(keyp2, alternative = "greater"), n = 10)

# it may be of interest to improve accuracy of the small p-values. 
# Think of this in terms of spending the computational budget mainly 
# on the terms for which higher accuracy matters most 

pvals <- p_value(keyp2, alternative = "greater")
table(pvals > 0.1)

small_p <- which(pvals < 0.1)

logratio_subset <- logratio[small_p]
reuters_ifl_subset <- create_ifl(tdm, subset_terms = small_p, corpus = corpus)

keyp2_subset <- keyperm(reuters_ifl_subset, logratio_subset, type = "logratio", 
                 laplace = 1, output = "counts", nperm = 9000)

# combine counts from both runs using the combiner

keyp2_combined <- combine_results(keyp2, keyp2_subset)

# smaller p-values are based on 1000, the larger ones on 10000 random permutations
# note that 10000 is still far too small for real applications

head(p_value(keyp2_combined, alternative = "greater"), n = 10)

