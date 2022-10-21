# keyperm: Keyword Analysis Using Permutation Tests

Implementation of permutation-based keyword analysis for corpus linguistics.

We propose a new approach for assessing keyness in corpus linguistics. Traditional approaches based on hypothesis tests (e.g. Likelihood Ratio Test) model the copora as independent identically distributed samples of **tokens**. This model does not account for the often observed uneven distribution of occurrences of a word across a corpus. When occurrences of a word are concentrated in few documents, large values of LLR and similar scores are in fact much more likely than accounted for by the token-by-token sampling model, leading to false positives.  

We replace the **token-by-token sampling model** by a model where corpora are *samples of documents* rather than tokens, which is much closer to the way corpora are actually assembled. We then use a permutation approach to approximate the distribution of a given keyness score under the null hypothesis of equal frequencies and obtain p-values for assessing significance. We do not need any assumption on how the tokens are organized within or across documents, and the approach works with basically **any** keyness score. The package currently implements **LLR**, **Chi-Square** and **Logratio scores**, the latter with a Laplace-type modification to avoid dividing by zero. 

The permutations test has been implemented in C++ using Rcpp, and the data are stored in a format that allows for fast calculations. 

The package is not yet on CRAN but can be installed by

``` 
library(devtools)
install_github("thmild/keyperm")
```

The package includes C++ code that needs to be compiled during installation, so compilers etc. are needed. Especially on Windows, these might be tedious to install. 

The main function ist `keyperm()`, see documentation. Before running `keyperm()`, the frequency counts have to be converted to a special format we call an *indexed frequency list* by calling `create_ifl()`. Currently, only term-document matrices in form of `tdm` objects from package `tm` are supported. Other input formats may be included in later versions. 

```
demo("example_reuters", package = "keyperm")
```

runs a commented toy example. 

A paper describing the methodology in detail is currently in preparation. 

