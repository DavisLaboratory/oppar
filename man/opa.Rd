% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opa-generic-methods.R
\docType{methods}
\name{opa}
\alias{opa}
\alias{opa,ExpressionSet-method}
\alias{opa,matrix-method}
\title{Outlier profile Analysis}
\usage{
opa(exprs.matrix, ...)

\S4method{opa}{matrix}(exprs.matrix, group, upper.quantile = 0.95,
  lower.quantile = 0.05)

\S4method{opa}{ExpressionSet}(exprs.matrix, group, upper.quantile = 0.95,
  lower.quantile = 0.05)
}
\arguments{
\item{exprs.matrix}{Gene expression data. Can be either a matrix or an
object of type ExpressionSet.}

\item{...}{Numeric. To supply values for upper.quantile and lower.quantile
arguments if default values are going to be override.}

\item{group}{A vector of factors representing the groups to which each sample belong.
This can be either a vector of 0s and 1s, or normal and cases.}

\item{upper.quantile}{Numeric. The cut-off for upper quantile when determining outliers.
Default to 0.95}

\item{lower.quantile}{Numeric. The cut-off for lower quantile when determining outliers.
Default to 0.05}
}
\value{
\code{opa} returns an object of type \code{OPPARList}. The outlier profiles
  are stored in \code{profileMatrix} and can be accessed using $. It it also
  possible to retrieve parameters used to run the outlier profile analysis, such
  as \code{upper.quantile}, \code{lower.quantile} via the $ operator.
}
\description{
Returns a matrix with 0, -1 and 1 entries that describe outlier profiles
in samples. The rows reperesent genes and the columns represent samples. -1 implies
that the gene is a down-regulated outlier, 1 indicates an up-regulate outlier and 0
 means that the gene is not an outlier in a sample.
}
\section{Methods (by class)}{
\itemize{
\item \code{matrix}: opa(exprs.matrix, group, lower.quantile = 0.05, upper.quantile = 0.95)

\item \code{ExpressionSet}: opa(eset, group, lower.quantile = 0.05, upper.quantile = 0.95)
}}
\examples{
# loading bcm object from GSE46141 dataset
data(GSE46141)
library(Biobase)
# defining the group variable. local breast tumors are the controls
# and the rest of the samples are the diseased samples
group <- sapply(pData(bcm)$source_name_ch1, function(x){ ifelse(x == "breast",0,1)})
group <- factor(group)
# running opa with default values (i.e upper.quantile = 0.95, lower.quantile = 0.05)
# the result is an object of type OPPARList
opa(bcm,group = group)
}
\seealso{
Wang, C., Taciroglu, A., Maetschke, S. R., Nelson, C. C., Ragan, M. A., & Davis, M. J. (2012).
    mCOPA: analysis of heterogeneous features in cancer expression data. Journal
    of Clinical Bioinformatics, 2, 22. http://doi.org/10.1186/2043-9113-2-22
}

