#' Tomlins et al. Prostate Cancer data (GEO: GSE6099)
#'
#'
#' An ExpressionSet object containing microarray gene expression measurements on normal tissue and metastatic
#' prostate cancer tumoures, and the corresponsing feature and phynotypic meta data.
#'
#'
#' @format An ExpressionSet object with 10945 features and 86 samples.
#' \describe{
#'   \item{title}{In pData(eset) -- Sample names}
#'   \item{geo_accession}{In pData(eset) --GEO accession numbers}
#'   \item{characteristics_ch1}{In pData(eset) -- sample description}
#'   \item{ID}{In fData(eset) -- probes Ids}
#'   \item{Gene.title}{...}
#'   \item{Gene.symbol}{In fData(eset) -- Gene Symbol}
#'   \item{Gene.ID}{In fData(eset) -- Entrez Gene IDs}
#' }
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6099}
"eset"
