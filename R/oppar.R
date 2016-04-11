#' oppar: A package for outlier profile and pathway analysis.
#'
#' The oppar package provides 3 main function for outlier profile analysis:
#' opa, getSampleOutlier and getSubtypeProbes, and 1 function for pathway and
#' gene set enrichment analysis, based on gsva function implemented in the GSVA package
#'
#' @section opa:
#' calculates the outlier profile matrix, using the method proposed in Wang et al. (2012) paper
#' @section getSampleOutlier:
#' extracts outlier profile for individual samples
#' @section getSubtypeProbes:
#' extracts outlier profile for a group of related samples, such as subtypes
#'
#' @import methods
#'
#' @docType package
#' @name oppar
NULL

#' Breast cancer metastases from different anatomical sites
#'
#' An ExpressionSet object containing trimmed GSE46141 data. The object contains
#' gene expression measurements on local breast tumours and liver, lymph node, skin local-
#' regional etc metastatic tumours. Contains 9503 features and 80 samples (ascite, bone, lung and skin
#' samples were removed).
#' @format Contains gene expression matrix, phenotype (pData) and feature (fData) data:
#' \describe{
#'   \item{ID}{In fData(e) -- probe IDs}
#'   \item{EntrezGeneID}{In fData(e) -- Entrez Ids}
#'   \item{GeneSymbol}{In fData(e) -- Gene Symbols}
#'   \item{geo_accession}{In pData(e) -- GEO accession IDs}
#'   \item{source_name_ch1}{In pData(e) -- Sample information}
#' }
#' @docType data
#' @name bcm
#' @usage data(GSE46141)
#' @keywords datasets
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46141}
#' @return An ExpressionSet object
NULL

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
#' @docType data
#' @name eset
#' @usage data(tomlins)
#' @keywords datasets
#' @return An ExpressionSet object
NULL

#' Maupin's TGFb data and a TGFb gene signature
#'
#' A list consisting of two components:
#' data: A matrix consisting of gene expression values on 3 control
#' and 3 TGFb treated samples.
#' sig: A TGFb gene signature. A dataframe containing a gene signature
#' and information on whether genes in the signature are up or dow regulated.
#'
#' @format A list of two components:
#' \describe{
#'   \item{M_Ctrl_R1}{In data; Control Sample - Replicate 1}
#'   \item{M_Ctrl_R2}{In data; Control Sample - Replicate 2}
#'   \item{M_Ctrl_R3}{In data; Control Sample - Replicate 3}
#'   \item{M_TGFb_R1}{In data; Control Sample - Replicate 1}
#'   \item{M_TGFb_R2}{In data; Control Sample - Replicate 2}
#'   \item{M_TGFb_R3}{In data; Control Sample - Replicate 3}
#'   \item{EntrezID}{In sig; The Entrez IDs}
#'   \item{Symbol}{In sig; Gene symbols}
#'   \item{upDown}{In sig; Direction of regulation e.g. up, down}
#'
#'   }
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23952}
#' @docType data
#' @name maupin
#' @usage data(Maupin)
#' @keywords datasets
#' @return A list of 2
NULL

