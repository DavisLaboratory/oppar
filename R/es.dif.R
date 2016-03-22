#' The cached results of the modified gsva output
#'
#' A matrix on one row and six columns. The row represents the gene signature
#' consisting of an up component and a down component. The columns are enrichment
#' scores for each samples in Maupin's data (see \code{?maupin})
#'
#' @format numeric matrix
#' \describe{
#'   \item{M_Ctrl_R1}{Enrichment score (es) for Control Sample - Replicate 1}
#'   \item{M_Ctrl_R2}{es for Control Sample - Replicate 2}
#'   \item{M_Ctrl_R3}{es for Control Sample - Replicate 3}
#'   \item{M_TGFb_R1}{es for Control Sample - Replicate 1}
#'   \item{M_TGFb_R2}{es for Control Sample - Replicate 2}
#'   \item{M_TGFb_R3}{es for Control Sample - Replicate 3}
#'   }
#'
"es.dif"
