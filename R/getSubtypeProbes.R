#' @include class_OPPARList.R class_nuchar.R
NULL

#' Retrieving outlier genes from a group of related samples
#'
#' Returns a list of genes that are outlier in a group of samples, such as samples from the same subtype.
#'
#' @param profileMatrix  A matrix of 0,1 and -1, representing outlier genes in samples.
#'    Also an object of type \code{OPPARList}.
#' @param sample.names  A character vector containing sample names, or a numeric vector
#'   containing the indices of the samples.
#' @return A list of lists. The sub-lists are up-regulated outlier genes, and down-regulated outlier genes.
#' @export
#' @examples
#' data(GSE46141g)
#' eset <- g[[1]]
#' library(Biobase)
#' group <- sapply(pData(eset)$source_name_ch1, function(x){ ifelse(x == "breast",0,1)})
#' group <- factor(group)
#' bcm.opa <- opa(eset,group=group)
#' # extracting liver samples
#' index <- which(pData(eset)$source_name_ch1 == "liver")
#' samples <- rownames(pData(eset)[index,])
#' samples <- match(samples, colnames(bcm.opa$profileMatrix))
#' samples <- Reduce(c,samples)
#' # liver subtype outlier profile
#' liver.subtype.outlier <- getSubtypeProbes(bcm.opa, samples)
setGeneric("getSubtypeProbes",
					 function(profileMatrix,sample.names){ standardGeneric("getSubtypeProbes")})



#' @describeIn getSubtypeProbes A method for getSubtypeProbes with signature profileMatrix = \code{matrix}
#'   and sample.names = \code{nuchar}
#' @export
setMethod("getSubtypeProbes",
					signature(profileMatrix = "matrix", sample.names = "nuchar"),
					function(profileMatrix, sample.names){
						up.outlr <- getOutlier(profileMatrix, sample.names)(1)
						down.outlr <- getOutlier(profileMatrix, sample.names)(-1)
						list( up = unique(unlist(up.outlr, use.names = FALSE)),
									down = unique(unlist(down.outlr, use.names = FALSE))
									)

					})


# method for OPPARList object setSubtypeProbes generic

#' @describeIn getSubtypeProbes A method for getSubtypeProbes with signature profileMatrix = \code{OPPARList}
#'   and sample.names = \code{nuchar}
#' @export
setMethod("getSubtypeProbes", signature(profileMatrix = "OPPARList", sample.names = "nuchar"),
					function(profileMatrix, sample.names){
						getSubtypeProbes(profileMatrix = profileMatrix@profileMatrix, sample.names = sample.names)
					})
